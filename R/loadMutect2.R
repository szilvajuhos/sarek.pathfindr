loadMutect2 <- function(mutect2_file) {
  PFconfig <- getEnvVariable('PFconfig')
  reference_genome<-getEnvVariable('reference_genome')
  
  cat("Processing Mutect2 SNV calls\n")
  # used to store final results
  mutect2_selected <- NULL
  # local variable
  allpass <- NULL
  # all the aux data used
  snptable <- getSNPtable(PFconfig$snptable)
  cosmic_coding <- getCountTable(PFconfig$coding_table,"coding_table")
  cosmic_noncoding <- getCountTable(PFconfig$noncoding_table,"noncoding_table")
  cosmic_fusions <- getCountTable(PFconfig$fusions_table,"fusions_table")
  hotspots_snv <- getSNVhotspots(PFconfig$hotspots_snv)
  
  # these two below have to be together always (TODO: find a safer way to get near_hotspots)
  hotspots_inframe <- getInframeHotspots(PFconfig$hotspots_inframe,hotspots_snv)
  near_hotspots <- get("near_hotspots",pfenv)

  # this function creates more tables, get them via get(tablename,pfenv)
  getTumorGenes(PFconfig$tumorgenes, PFconfig$local_tumorgenes)
  tumorgenes <- get("tumorgenes",pfenv)
  alltumorgenes <- get("alltumorgenes",pfenv)
  # get TSGs
  alltsg <- get("alltsg",pfenv)
  # get tier data
  alltier1 <- get("alltier1",pfenv)
  alltier2 <- get("alltier2",pfenv)

  if (exists('mutect2_file')) if (length(mutect2_file)>0) {
    for (s in 1:length(mutect2_file)) {
      vcf=readVcf(file = mutect2_file[s],genome = reference_genome)
      pass=rowRanges(vcf)$FILTER=='PASS'
      allpass=c(allpass,names(vcf)[pass])
    }
    cat("Collecting variants...\n")
    mutect2_table=NULL
    for (s in 1:length(mutect2_file)) {
      sample=strsplit(basename(mutect2_file[s]),'[.]')[[1]][1]
      
      vcf=readVcf(file = mutect2_file[s],genome = reference_genome)
      vcf=vcf[names(vcf) %in% allpass]
      if (length(vcf)>0) {
        rr=DelayedArray::rowRanges(vcf)
        g=geno(vcf)
        inf=info(vcf)
        
        # manipulate into a data frame with relevant data
        mutations=as.data.table(ranges(rr))
        mutations$chr <- as.character(seqnames(rr))
        mutations$sample=sample
        mutations=mutations[,.(ID=names,sample,chr,start,end,width)]
        
        mutations$ref=as.character(ref(vcf))
        mutations$alt=as.data.table(alt(vcf))[, .(values = list(value)), by = group][,values]
        mutations$type=paste0(mutations$ref,'>',mutations$alt)
        mutations$type[nchar(mutations$ref)>nchar(mutations$alt)]='del'
        mutations$type[nchar(mutations$ref)<nchar(mutations$alt)]='ins'
        mutations$type[nchar(mutations$type)>3]='other'
        mutations$type[mutations$type=='T>G']='A>C'
        mutations$type[mutations$type=='T>C']='A>G'
        mutations$type[mutations$type=='T>A']='A>T'
        mutations$type[mutations$type=='G>T']='C>A'
        mutations$type[mutations$type=='G>C']='C>G'
        mutations$type[mutations$type=='G>A']='C>T'
        
        # Headers are sometimes sample names instead of TUMOR-NORMAL
        headers=colnames(g$AD)
        if (!'TUMOR' %in% headers) {
          ix=grep('-02$',headers)
          if (length(ix)==0) ix=grep('[TR]',headers)
          headers[ix]='TUMOR'
          headers[-ix]='NORMAL'
          colnames(g$AD)=headers
          colnames(g$AF)=headers # <-- assumes the same header order in AF
        }
        mutations$AFreq=round(sapply(g$AF[,'TUMOR'], "[", 1),2)
        ad=as.data.table(g$AD)
        colnames(ad)=paste0('AD_',colnames(ad))
        mutations=cbind(mutations,ad)
        # mutations$DP=sapply(g$AD[,'TUMOR'], sum)
        # mutations$AD=sapply(g$AD[,'TUMOR'], "[[", 2)
        # mutations$DPn=sapply(g$AD[,'NORMAL'], sum)
        # mutations$ADn=sapply(g$AD[,'NORMAL'], "[[", 2)
        mutations$reads=''
        temp=sapply(g$AD[,'TUMOR'], "[[", 2)
        mutations$reads[temp<5]='<5'
        mutations$reads[temp>=5]='≥5'
        
        # add counts from the new sources
        suppressWarnings( {
          t=unlist(lapply(X=inf$CAF,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          mutations$CAF=t
          t=unlist(lapply(X=inf$SWAF,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          mutations$SWAF=t
          t=unlist(lapply(X=inf$TOPMED,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          mutations$TOPMED=t
        } )
        
        # Add Swegen counts
        mutations$Swegen_count=snptable[names(vcf),value]
        mutations$Swegen_count[is.na(mutations$Swegen_count)]=0
        # # Some have multiple rsIDs (diminishingly few)
        # ix=grep(';',mutations$names)
        # ids=data.table(ix,id=strsplit(mutations$names[ix],';'))
        # ids$counts=unlist(lapply(ids$id,function(x) return(max(snptable[x,value],na.rm=T))))
        # ids$counts[is.infinite(ids$counts)]=0
        # mutations$Swegen_count[ix]=ids$counts
        # Also double check Swegen in case pos in reference but rsID or not properly matched in data
        pos=paste0(mutations$chr,':',mutations$start,'_',mutations$ref,'/',mutations$alt)
        snps=snptable[c(pos)][!is.na(value)] # extract from reference[those present]
        if (nrow(snps)>0) mutations$Swegen_count[match(snps$name,pos)]=snps$value # put
        # # Some have multiple alt alleles (diminishingly few)
        # ix=grep(',',mutations$alt)
        # for (i in ix) {
        #   alts=mutations$alt[[i]]
        #   snps=snptable[paste0('chr',mutations$chr[i],':',mutations$start[i],'_',mutations$ref[i],'/',alts)]
        #   mutations$Swegen_count[i]=max(c(0,snps$value),na.rm = T)
        # }
        
        cat("Cosmic annotations\n")
        ## Annotate by cosmic (for possibly retaining non PASS hotspots, which is not necessary smart)
        key=paste0(substr(mutations$chr,4,6),':',mutations$start,'-',mutations$end)
        key=str_replace(key,'X:','23:')
        key=str_replace(key,'Y:','24:')
        
        counts=cbind(cosmic_coding[key,value],cosmic_noncoding[key,value],0)
        max_=apply(counts,1,max,na.rm=T)
        mutations$Cosmic_count=max_
        
        mutations$rank_score=0
        mutations$rank_terms=''
        mutations$LOH=''
        
        # Add Control Freec LOH.
        cat("Adding LOH from Control-FREEC")
        try( {
          if (nrow(mutations)>0 & !is.null(freec_loh)) for (i in 1:nrow(mutations)) {
            ix=freec_loh$chr==mutations$chr[i] &
              freec_loh$start<mutations$start[i] &
              freec_loh$end>mutations$end[i]
            ix=ix[!is.na(ix)]
            if (sum(ix) >0) mutations$LOH[i]='Y'
          }},silent=T)
        
        # "info"" annotations also need to be added
        mutations=cbind(mutations,as.data.table(info(vcf))[,.(CSQ)])
        
        # keep only PASS
        # mutations=mutations[filter=='PASS',]
        
        # get snpEff headers
        # snpeff_header=strsplit(info(header(vcf))['ANN',][,3],'ions: \'')[[1]][2]
        # snpeff_header=strsplit(snpeff_header,' \\| ')[[1]]
        ## get VEP headers
        vep_header=strsplit(info(header(vcf))['CSQ',][,3],'Format: ')[[1]][2]
        vep_header=strsplit(vep_header,'\\|')[[1]]
        
        cat("VEP into annotation table\n")
        # VEP annotation is put in annotation_table
        annotation_table=matrix(data = NA,nrow = length(unlist(mutations$CSQ)),ncol = length(vep_header)+1)
        colnames(annotation_table)=c('ID',vep_header)
        row=1
        for (i in 1:nrow(mutations)) { # for each variant
          # for each VEP annotation:
          for (j in 1:length(mutations$CSQ[[i]])) {
            line=strsplit(mutations$CSQ[[i]][j],'\\|')[[1]]
            annotation_table[row,1]=mutations$ID[i]
            annotation_table[row,1+(1:length(line))]=line
            row=row+1
          }
        }
        annotation_table=as.data.table(annotation_table)
        
        # Annotation table then concatenated
        mutect2_table=rbind(mutect2_table,merge(mutations,annotation_table,by='ID'))
      }
    } # done collecting from vcf files
    
    
    mutect2_table$cumstart=NA
    mutect2_table$cumend=NA
    cat("For each chromosome get cumulative values\n")
    for (i in 1:nrow(chrsz)) {
      ix=mutect2_table$chr==chrsz$chr[i]
      mutect2_table$cumstart[ix]=mutect2_table$start[ix]+chrsz$starts[i]
      mutect2_table$cumend[ix]=mutect2_table$end[ix]+chrsz$starts[i]
    }
    setkey(mutect2_table,'sample')
    selection=mutect2_table[,-c('CSQ')] # not needed after parsing/merging the annotations

    if (nrow(selection)>0) {
      
      # Cosmic/local Tier2 :
      ix=selection$SYMBOL %in% alltier2
      selection$rank_score[ix]=2
      selection$rank_terms[ix]='T2_gene'
      # tier 1 priority:
      ix=selection$SYMBOL %in% alltier1
      selection$rank_score[ix]=2
      selection$rank_terms[ix]='T1_gene'
      
      cat("Add high impact\n")
      ix=selection$IMPACT=='HIGH'
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high_impact')
      }
      
      cat("Add moderate impact\n")
      ix=selection$IMPACT=='MODERATE'
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'moderate_impact')
      }
      
      cat("Additional +1 if high impact and TSG\n")
      ix=selection$IMPACT=='HIGH' & selection$SYMBOL %in% alltsg
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high+TSG')
      }
      
      cat("Add clinvar pathogenic\n")
      ix=grep('pathogenic',selection$CLIN_SIG)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'clinvar')
      }
      
      # Add Polyphen/SIFT damaging/deleterious
      ix=union(grep('damaging',selection$PolyPhen),grep('deleterious',selection$SIFT))
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'polyphen/SIFT')
      }
      
      # Add hotspots and near hotspots (within 2 residues of a hotspot)
      key=paste(selection$SYMBOL,
                str_replace(string = selection$Protein_position,
                            pattern = '/.*',replacement = ''))
      ix <-
        key %in% paste(hotspots_snv[,Hugo_Symbol],hotspots_snv[,Amino_Acid_Position]) |
        key %in% paste(hotspots_inframe[,Hugo_Symbol],hotspots_inframe[,Amino_Acid_Position])  ## Warning: Exact match used with inframes
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'hotspot')
      }
      ix = !ix & key %in% near_hotspots # the near_hotspot excludes those that were a hotspot
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'near_hotspot')
      }
      
      # Add cosmic counts
      ix=which(selection$Cosmic_count>50)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'cosmic_>50')
      }
      ix=which(selection$Cosmic_count>5 & selection$Cosmic_count<=50)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'cosmic_>5')
      }
      
      # Special case for TERT promoter
      ix=which(selection$SYMBOL=='TERT' & selection$Consequence %in% c('5_prime_UTR_variant','upstream_gene_variant'))
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+4
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'TERT_upstr')
      }
      
      # Add TF binding variants near (100kb) cancer genes
      ix=grep('TF',selection$Consequence)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'TFBS')
        for (i in 1:length(ix)) {
          near=which(selection$cumend[ix[i]] > (tumorgenes$cumstart-100e3) & selection$cumstart[ix[i]] < (tumorgenes$cumend + 100e3) &
                       !selection$SYMBOL[ix[i]] %in% alltumorgenes)
          if (length(near)>0) {
            genes=paste(tumorgenes$`Gene Symbol`[near],collapse = ',')
            selection$rank_score[ix[i]]=selection$rank_score[ix[i]]+2
            selection$rank_terms[ix[i]]=paste(selection$rank_terms[ix[i]],paste0('near_',genes))
          }
        }
      }
      
      ix=which(selection$CANONICAL!='YES')
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]-0
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'not_canonical')
      }
      
      ix=which(as.numeric(selection$CADD_PHRED) > 30)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'cadd_>30')
      }
    }
    
    firstcols=c('ID','sample','SYMBOL','rank_score','rank_terms','LOH','AFreq','Consequence','IMPACT','CADD_PHRED','SWAF','TOPMED','Swegen_count')
    cols=colnames(selection)
    setcolorder(x = selection,neworder = c(firstcols,cols[!cols %in% firstcols]))
    selection <- selection[order(cumstart,Allele)][order(rank_score,decreasing = T)]
    
    mutect2_selected <- selection[rank_score>3]
    m2fname <- paste0(csv_dir,'/',sample,'_mutect2_tumor.csv')
    fwrite(mutect2_selected,file=m2fname)
    cat("Ranks written to ",m2fname,"\n")
    
#    if (write_tables) 
#      fwrite(mutect2_selected[rank_score>3],file=paste0(csv_dir,'/',sampleData$name,'_mutect2_tumor.csv'))
#    tableWrapper(mutect2_selected[,-c('cumstart','cumend','DOMAINS')][rank_score>3])
  } else {
    cat("No variants considered for scoring")
  }
  mutect2_selected
}