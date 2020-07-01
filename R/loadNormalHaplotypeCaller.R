loadNormalHaplotypeCaller <- function(haplotypecaller_N_file) {
  haplotypecaller_table=NULL
  haplotypecaller_ids=NULL
  
  cat("Reading large germline file, this can take a while ... \n")
  tic("Collecting variants from VCF")  
  vcf = VariantAnnotation::readVcf(file = haplotypecaller_N_file, genome = reference_genome)
  sample = strsplit(basename(haplotypecaller_N_file), '[.]')[[1]][1]
  haplotypecaller_ids[[sample]] = names(vcf)
  toc()
  
  # if (s>1) # if not the normal, keep only IDs already present (the normal is first)  <--- this is why somatic variants are not shown
  #   vcf=vcf[names(vcf) %in% haplotypecaller_table$ID]
  
  tic("pre-filtering by SWAF / TOPMED to avoid extreme number of variants")
  swaf = as.data.table(info(vcf)$SWAF)
  topmed = as.data.table(info(vcf)$TOPMED)
  values = data.table(swaf$value, topmed$value, 0) # where both are NA, 0 will be kept as AF
  values$max = apply(values, 1, max, na.rm = T)
  keep = unique(swaf$group[which(values$max < 0.01)])  # variant IXs with all alleles <1% in SWAF and TOPMED (swaf$group equals topmed$group)
  vcf = vcf[keep]
  toc()
  
  tic("Get germline mutations mutations from file")
  g = geno(vcf)
  inf = (info(vcf))
  rr = DelayedArray::rowRanges(vcf)
  
  mutations = as.data.table(ranges(rr))
  mutations$chr <- as.character(seqnames(rr))
  mutations$sample = sample
  mutations = mutations[, .(ID = names, sample, chr, start, end, width)]

  mutations$ref = as.character(ref(vcf))
  mutations$alt = as.data.table(alt(vcf))[, .(values = list(value)), by = group][, values]
  mutations$type = paste0(mutations$ref, '>', mutations$alt)
  mutations$type[nchar(mutations$ref) > nchar(mutations$alt)] = 'del'
  mutations$type[nchar(mutations$ref) < nchar(mutations$alt)] = 'ins'
  mutations$type[nchar(mutations$type) > 3] = 'other'
  mutations$type[mutations$type == 'T>G'] = 'A>C'
  mutations$type[mutations$type == 'T>C'] = 'A>G'
  mutations$type[mutations$type == 'T>A'] = 'A>T'
  mutations$type[mutations$type == 'G>T'] = 'C>A'
  mutations$type[mutations$type == 'G>C'] = 'C>G'
  mutations$type[mutations$type == 'G>A'] = 'C>T'

  mutations$AFreq = round(sapply(g$AD, "[[", 2) / unlist(g$DP), 2) # not a perfect estimate...
  ad = as.data.table(g$AD)
  colnames(ad) = 'AD'
  mutations = cbind(mutations, ad)
  # add counts from the new sources
  suppressWarnings({
    t = unlist(lapply(
      X = inf$CAF,
      FUN = max,
      na.rm = T
    ))
    t[is.infinite(t)] = 0
    mutations$CAF = t
    t = unlist(lapply(
      X = inf$SWAF,
      FUN = max,
      na.rm = T
    ))
    t[is.infinite(t)] = 0
    mutations$SWAF = t
    t = unlist(lapply(
      X = inf$TOPMED,
      FUN = max,
      na.rm = T
    ))
    t[is.infinite(t)] = 0
    mutations$TOPMED = t
  })
  toc()
  tic("Add INFO annotations and SweGen") 
  # "info"" annotations also need to be added
  mutations = cbind(mutations, as.data.table(info(vcf))[, .(CSQ)])

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
  
  # Add Swegen counts
  mutations$Swegen_count = snptable[names(vcf), value]
  mutations$Swegen_count[is.na(mutations$Swegen_count)] = 0

  # Also double check Swegen in case pos in reference but rsID or not properly matched in data
  pos = paste0(mutations$chr,
               ':',
               mutations$start,
               '_',
               mutations$ref,
               '/',
               mutations$alt)
  snps = snptable[c(pos)][!is.na(value)] # extract from reference[those present]
  if (nrow(snps) > 0)
    mutations$Swegen_count[match(snps$name, pos)] = snps$value # put
  # # Some have multiple alt alleles
  # ix=grep(',',mutations$alt)
  # alts=data.table(ix,alt=mutations$alt[ix])
  # alts$counts[ix]=unlist(lapply(alts$pos,function(x) return(max(snptable[x,value],na.rm=T))))

  # More thorough filtering
  if (any(mutations$Swegen_count >= 10))
    mutations = mutations[-which(Swegen_count >= 10)]
  toc()

  tic("Annotate by cosmic (for possibly retaining non PASS hotspots, which is not necessarily smart)")
  key = paste0(substr(mutations$chr, 4, 6),
               ':',
               mutations$start,
               '-',
               mutations$end)
  key = str_replace(key, 'X:', '23:')
  key = str_replace(key, 'Y:', '24:')

  counts = cbind(cosmic_coding[key, value], cosmic_noncoding[key, value], 0)
  max_ = apply(counts, 1, max, na.rm = T)
  mutations$Cosmic_count = max_
  toc()
  
  tic("Add Control Freec LOH.")
  mutations$rank_score = 0
  mutations$rank_terms = ''
  mutations$LOH = ''
  try({
    if (nrow(mutations) > 0 &
        !is.null(freec_loh))
      for (i in 1:nrow(freec_loh)) {
        ix = freec_loh$chr[i] == mutations$chr &
          freec_loh$start[i] < mutations$start &
          freec_loh$end[i] > mutations$end
        ix = ix[!is.na(ix)]
        if (sum(ix) > 0)
          mutations$LOH[ix] = 'Y'
      }
  }, silent = T)
  toc()

  vep_header = strsplit(info(header(vcf))['CSQ', ][, 3], 'Format: ')[[1]][2]
  vep_header = strsplit(vep_header, '\\|')[[1]]

  tic("Process VEP annotations")
  # VEP annotation is put in annotation_table
  annotation_table = matrix(
    data = NA,
    nrow = length(unlist(mutations$CSQ)),
    ncol = length(vep_header) + 1
  )
  colnames(annotation_table) = c('ID', vep_header)
  row = 1
  for (i in 1:nrow(mutations)) {
    # for each variant
    # for each VEP annotation:
    for (j in 1:length(mutations$CSQ[[i]])) {
      line = strsplit(mutations$CSQ[[i]][j], '\\|')[[1]]
      annotation_table[row, 1] = mutations$ID[i]
      annotation_table[row, 1 + (1:length(line))] = line
      row = row + 1
    }
  }
  annotation_table = as.data.table(annotation_table)

  # Annotation table then concatenated
  haplotypecaller_table = rbind(haplotypecaller_table,
                                merge(mutations, annotation_table, by = 'ID'))

  # done collecting from vcf files
  toc()
  haplotypecaller_table$cumstart=haplotypecaller_table$start
  haplotypecaller_table$cumend=haplotypecaller_table$end
  # # for each chr get cumulative pos
  for (i in 1:nrow(chrsz)) {
    ix=haplotypecaller_table$chr==chrsz$chr[i]
    haplotypecaller_table$cumstart[ix]=haplotypecaller_table$start[ix]+chrsz$starts[i]
    haplotypecaller_table$cumend[ix]=haplotypecaller_table$end[ix]+chrsz$starts[i]
  }
  setkey(haplotypecaller_table,'sample')

  tic("Due to size remove intron/intergenic variants")
  # CSQ is not needed after parsing/merging the annotations
  selection=haplotypecaller_table[!Consequence %in% c('intron_variant','intergenic_variant'),-c('CSQ')] 
  toc()
  if (nrow(selection)>0) {

    tic("Cosmic/local Tier2 :")
    ix=selection$SYMBOL %in% alltier2
    selection$rank_score[ix]=2
    selection$rank_terms[ix]='T2_gene'
    toc()
    
    tic("tier 1 priority:")
    ix=selection$SYMBOL %in% alltier1
    selection$rank_score[ix]=2
    selection$rank_terms[ix]='T1_gene'
    toc()
    
    tic("Add high impact")
    ix=selection$IMPACT=='HIGH'
    if (any(ix)) {
      selection$rank_score[ix]=selection$rank_score[ix]+2
      selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high_impact')
    }
    toc()
    
    tic("Add moderate impact")
    ix=selection$IMPACT=='MODERATE'
    if (any(ix)) {
      selection$rank_score[ix]=selection$rank_score[ix]+1
      selection$rank_terms[ix]=paste(selection$rank_terms[ix],'moderate_impact')
    }
    toc()
    
    tic("additional +1 if high impact and TSG")
    ix=selection$IMPACT=='HIGH' & selection$SYMBOL %in% alltsg
    if (any(ix)) {
      selection$rank_score[ix]=selection$rank_score[ix]+1
      selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high+TSG')
    }
    toc()
    
    tic("Add clinvar pathogenic")
    ix=grep('pathogenic',selection$CLIN_SIG)
    if (length(ix)>0) {
      selection$rank_score[ix]=selection$rank_score[ix]+2
      selection$rank_terms[ix]=paste(selection$rank_terms[ix],'clinvar')
    }
    toc()
    
    tic("Add Polyphen/SIFT damaging/deleterious")
    ix=union(grep('damaging',selection$PolyPhen),grep('deleterious',selection$SIFT))
    if (length(ix)>0) {
      selection$rank_score[ix]=selection$rank_score[ix]+1
      selection$rank_terms[ix]=paste(selection$rank_terms[ix],'polyphen/SIFT')
    }
    toc()
    
    tic("Add hotspots and near hotspots (within 2 residues of a hotspot)")
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
    toc()
    
    tic("Add cosmic counts")
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
    toc()
    
    tic("Special case for TERT promoter")
    ix=which(selection$SYMBOL=='TERT' & selection$Consequence %in% c('5_prime_UTR_variant','upstream_gene_variant'))
    if (length(ix)>0) {
      selection$rank_score[ix]=selection$rank_score[ix]+4
      selection$rank_terms[ix]=paste(selection$rank_terms[ix],'TERT_upstr')
    }
    toc()
    cat("Finding TF binding variants also takes a while ...\n")
    tic("Add TF binding variants near (100kb) cancer genes")
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
    toc()
    
    tic("Filter for canonical and CADD")
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
    toc()
  }

  firstcols=c('ID','sample','SYMBOL','rank_score','rank_terms','LOH','AFreq','Consequence','IMPACT','CADD_PHRED','SWAF','TOPMED','Swegen_count')
  cols=colnames(selection)
  setcolorder(x = selection,neworder = c(firstcols,cols[!cols %in% firstcols]))
  haplotypecaller_selected <- selection[order(cumstart,Allele)][order(rank_score,decreasing = T)]
  haplotypecaller_selected <- haplotypecaller_selected[rank_score>3]
  hc_filename <- paste0(csv_dir,'/',sample,'_haplotypecaller.csv')
  fwrite(haplotypecaller_selected,file=hc_filename)
  
  cat("Ranks written to ",hc_filename,"\n")  
  ix=haplotypecaller_selected$sample==strsplit(basename(haplotypecaller_N_file),'[.]')[[1]][1] # first file is the normal
  #tableWrapper(haplotypecaller_selected[ix][,-c('cumstart','cumend','DOMAINS')][rank_score>3])
  haplotypecaller_selected
}