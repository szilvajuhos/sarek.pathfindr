scoreStrelkaSNVs <- function(strelka_snv_file,sample) {
  allpass=NULL
  table_snvs=NULL
  csq=NULL
  vcf=NULL
  vcfs=list()
  tic("Reading SNVs VCF")
  reference_genome<-getEnvVariable('reference_genome')
  if (length(strelka_snv_file)>0) for (s in 1:length(strelka_snv_file)) {
    # read snvs
    vcfs[[strelka_snv_file]]=VariantAnnotation::readVcf(file = strelka_snv_file,genome=reference_genome)
    pass=DelayedArray::rowRanges(vcfs[[strelka_snv_file]])$FILTER=='PASS'
    allpass=c(allpass,names(vcfs[[strelka_snv_file]])[pass])
  }
  toc()
  tic("Processing SNV annotations")
  if (length(strelka_snv_file)>0 ) {
    # first snvs
    vcf=vcfs[[strelka_snv_file]]
             vcf=vcf[names(vcf) %in% allpass] # only those with PASS in either sample
             if (!'TUMOR' %in% colnames(vcf)) 
               colnames(vcf)=c('NORMAL','TUMOR') # assumption that normal comes first
             # Collect sample-unspecific data for the [PASS in any] IDs into a table
             g=geno(vcf)
             rr=DelayedArray::rowRanges(vcf)
             inf=info(vcf)
             if (length(vcf)>0) {
               table_snvs=data.table(ID=as.character(names(vcf)),
                                     sample=sample,
                                     chr=as.character(seqnames(rr)),
                                     start=start(rr),
                                     end=end(rr),
                                     ref=as.data.table(rr$REF)$x,
                                     alt=as.data.table(rr$ALT)$value,
                                     AFreq=NA,AD=NA,DP=as.data.table(g$DP)$TUMOR,
                                     AD_normal=NA,DP_normal=as.data.table(g$DP)$NORMAL)
               # Collect variant allele depths for SNVs:
               for (this_ref in unique(table_snvs$ref)) {
                 for (this_alt in unique(table_snvs[ref==this_ref,alt])) {
                   which_ones=table_snvs$ref==this_ref & table_snvs$alt==this_alt
                   refcounts=as.data.table(data.frame(g[[paste0(this_ref,'U')]]))
                   altcounts=as.data.table(data.frame(g[[paste0(this_alt,'U')]]))
                   table_snvs$AFreq[which_ones]=round(altcounts[which_ones,TUMOR.1] /
                                                        (altcounts[which_ones,TUMOR.1]+refcounts[which_ones,TUMOR.1]),2)
                   table_snvs$AD[which_ones]=altcounts[which_ones,TUMOR.1]
                   table_snvs$AD_normal[which_ones]=altcounts[which_ones,NORMAL.1]
                 }
               }
               table_snvs$reads='≥5'
               table_snvs$reads[table_snvs$AD<5]='<5'
               # CAF,SWAF,TOPMED
               suppressWarnings({
                 t=unlist(lapply(X=inf$CAF,FUN=max,na.rm=T))
                 t[is.infinite(t)]=0
                 table_snvs$CAF=t
                 t=unlist(lapply(X=inf$SWAF,FUN=max,na.rm=T))
                 t[is.infinite(t)]=0
                 table_snvs$SWAF=t
                 t=unlist(lapply(X=inf$TOPMED,FUN=max,na.rm=T))
                 t[is.infinite(t)]=0
                 table_snvs$TOPMED=t
               } )
               
               # The snv annotations will be used later:
               csq=as.data.table(info(vcf))[,CSQ]
             }
  
  }
  rm(vcfs)
  toc()
  # return with the collected SNVs and VEP annotations
  c("table_snvs"=list(table_snvs),"snv_csq"=list(csq),"snv_vcf"=vcf)
}
