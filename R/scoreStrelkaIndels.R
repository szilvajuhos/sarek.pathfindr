scoreStrelkaIndels <- function(strelka_indel_file,sample) {
  allpass=NULL
  table_indels=NULL
  reference_genome<-getEnvVariable('reference_genome')
  vcfs=list()
  if (length(strelka_indel_file)>0) {
    # read indels
    vcfs[[strelka_indel_file]]=VariantAnnotation::readVcf(file = strelka_indel_file,genome=reference_genome)
    pass=DelayedArray::rowRanges(vcfs[[strelka_indel_file]])$FILTER=='PASS'
    allpass=c(allpass,names(vcfs[[strelka_indel_file]])[pass])
  }
  vcf=vcfs[[strelka_indel_file]]
  vcf=vcf[names(vcf) %in% allpass] # only those with PASS in either sample
  if (!'TUMOR' %in% colnames(vcf)) 
    colnames(vcf)=c('NORMAL','TUMOR') # assumption that normal comes first
  # Collect sample-unspecific data for the [PASS in any] IDs into a table
  g=geno(vcf)
  rr=DelayedArray::rowRanges(vcf)
  inf=info(vcf)
  csq=NULL
  if (length(vcf)>0) {
    table_indels=data.table(ID=as.character(names(vcf)),
                            sample=sample,
                            chr=as.character(seqnames(rr)),
                            start=start(rr),
                            end=end(rr),
                            ref=as.data.table(rr$REF)$x,
                            alt=as.data.table(rr$ALT)$value,
                            AFreq=NA,AD=NA,DP=as.data.table(g$DP)$TUMOR,
                            AD_normal=NA,DP_normal=as.data.table(g$DP)$NORMAL)
    refcounts=as.data.table(data.frame(g$TAR))
    altcounts=as.data.table(data.frame(g$TIR))
    table_indels$AFreq=round(altcounts$TUMOR.1 / (altcounts$TUMOR.1 + refcounts$TUMOR.1),2)
    table_indels$AD_normal=altcounts$NORMAL.1
    table_indels$AD=altcounts$TUMOR.1
    table_indels$reads='≥5'
    table_indels$reads[table_indels$AD<5]='<5'
    # CAF,SWAF,TOPMED
    suppressWarnings({
      t=unlist(lapply(X=inf$CAF,FUN=max,na.rm=T))
      t[is.infinite(t)]=0
      table_indels$CAF=t
      t=unlist(lapply(X=inf$SWAF,FUN=max,na.rm=T))
      t[is.infinite(t)]=0
      table_indels$SWAF=t
      t=unlist(lapply(X=inf$TOPMED,FUN=max,na.rm=T))
      t[is.infinite(t)]=0
      table_indels$TOPMED=t
    })
    
    # indel annotations added to snv annotations
    csq=c(csq,as.data.table(info(vcf))[,CSQ])
  }
  rm(vcfs)
  c("table_indels"=list(table_indels),"indel_csq"=list(csq),"indel_vcf"=vcf)
}
