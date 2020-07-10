mantaFiles <- function() {
  PFconfig<-getEnvVariable('PFconfig')
  # OK, this part has to be at a higher level at init
  result_files <- dir(path = PFconfig$manta_directory, recursive = F,full.names = T)
  
  # Manta structural variant files
  manta_tumor_file<-grep(pattern = ".*/Manta_.*vs.*somaticSV.vcf.*.vcf$",result_files,value = T)
  manta_normal_file<-grep(pattern = ".*/Manta_.*vs.*diploidSV.vcf.*.vcf$",result_files,value = T)[1]
  tic("Read SweGen SV counts")
  # TODO: get rid of hardcoded allele-frequency database
  swegen_manta_all<-list(data.table::fread('~/reports/reference_data/swegen_sv_counts.csv',key='name'))
  toc()
  
  # return with a file hash
  c("manta_tumor_file" = manta_tumor_file,
    "manta_normal_file" = manta_normal_file,
    "swegen_manta_all" = swegen_manta_all
    )
}