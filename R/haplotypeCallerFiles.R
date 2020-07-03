haplotypeCallerFiles <- function(PFconfig) {
  result_files <- dir(path = PFconfig$haplotyper_directory, recursive = F,full.names = T)
  tic("reading HaplotypeCaller variant files")
  haplotypecaller_T_file <- grep(pattern = ".*haplotypecaller.*[TR2][.].*AF.*vep.ann.vcf$",result_files,value = T)
  haplotypecaller_N_file <- grep(pattern = ".*haplotypecaller.*[NB1][.].*AF.*vep.ann.vcf$",result_files,value = T)
  toc()
  c(
    hc_normal=haplotypecaller_N_file,
    hc_tumor=haplotypecaller_T_file
    )
}