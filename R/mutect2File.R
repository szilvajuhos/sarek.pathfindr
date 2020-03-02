mutect2File <- function(PFconfig) {
  result_files <- dir(path=PFconfig, recursive = F,full.names = T)
  mutect2_result_file <- grep(pattern = "mutect2.*_vs_.*vep.ann.vcf$",result_files,value = T)
}