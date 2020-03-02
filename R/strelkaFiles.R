strelkaFiles <- function(PFconfig) {
  # OK, this part has to be at a higher level at init
  result_files <- dir(path = PFconfig$strelka_directory, recursive = F,full.names = T)

  strelka_snv_file   <- grep(pattern = ".*/Strelka_.*_somatic_snvs.*vep.ann.vcf$",result_files,value = T)
  strelka_indel_file <- grep(pattern = ".*/Strelka_.*_somatic_indels.*vep.ann.vcf$",result_files,value = T)

  # return with a file hash
  c("strelka_snv_file" = strelka_snv_file,
    "strelka_indel_file" = strelka_indel_file
  )
}