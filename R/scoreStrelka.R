scoreStrelka <- function() {
  cat(" ------------ Score Strelka function -------------\n")
  cat("Reference is: ",reference_genome,'\n')
  
  strelka_result_files <- strelkaFiles()
  cat('Selected files for Strelka in',getwd(),':\n\n')
  cat("Strelka Somatic SNVs file:  ", strelka_result_files["strelka_snv_file"], '\n')
  cat("Strelka Somatic indels file:  ", strelka_result_files["strelka_indel_file"], '\n')
  indels <- scoreStrelkaIndels(strelka_result_files["strelka_indel_file"])
  snvs <- scoreStrelkaSNVs(strelka_result_files["strelka_snv_file"])
}