scoreMutect2 <- function() {
  PFconfig<-getEnvVariable('PFconfig')
  mutect2_result_file <- PFconfig$mutect2file
  cat("Mutect2 Somatic Variants file:  ", mutect2_result_file, '\n')
  loadMutect2(mutect2_result_file)
}