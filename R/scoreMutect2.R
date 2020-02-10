scoreMutect2 <- function() {
  cat(" ------------ Score GATK 4X Mutect2 function -------------\n")
  cat("Reference is: ",reference_genome,'\n')
  
  mutect2_result_file <- mutect2File()
  cat('Selected files for Mutect2 in',getwd(),':\n\n')
  cat("Mutect2 Somatic Variants file:  ", mutect2_result_file, '\n')
  loadMutect2(mutect2_result_file)
}