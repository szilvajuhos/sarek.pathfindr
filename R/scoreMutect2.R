scoreMutect2 <- function(PFconfig) {
  cat(" ------------ Score GATK 4X Mutect2 function -------------\n")
  cat("Reference is: ",reference_genome,'\n')
  
  mutect2_result_file <- PFconfig$mutect2file
  cat("Mutect2 Somatic Variants file:  ", mutect2_result_file, '\n')
  loadMutect2(mutect2_result_file)
}