scoreManta <- function(PFconfig) {
  cat(" ------------ Score MANTA function -------------\n")
  cat("Reference is: ",reference_genome,'\n')
  
  manta_result_files <- mantaFiles(PFconfig)
  cat('Selected files for manta: \n')
  cat("manta Somatic Structural Variants file:  ", manta_result_files$manta_tumor_file, '\n')
  cat("manta Normal Structural Variants file: ", manta_result_files$manta_normal_file, '\n')
  
  tic("************** Calculating Manta germline scores **************")
  loadGermlineManta(manta_result_files)
  toc()  
  
  tic("************** Calculating Manta somatic scores **************")
  loadSomaticManta(manta_result_files)
  toc()
}
