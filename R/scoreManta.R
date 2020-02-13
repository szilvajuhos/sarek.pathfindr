scoreManta <- function() {
  cat(" ------------ Score MANTA function -------------\n")
  cat("Reference is: ",reference_genome,'\n')
  
  manta_result_files <- mantaFiles()
  cat('Selected files for manta in',getwd(),':\n\n')
  cat("manta Somatic Structural Variants file:  ", manta_result_files$manta_tumor_file, '\n')
  cat("manta Normal Structural Variants file: ", manta_result_files$manta_normal_file, '\n')
  #cat("manta SweGen file: ", manta_result_files["swegen_manta_all"], '\n')

  tic("Calculating Manta germline scores")
  loadGermlineManta(manta_result_files)
  toc()  
  
  tic("Calculating Manta somatic scores")
  #loadSomaticManta(manta_result_files)
  toc()
}
