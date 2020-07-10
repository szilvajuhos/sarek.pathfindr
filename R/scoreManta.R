scoreManta <- function() {
  cat(" ------------ Score MANTA function -------------\n")
  
  manta_result_files <- mantaFiles()
  cat('Selected files for manta: \n')
  cat("manta Somatic Structural Variants file:  ", manta_result_files$manta_tumor_file, '\n')
  cat("manta Normal Structural Variants file: ", manta_result_files$manta_normal_file, '\n')
  
  cat("************** Calculating Manta germline scores **************")
  tic("Manta Germline")
  loadGermlineManta(manta_result_files)
  toc()  
  
  cat("************** Calculating Manta somatic scores **************")
  tic("Manta Somatic")
  loadSomaticManta(manta_result_files)
  toc()
}
