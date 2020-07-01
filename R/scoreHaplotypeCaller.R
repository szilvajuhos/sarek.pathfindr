scoreHaplotypeCaller <- function(PFconfig) {
  hc_file <- haplotypeCallerFiles(PFconfig)
  cat("HaplotypeCaller NORMAL file: ",hc_file["hc_normal"],"\n")
  cat("HaplotypeCaller TUMOR file: ", hc_file["hc_tumor"],"\n")
  
  tic("Ranking NORMAL calls")
  hc_selected <- loadNormalHaplotypeCaller(hc_file["hc_normal"])
  toc()
  
  # tic("Processing TUMOR calls")
  # loadTumorHaplotypeCaller(hc_file["hc_tumor"])
  # toc()
  
  hc_selected
}