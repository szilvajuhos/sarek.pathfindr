scoreHaplotypeCaller <- function(PFconfig) {
  hc_selected <- loadHaplotypeCaller( haplotypeCallerFiles(PFconfig) )
  hc_selected
}