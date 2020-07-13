getRankThreshold <- function(caller_str) {
  threshold <- NULL
  PFconfig <- getEnvVariable('PFconfig')
  if(caller_str %in% names(PFconfig)) {
    threshold <- PFconfig[caller_str]
  } else if("default_threshold" %in% names(PFconfig)) {
    threshold <- PFconfig["default_threshold"]
  }
  
  threshold
}