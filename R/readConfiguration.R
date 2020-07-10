readConfiguration <- function(config_file) {
  cat("Reading configuration and auxiliary datafiles from",config_file,"\n")
  PFconfig <- config::get(file = config_file)
  reference_genome <- PFconfig$reference
  assign(x = "PFconfig", PFconfig, envir = pfenv)
  assign(x = "reference_genome", reference_genome, envir = pfenv)    
  PFconfig
}