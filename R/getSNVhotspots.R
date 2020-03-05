getSNVhotspots <- function(shs_filename) {
  # get the hotspots_snv variable from the package-wide environment
  shs = get("hotspots_snv",envir = pfenv)
  # re-read the file only if it is not defined yet
  if(is.null(shs)) {
    tic(paste0("Reading SNV hotspots from file ",shs_filename))
    shs <- unique(
      data.table::fread(shs_filename)[,.(Hugo_Symbol,Amino_Acid_Position)])[-grep('splice',Amino_Acid_Position)]
    shs$pos <- as.numeric(shs$Amino_Acid_Position)
    assign("hotspots_snv",shs,envir = pfenv)
    toc()
  }
  get("hotspots_snv",pfenv)
}