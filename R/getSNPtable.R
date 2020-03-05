getSNPtable <- function(st_filename) {
  # get the snptable variable from the package-wide environment
  st = get("snptable",envir = pfenv)
  # re-read the file only if it is not defined yet
  if(is.null(st)) {
    tic(paste0("Reading population-specific SNV frequencies from file ",st_filename))
    st = data.table::fread(st_filename, key = 'name')
    assign("snptable",st,envir = pfenv)
    toc()
  }
  get("snptable",pfenv)
}