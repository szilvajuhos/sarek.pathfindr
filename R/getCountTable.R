getCountTable <- function(t_filename,table_name) {
  # get the variable from the package-wide environment
  t = get(table_name, envir = pfenv)
  # re-read the file only if it is not defined yet
  if(is.null(t)) {
    tic(paste0("Reading region counts for ",table_name," from file ",t_filename))
    t = data.table::fread(t_filename, key = 'name')
    assign(table_name,t,envir = pfenv)
    toc()
  }
  get(table_name,pfenv)
}