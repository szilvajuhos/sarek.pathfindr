getFusionTable <- function(ft_filename) {
  # get the snptable variable from the package-wide environment
  ft = get("fusion_table",envir = pfenv)
  # re-read the file only if it is not defined yet
  if(is.null(ft)) {
    tic(paste0("Reading fusions from file ",ft_filename))
    ft = data.table::fread(ft_filename, key = 'name')
    assign("fusion_table",ft,envir = pfenv)
    toc()
  }
  get("fusion_table",pfenv)
}