haplotypeCallerOverlap <- function(haplotypecaller_N_file) {
  table <- NULL
  haplotypecaller_ids <- getEnvVariable('haplotypecaller_ids')
  
  if ( exists('haplotypecaller_N_file') && !is.null(haplotypecaller_ids) ) {
      names = names(haplotypecaller_ids)
      nsamples = length(names)
      table = matrix(
        '',
        ncol = nsamples,
        nrow = nsamples,
        dimnames = list(`Variants in` = names,
                        `Also in` = names)
      )
      for (i in 1:nsamples)
        for (j in 1:nsamples) {
          s = sum(haplotypecaller_ids[[names[i]]] %in% haplotypecaller_ids[[names[j]]])
          table[i, j] = paste(round(s / 1e6, 2), 'M')
        }
    }
  table
}