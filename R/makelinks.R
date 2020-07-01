makelinks <- function(strings, sep = '&') {
  alltumorgenes <- sarek.pathfindr::getEnvVariable("alltumorgenes")
  if (length(strings) > 0)
    for (i in 1:length(strings))
      if (strings[i] != '') {
        t = strsplit(strings[i], sep)[[1]]
        for (j in 1:length(t)) {
          if (t[j] %in% alltumorgenes) {
            t[j] = kableExtra::cell_spec(
              t[j],
              link = paste0(
                "https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=",
                t[j]
              ),
              format = 'html'
            )
          } else if (startsWith(t[j], 'rs')) {
            t[j] = kableExtra::cell_spec(
              t[j],
              link = paste0(
                "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",
                substr(t[j], 3, 100)
              ),
              format = 'html'
            )
          } else if (startsWith(t[j], 'COSM')) {
            t[j] = kableExtra::cell_spec(
              t[j],
              link = paste0(
                "https://cancer.sanger.ac.uk/cosmic/mutation/overview?id=",
                substr(t[j], 5, 100)
              ),
              format = 'html'
            )
          }
        }
        strings[i] = paste(t, collapse = ' ')
      }
  return(strings)
}
