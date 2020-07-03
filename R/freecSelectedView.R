freecSelectedView <- function(freec_cnvs) {
  cat('Control Freec')
  if (!is.null(freec_cnvs) && nrow(freec_cnvs) > 0) {
      freec_cnvs = data.table(
        Sample = freec_cnvs$sample,
        Gene = freec_cnvs$cancer_genes,
        Mutation = freec_cnvs$effect,
        'Rank score' = freec_cnvs$rank_score,
        'Rank Terms' = freec_cnvs$rank_terms,
        Allele_ratio = ''
      )
      kable(freec_cnvs[Gene != ''], "html") %>%
        kable_styling(bootstrap_options = "striped", position = "left") %>%
        scroll_box(width = "900px", height = "300px")
    }
}