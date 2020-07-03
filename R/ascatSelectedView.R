ascatSelectedView <- function(ascat_cnvs) {
  cat('Ascat')
  if (nrow(ascat_cnvs) > 0) {
    ascat_cnvs = data.table(
      Sample = ascat_cnvs$sample,
      Gene = ascat_cnvs$cancer_genes,
      Mutation = paste(ascat_cnvs$nMajor, ascat_cnvs$nMinor),
      'Rank score' = ascat_cnvs$rank_score,
      'Rank Terms' = ascat_cnvs$rank_terms,
      Allele_ratio = ''
    )
    kable(ascat_cnvs[Gene != ''], "html") %>%
      kable_styling(bootstrap_options = "striped", position = "left") %>%
      scroll_box(width = "900px", height = "300px")
  }
}

