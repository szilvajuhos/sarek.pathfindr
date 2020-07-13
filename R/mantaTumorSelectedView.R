mantaTumorSelectedView <- function(manta_tumor_selected) {
  cat('Manta somatic')
  try({
    t = manta_tumor_selected[rank_score > getRankThreshold("manta_threshold")]
    t$Gene_Name[nchar(t$Gene_Name) > 100] = 'many'
    if (nrow(t) > 0) {
      t = unique(
        data.table(
          Sample = t$sample,
          Gene = t$Gene_Name,
          Mutation = t$SVTYPE,
          'Rank score' = t$rank_score,
          'Rank Terms' = t$rank_terms,
          Allele_ratio = t$AFreq
        )
      )
      kable(t[Gene != ''], "html") %>%
        kable_styling(bootstrap_options = "striped", position = "left") %>%
        scroll_box(width = "900px", height = "300px")
    }
  }, silent = T)
}