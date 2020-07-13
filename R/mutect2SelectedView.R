mutect2SelectedView <- function(mutect2_selected) {
  cat('Mutect2')
  t = mutect2_selected[rank_score > getRankThreshold("mutect2_threshold")]
  if (nrow(t) > 0) {
    t = data.table(
      Sample = t$sample,
      Gene = t$SYMBOL,
      Mutation = paste(t$Protein_position, t$Amino_acids),
      'Rank score' = t$rank_score,
      'Rank Terms' = t$rank_terms,
      Allele_ratio = t$AFreq
    )
    kable(t[Gene != ''], "html") %>%
      kable_styling(bootstrap_options = "striped", position = "left") %>%
      scroll_box(width = "900px", height = "300px")
    
  }
}