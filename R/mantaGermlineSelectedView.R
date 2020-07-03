mantaGermlineSelectedView <- function(manta_normal_selected) {
  cat('Manta germline')
  try({
    manta_normal_selected$Gene_Name[nchar(manta_normal_selected$Gene_Name) >
                                      100] = 'many'
    if (nrow(manta_normal_selected) > 0) {
      manta_normal_selected = unique(
        data.table(
          Sample = manta_normal_selected$sample,
          Gene = manta_normal_selected$Gene_Name,
          Mutation = manta_normal_selected$SVTYPE,
          'Rank score' = manta_normal_selected$rank_score,
          'Rank Terms' = manta_normal_selected$rank_terms,
          Allele_ratio = manta_normal_selected$AFreq
        )
      )
      kable(manta_normal_selected[Gene != ''], "html") %>%
        kable_styling(bootstrap_options = "striped", position = "left") %>%
        scroll_box(width = "900px", height = "300px")
    }
  }, silent = T)
}