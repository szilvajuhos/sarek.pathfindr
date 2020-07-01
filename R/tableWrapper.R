tableWrapper <- function(table) {
  if ('Existing_variation' %in% colnames(table)) {
    table$Existing_variation = makelinks(table$Existing_variation)
  }
  if ('Gene_Name' %in% colnames(table)) {
    table$Gene_Name[nchar(table$Gene_Name) > 200] = 'many'
  }
  if ('SYMBOL' %in% colnames(table)) {
    table$SYMBOL = makelinks(table$SYMBOL)
  }
  if ('cancer_genes' %in% colnames(table)) {
    table$cancer_genes = makelinks(table$cancer_genes, sep = ' ')
  }
  if ('AD_TUMOR' %in% colnames(table)) {
    table$AD_TUMOR = unlist(lapply(table$AD_TUMOR, paste, collapse = ', '))
  }
  if ('AD_NORMAL' %in% colnames(table)) {
    table$AD_NORMAL = unlist(lapply(table$AD_NORMAL, paste, collapse = ', '))
  }
  if ('AD' %in% colnames(table)) {
    if (is.list(table$AD))
      table$AD = unlist(lapply(table$AD, paste, collapse = ', '))
  }
  
  samples = sort(unique(table$sample))
  
  if (T)
    if (nrow(table) > 0)
      htmlTable::htmlTable(
        table,
        col.rgroup = c("none", "#F7F7F7"),
        tspanner = paste('Rank score:', unique(table$rank_score)),
        n.tspanner = rev(table(table$rank_score))
      )
  # if (length(samples)>1) htmlTable(table,
  #           col.rgroup = c("none", "#F7F7F7"),
  #           n.rgroup='',
  #           tspanner=paste('Rank score:',unique(table$rank_score)),
  #           n.tspanner=table(table$rank_score))
}
