saveToBed <-
  function(sample,
           manta_tumor_table,
           manta_normal_table,
           mutect2_table,
           strelka_table,
           haplotypecaller_table) {
    x = unique(
      rbind(
        manta_tumor_table[, .(chr, start, end = start)],
        manta_tumor_table[, .(chr, start = end, end)],
        manta_tumor_table[!is.na(altpos), .(chr = altchr, start = altpos, end =
                                              altpos)],
        manta_normal_table[, .(chr, start, end = start)],
        manta_normal_table[, .(chr, start = end, end)],
        manta_normal_table[!is.na(altpos), .(chr = altchr, start = altpos, end =
                                               altpos)],
        mutect2_table[, .(chr, start, end)],
        strelka_table[, .(chr, start, end)],
        haplotypecaller_table[, .(chr, start, end)]
      )
    )
    bed_name <- paste0(csv_dir, '/', sample, '_PF.bed')
    write.table(
      x,
      file = bed_name,
      sep = '\t',
      row.names = F,
      col.names = F,
      quote = F
    )
    cat("Regions written to ",bed_name)
  }