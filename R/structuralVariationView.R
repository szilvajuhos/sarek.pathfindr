structuralVariationView <- function(manta_vcfs) {
  # manta_tumor_selected
  # manta_tumor_table
  # cf and ascat tratio
  # 
  
  g <- NULL # the plot we are returning with
  manta_tumor_file <- manta_vcfs['manta_tumor_file'][[1]]
  manta_normal_file <- manta_vcfs['manta_normal_file'][[1]]
  cat(manta_tumor_file, manta_normal_file, sep = '\n')
  #browser()
  
  # to get Control-FREEC Tratio data run load_freec()
  cf_tratio <- NULL
  try( cf_tratio <- get('cf_tratio',pfenv) )
  # ditto to get ASCAT tratio from env, run load_ascat()
  # we actually do not need ascat_tratio if CF is already there
  ascat_tratio <- NULL
  try(ascat_tratio <- get('ascat_tratio',pfenv),silent=T)
  
  manta_tumor_selected <- get('manta_tumor_selected',pfenv)
  sv_table = unique(manta_tumor_selected[, .(sample, chr, altchr, cumstart, altcumpos, IMPRECISE, SVTYPE, plot)])
  setkey(sv_table, 'sample')
  sv_samples = sort(unique(sv_table$sample))
  if (exists('cf_tratio')) {
    cnv_table = cf_tratio[seq(1, nrow(cf_tratio), 10), .(sample, x = cumstart, y =
                                                     log2(Ratio))]
    cnv_samples = sort(unique(cnv_table$sample))
  } else if (exists('ascat_tratio')) {
    cnv_table = ascat_tratio[seq(1, nrow(ascat_tratio), 100), .(sample, x =
                                                                  cumstart, y = log2(smoothed))]
    cnv_samples = sort(unique(cnv_table$sample))
  }
  
  for (s in 1:length(sv_samples)) {
    cat('Manta:', sv_samples[s])
    
    par(mai = c(0, 4, 0, 2))
    g = ggplot() + ylab('Structural variants') + expand_limits(x = c(0, 3210e6)) +
      ylim(-1.5, 1.5) +
      xlab(paste0('SVs: ', sv_samples[s], ',  CNVs: ', cnv_samples[s])) +
      ylab(sv_samples[s]) +
      scale_color_manual(
        values = c(
          "BND" = "#E41A1C",
          "DEL" = "#377EB8",
          "INV" = "#4DAF4A",
          "INS" = "#984EA3",
          "DUP" = "#FF7F00"
        )
      )#+
    #scale_linetype_manual(values=c('TRUE'=2,'FALSE'=1))
    # fix theme
    g = g + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
    # add chromosome lines
    g = g + geom_segment(
      mapping = aes(
        x = starts,
        xend = starts,
        y = -1.5,
        yend = 1.5
      ),
      data = chrsz,
      inherit.aes = F,
      alpha = 1
    )
    # add copy data if present
    ### TUMOR LOGRATIO
    if (exists('cnv_table')) {
      g = g + geom_line(
        aes(x = x, y = y),
        data = cnv_table[cnv_samples[s]],
        alpha = 1 / 3,
        stat = "identity",
        colour = "black"
      )
    }
    
    # add chromosome labels
    g = g + geom_text(
      aes(
        x = starts + 0.5 * length,
        y = -0.7,
        label = label
      ),
      data = chrsz,
      inherit.aes = F
    )
    # add the TUMOR structural variants
    if (exists('manta_tumor_table')) {
      temp = manta_tumor_table[sv_samples[s]][plot == T][IMPRECISE == FALSE]
      if (nrow(temp) > 0)
        g = g + geom_curve(
          aes(
            x = cumstart,
            xend = altcumpos,
            y = 0.5,
            yend = 0.5001,
            col = SVTYPE
          ),
          #,linetype=IMPRECISE),
          data = temp,
          alpha = 1 / 2,
          stat = "identity",
          curvature = -0.2,
          size = 0.5,
          position = "identity",
          angle = 90,
          ncp = 10,
          arrow = arrow(
            length = unit(0.10, "cm"),
            ends = "both",
            type = "open"
          ),
          lineend = "butt",
          na.rm = FALSE,
          show.legend = NA,
          inherit.aes = FALSE
        )
    }
  }
  g
}