ascatGenomeView <- function(ascat_results) {
  ascat_tratio <- ascat_results['tratio'][[1]]
  ascat_tbaf <- ascat_results['tbaf'][[1]]
  ascat_nbaf <- ascat_results['nbaf'][[1]]
  g <- NULL
  if (exists('ascat_tratio'))
    for (sample in unique(ascat_tratio$sample)) {
      cat('ASCAT: ', sample)
      par(mai = c(0, 4, 0, 2))
      temp <- ascat_tratio[sample]
      temp$smoothed[temp$smoothed > 3] = 3
      g = ggplot2::ggplot() +
        ggplot2::ylab('B allele ratio and Coverage Ratio') +
        ggplot2::scale_color_gradientn(colours = c('violet', 'blue', 'black', 'orange', 'red')) +
        ggplot2::scale_y_continuous(
          limits = c(-1, 3),
          breaks = 0:2,
          minor_breaks = seq(0, 2.5, 0.5),
          labels = 0:2
        ) +
        ggplot2::expand_limits(x = c(0, 3200e6))
      # using smoothed logR   <--- downsampled by 10
      g = g + ggplot2::geom_point(
        ggplot2::aes(x = cumstart, y = smoothed, color = CopyNumber),
        data = temp[c(1, seq(2, nrow(temp), 10))],
        alpha = 1 / 5,
        shape = '.'
      )
      
      g = g + ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
      # add tBAF  <--- downsampled by 10
      g = g + ggplot2::geom_point(
        mapping = ggplot2::aes(x = cumstart, y = -0.5 + 0.5 * BAF),
        data = ascat_tbaf[seq(1, nrow(ascat_tbaf), 10)][sample],
        alpha = 1 / 10,
        shape = '.'
      )
      # add nBAF  <--- downsampled by 10
      g = g + ggplot2::geom_point(
        mapping = ggplot2::aes(x = cumstart, y = -1 + 0.5 * BAF),
        data = ascat_nbaf[seq(1, nrow(ascat_nbaf), 10)],
        alpha = 1 / 10,
        shape = '.'
      )
      # add chromosome lines
      g = g + ggplot2::geom_segment(
        mapping = ggplot2::aes(
          x = starts,
          xend = starts,
          y = -1,
          yend = 3
        ),
        data = chrsz,
        inherit.aes = F,
        alpha = 1 / 5
      )
      # add chromosome labels
      g = g + ggplot2::geom_text(
        ggplot2::aes(
          x = starts + 0.5 * length,
          y = 0.1,
          label = label
        ),
        data = chrsz,
        inherit.aes = F
      )
    }
  g
}