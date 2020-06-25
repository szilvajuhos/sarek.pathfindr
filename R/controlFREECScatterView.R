controlFREECScatterView <- function(binned) {
  g<-NULL
  if (exists('binned')) for (sample in unique(binned$sample)) {
    cat('Control Freec: ',sample)
    binned$tmaf[binned$tmaf<0]=1
    g=ggplot2::ggplot()+
      ggplot2::expand_limits(y=c(.5,1))+
      ggplot2::expand_limits(x=c(.3,3))+
      ggplot2::geom_point(ggplot2::aes(x = tratio,y = tmaf,col=chr), data=binned[sample][tratio<3])+
      ggplot2::xlab('Tumor (1Mb mean) DNA ratio')+
      ggplot2::ylab('Tumor (1Mb mean) major allele ratio')+
      ggplot2::scale_y_continuous(breaks = seq(0.5,1,0.1)) +
      ggplot2::scale_x_continuous(breaks = seq(0.5,2.5,0.1))+
      ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_line(colour="grey", size=0.2),
            panel.grid.major.y = ggplot2::element_line(colour="grey", size=0.2),
            panel.grid.minor.y = ggplot2::element_blank()
      )
  }
  g
}