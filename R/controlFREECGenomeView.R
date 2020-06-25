controlFREECGenomeView <- function(freec_results) {
  tratio <- freec_results['tratio'][[1]]
  nratio <- freec_results['nratio'][[1]]
  tbaf <- freec_results['tbaf'][[1]]
  nbaf <- freec_results['nbaf'][[1]]
  g = NULL
  if (exists('tratio')) for (sample in unique(tratio$sample)) {
    cat('Control Freec: ',sample)
    par(mai=c(0,4,0,2))
    temp <- tratio[sample]
    temp$Ratio[temp$Ratio>3]=3
    
    g=ggplot2::ggplot()+
      ggplot2::ylab('B allele ratio and Coverage Ratio')+
      ggplot2::xlab(sample)+
      ggplot2::scale_color_gradientn(colours = c('violet','blue','black','orange','red'))+
      ggplot2::expand_limits(x=c(0,3200e6))+
      ggplot2::scale_y_continuous(limits=c(-1,3),
                         breaks = 0:2,
                         minor_breaks = seq(0,2.5,0.5),
                         labels = 0:2
      )
    # normal logR below
    g=g+ggplot2::geom_point(ggplot2::aes(x=cumstart,y=Ratio-0.3),col='darkgrey',data=nratio,alpha=1/5,shape='.')
    # tumor logR
    g=g+ggplot2::geom_point(ggplot2::aes(x=cumstart,y=Ratio,col=CopyNumber),data=temp,alpha=1/5,shape='.')
    # add smoothed
    #g=g+geom_line(aes(x=cumpos,y=tratio),stat="identity",colour="green",size=0.05,data=binned[sample])
    
    g=g+ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                       panel.grid.minor = ggplot2::element_blank(),
                       axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank())
    # add tBAF  <--- downsampled by 10
    g=g+ggplot2::geom_point(data = tbaf[seq(1,nrow(tbaf),10)][sample],
                            ggplot2::aes(x=cumstart,y=-0.5+0.5*BAF),
                            alpha=1/10,shape='.')
    # add nBAF  <--- downsampled by 10
    g=g+ggplot2::geom_point(data = nbaf[seq(1,nrow(nbaf),10)],
                            ggplot2::aes(x=cumstart,y=-1+0.5*BAF),
                            alpha=1/10,shape='.')
    # add chromosome lines
    g=g+ggplot2::geom_segment(mapping = ggplot2::aes(x = starts,xend=starts,y= -1,yend=3),
                              data=chrsz,
                              inherit.aes = F,
                              alpha=1/5)
    # add chromosome labels
    g=g+ggplot2::geom_text(ggplot2::aes(x=starts+0.5*length,y=0.1,label=label),
                           data = chrsz,
                           inherit.aes = F)
  }
  g
}