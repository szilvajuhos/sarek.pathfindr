strelkaView<-function(strelka_table) {
  cat(PFconfig$strelka_indel_file,PFconfig$strelka_snv_file,sep = '\n')
  g_strelka <- NULL
  alltier1 <- get("alltier1",pfenv)
  alltier2 <- get("alltier2",pfenv)
  # TODO: DRY out with mutect2View
  if (exists('strelka_table')) for (sample in sort(unique(strelka_table$sample))) {
    par(mai=c(1,4,0,2))
    
    # Basic plot element
    g_basic=
      ggplot()+
      ylab('Somatic mutation allele ratio')+
      expand_limits(y=c(0,1))+
      expand_limits(x=c(0,3210e6))+
      scale_fill_manual(values = c("A>C"="#1B9E77","A>G"="#D95F02", "A>T"="#7570B3", "C>A"="#E7298A", "C>G"="#66A61E", "C>T"="#E6AB02", "del"="#A6761D", "ins"="#666666"))+
      scale_shape_manual(values=c("A>C"=21,"A>G"=21, "A>T"=21, "C>A"=21, "C>G"=21, "C>T"=21, "del"=24, "ins"=25))+
      scale_color_manual(values=c('<5'='lightgrey','≥5'='black'))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      geom_segment(mapping = aes(x = starts,xend=starts,y=0,yend=1),data=chrsz,inherit.aes = F,alpha=1/5)+
      geom_text(aes(x=starts+0.5*length,y=0,label=label),data = chrsz,inherit.aes = F)
    
    
    # Strelka plot
    cat('Strelka: ',sample)
    temp=unique(strelka_table[on=sample][,.(cumstart,reads,AFreq,type,SYMBOL,Protein_position,Amino_acids,IMPACT)])
    g_strelka=g_basic+
      xlab(paste('Strelka: ',sample))+
      geom_point(data = temp,mapping = aes(x=cumstart,y=AFreq,fill=type,shape=type,col=reads),alpha=1)+
      geom_label_repel(data = temp[SYMBOL %in% alltier1 & IMPACT %in% c('MODERATE','HIGH'),],
                       mapping = aes(x=cumstart,y=AFreq,label=paste(SYMBOL,Protein_position,Amino_acids)),
                       inherit.aes = F,size=3,nudge_y = 0.2,col='black')+
      geom_label_repel(data = temp[SYMBOL %in% alltier2 & IMPACT %in% c('MODERATE','HIGH'),],
                       mapping = aes(x=cumstart,y=AFreq,label=paste(SYMBOL,Protein_position,Amino_acids)),
                       inherit.aes = F,size=3,nudge_y = 0.2,col='grey')
  }
  g_strelka
}