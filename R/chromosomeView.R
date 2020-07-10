chromosomeView <- function(thischr, ascat_results, freec_results, manta_tumor_table, strelka_table, haplotypecaller_selected) {
  ascat_binned <- ascat_results['binned'][[1]]
  
  samplenames=sort(unique(ascat_binned$sample))
  nsamples=length(samplenames)
  g<-NULL
  for (s in 1:nsamples) {
    cat(samplenames[s])
    par(mai=c(0,4,0,2))
    g=ggplot()+
      scale_y_continuous(limits=c(-2.7,3),
                         breaks = seq(1.5,-2.5,-0.5),
                         minor_breaks = seq(1.25,-2.25,-0.5),
                         labels =
                           c('Tumor DNA=1.5',
                             'Tumor Copy Ratio\n(median centered)',
                             'Tumor DNA=0.5\nMutation AF=1',
                             'Mutation AF=0.5',
                             'Mutation AF=0\nTumor SNP AF=1',
                             'Tumor\nSNPs',
                             'Tumor SNP AF=0\nNormal SNP AF=1',
                             'Normal\nSNPs',
                             'Normal SNP AF=0'))+
      scale_x_continuous(
        breaks = 10e6*seq(0,ceiling(chrsz$length[chrsz$chr==thischr]/10e6)),
        labels= c('0',paste0(10*seq(1,ceiling(chrsz$length[chrsz$chr==thischr]/10e6)),'M')),
        minor_breaks = 1e6*seq(0,ceiling(chrsz$length[chrsz$chr==thischr]/1e6)))+
      theme(panel.grid.minor.x = element_line(colour="lightgrey", size=0.05),
            panel.grid.major.x = element_line(colour="grey", size=0.2),
            panel.grid.major.y = element_line(colour="black", size=0.2),
            panel.grid.minor.y = element_line(colour="lightgrey", size=0.05),
            axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank())#+
    
    # lines to separate parts of fig
    g=g+geom_segment(mapping = aes(x = 0,xend=chrsz$length[chrsz$chr==thischr],y=0.5,yend=0.5),alpha=1)
    g=g+geom_segment(mapping = aes(x = 0,xend=chrsz$length[chrsz$chr==thischr],y=-0.5,yend=-0.5),alpha=1)
    g=g+geom_segment(mapping = aes(x = 0,xend=chrsz$length[chrsz$chr==thischr],y=-1.5,yend=-1.5),alpha=1)
    g=g+geom_segment(mapping = aes(x = 0,xend=chrsz$length[chrsz$chr==thischr],y=-2.5,yend=-2.5),alpha=1)

    ### TUMOR (& normal) LOGRATIO
    if( !is.null(freec_results) ) {
      freec_tratio <- freec_results['tratio'][[1]]
      freec_nratio <- freec_results['nratio'][[1]]

      temp <- freec_nratio
      temp$Ratio[temp$Ratio>3]=3
      g=g+geom_line(aes(x=Start,y=Ratio),
                    data = temp[temp$Chromosome==thischr,],alpha=1/3,
                    stat="identity",colour="red")
      samples=sort(unique(freec_tratio$sample))
      temp <- freec_tratio[samples[s]]
      temp$Ratio[temp$Ratio>3]=3
      g=g+geom_line(aes(x=Start,y=Ratio),
                    data = temp[temp$Chromosome==thischr,],alpha=2/3,
                    stat="identity",colour="black")
    } else if (!is.null(ascat_results)) {
      ascat_tratio <- ascat_results['tratio'][[1]]
      samples=sort(unique(ascat_tratio$sample))
      temp <- ascat_tratio[samples[s]]
      temp$smoothed[temp$smoothed>3]=3
      g=g+geom_line(aes(x=Position,y=smoothed),
                    data = temp[temp$Chr==thischr,],alpha=1/3,
                    stat="identity",colour="black")
    }

    ### TUMOR STRVAR
    if (!is.null(manta_tumor_table)) if (nrow(manta_tumor_table)>0) {
      # rm(chr) why is this here?
      samples=sort(unique(manta_tumor_table$sample))
      # first same chrom SVs
      temp=unique(manta_tumor_table[sample==samples[s] & chr==thischr & chr==altchr,
                                    .(sample,chr,altchr,start,altpos,IMPRECISE,SVTYPE)])
      if (nrow(temp)>0) g=g+geom_curve(aes(x = start,xend = altpos,y=1.99,yend=2,col=SVTYPE,linetype=IMPRECISE),
                                       data = temp,alpha=1/2,
                                       curvature=-0.2,size=1,
                                       position = "identity", angle = 90, ncp = 50,
                                       arrow = arrow(length=unit(0.10,"cm"), ends="both", type = "open"),
                                       lineend = "butt", na.rm = FALSE, show.legend = NA,
                                       inherit.aes = FALSE)+
        scale_color_manual(values=c("BND"="#E41A1C", "DEL"="#377EB8", "INV"="#4DAF4A", "INS"="#984EA3", "DUP"="#FF7F00"))+
        scale_linetype_manual(values=c('TRUE'=2,'FALSE'=1))

      # then diff chrom SVs with other end towards left
      temp=unique(manta_tumor_table[sample==samples[s] & chr==thischr & altcumpos<cumstart & chr!=altchr,
                                    .(sample,chr,altchr,start,IMPRECISE,SVTYPE)])
      if (nrow(temp)>0) {
        g=g+geom_curve(aes(x = 0,xend = start,y=3,yend=2,col=SVTYPE,linetype=IMPRECISE),
                       data = temp,alpha=1/2,
                       curvature=-0.2,size=1,
                       position = "identity", angle = 90, ncp = 50,
                       arrow = arrow(length=unit(0.10,"cm"), ends="both", type = "open"),
                       lineend = "butt", na.rm = FALSE, show.legend = NA,
                       inherit.aes = FALSE)
      }
      # then diff chrom SVs with other end towards right
      temp=unique(manta_tumor_table[sample==samples[s] & chr==thischr & altcumpos>cumstart & chr!=altchr,
                                    .(sample,chr,altchr,start,IMPRECISE,SVTYPE)])

      if (nrow(temp)>0) {
        g=g+geom_curve(aes(x = start,xend = chrsz$length[chrsz$chr==thischr],y=2,yend=3,col=SVTYPE,linetype=IMPRECISE),
                       data = temp,alpha=1/2,
                       curvature=-0.2,size=1,
                       position = "identity", angle = 90, ncp = 50,
                       arrow = arrow(length=unit(0.10,"cm"), ends="both", type = "open"),
                       lineend = "butt", na.rm = FALSE, show.legend = NA,
                       inherit.aes = FALSE)
      }
    }

     
    ### TUMOR MUTATIONS
    alltier1 <- getEnvVariable('alltier1')
    alltier2 <- getEnvVariable('alltier2')
    
    samples=sort(unique(strelka_table$sample))
    temp=unique(strelka_table[on=samples[s]][chr==thischr & type!='other',][,.(start,AFreq,type,Amino_acids,SYMBOL,IMPACT,Protein_position)])
    if (nrow(temp)>0) {
      g=g+geom_point(aes(x=start,y=AFreq -0.5,fill=type,shape=type),data = temp,alpha=1)+
        scale_fill_manual(values = c("A>C"="#1B9E77","A>G"="#D95F02", "A>T"="#7570B3", "C>A"="#E7298A", "C>G"="#66A61E", "C>T"="#E6AB02", "del"="#A6761D", "ins"="#666666"))+
        scale_shape_manual(values=c("A>C"=21,"A>G"=21, "A>T"=21, "C>A"=21, "C>G"=21, "C>T"=21, "del"=24, "ins"=25))#+

      # with labels
      g=g+geom_label_repel(data = temp[IMPACT %in% c('MODERATE','HIGH') & SYMBOL %in% alltier2,],
                           mapping = aes(x=start,y=AFreq -0.5,
                                         label=paste(SYMBOL,Protein_position,Amino_acids)),
                           inherit.aes = F,size=3,min.segment.length = 0,nudge_y = 0.3,col='grey')
      g=g+geom_label_repel(data = temp[IMPACT %in% c('MODERATE','HIGH') & SYMBOL %in% alltier1,],
                           mapping = aes(x=start,y=AFreq -0.5,
                                         label=paste(SYMBOL,Protein_position,Amino_acids)),
                           inherit.aes = F,size=3,min.segment.length = 0,nudge_y = 0.3,col='black')
    }

    ### SNPS
    if (!is.null(freec_results)) {
      tbaf <- freec_results['tbaf'][[1]]
      nbaf <- freec_results['nbaf'][[1]]
      
      samples=sort(unique(tbaf$sample))
      g=g+geom_point(aes(x=Position,y=BAF -1.5),data = tbaf[sample==samples[s]][Chromosome==thischr],alpha=1/8,shape='.')
      g=g+geom_point(aes(x=Position,y=BAF -2.5),data = nbaf[nbaf$Chromosome==thischr,],alpha=1/8,shape='.')
    } else if (!is.null(ascat_results) )  {
      ascat_tbaf <- ascat_results['tbaf'][[1]]
      ascat_nbaf <- ascat_results['nbaf'][[1]]
      
      samples=sort(unique(ascat_tbaf$sample))
      g=g+geom_point(aes(x=Position,y=BAF -1.5),data = ascat_tbaf[sample==samples[s]][Chr==thischr],alpha=1/8,shape='.')
      g=g+geom_point(aes(x=Position,y=BAF -2.5),data = ascat_nbaf[ascat_nbaf$Chr==thischr,],alpha=1/8,shape='.')
    }

    ## Germline muts
    if (!is.null(haplotypecaller_selected)) {
      alltumorgenes <- getEnvVariable('alltumorgenes')
      hc_files <- haplotypeCallerFiles()
      haplotypecaller_N_file <- hc_files[1]
      if (!is.null(hc_files)) {
        haplotypecaller_T_file <- hc_files[2]
        sample=strsplit(basename(haplotypecaller_T_file[s]),'[.]')[[1]][1]
        temp=haplotypecaller_selected[sample==sample][chr==thischr][rank_score>1][type!='other']
        temp=unique(temp[,.(start,type,AFreq,Amino_acids,Protein_position,SYMBOL)])
        g=g+
          geom_point(aes(x=start,y=AFreq -1.5,fill=type,shape=type),
                     data = temp,alpha=1)+
          geom_label_repel(data = temp[Amino_acids!='' & SYMBOL %in% alltumorgenes],
                           mapping = aes(x=start,y=AFreq -1.5,label=paste(SYMBOL,Protein_position,Amino_acids)),
                           inherit.aes = F,size=3,min.segment.length = 0)
        normal=strsplit(basename(haplotypecaller_N_file),'[.]')[[1]][1]
        temp=haplotypecaller_selected[sample==normal][chr==thischr][rank_score>1][type!='other']
        temp=unique(temp[,.(start,type,AFreq,Amino_acids,Protein_position,SYMBOL)])
        g=g+
          geom_point(aes(x=start,y=AFreq -2.5,fill=type,shape=type),
                     data = temp,alpha=1)+
          geom_label_repel(data = temp[Amino_acids!='' & SYMBOL %in% alltumorgenes],
                           mapping = aes(x=start,y=AFreq -2.5,label=paste(SYMBOL,Protein_position,Amino_acids)),
                           inherit.aes = F,size=3,min.segment.length = 0)
        
      }
    }

    # add gene names
    PFconfig<-getEnvVariable('PFconfig')
    getTumorGenes(PFconfig$tumorgenes, PFconfig$local_tumorgenes)
    tumorgenes <- getEnvVariable('tumorgenes')
    temp=tumorgenes[chr==thischr]
    if (nrow(temp) > 0)
      g = g + geom_text_repel(
        mapping = aes(
          x = (start + end) / 2,
          y =  -2.6,
          label = `Gene Symbol`
        ),
        data = temp,
        nudge_y = -0.1,
        size = 2
      )
    print(g) 
  }
}
