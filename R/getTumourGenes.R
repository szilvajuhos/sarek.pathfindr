getTumourGenes <- function(tg_filename, ltg_filename) {
    # get the variable from the package-wide environment
    tg = get("tumorgenes", envir = pfenv)
    # re-read the file only if it is not defined yet
    if(is.null(tg)) {
      tic(paste0("Reading tumor genes from file ",tg_filename))
      tg = data.table::fread(tg_filename,key='Gene Symbol')
      tg$chr=paste0('chr',str_replace(string=tg$`Genome Location`,pattern = ':.*',replacement = ''))
      temp=str_replace(string=tg$`Genome Location`,pattern = '.*:',replacement = '')
      tg$start=as.numeric(str_replace(string=temp,pattern = '-.*',replacement = ''))
      tg$end=as.numeric(str_replace(string=temp,pattern = '.*-',replacement = ''))
      tg$cumstart=NA; 
      tg$cumend=NA
      for (i in 1:nrow(chrsz)) {
        ix <- tg$chr==chrsz$chr[i]
        tg$cumstart[ix] <- tg$start[ix]+chrsz$starts[i]
        tg[ix] <- tg$end[ix]+chrsz$starts[i]
      }
      toc()
      
      tic("Merging genes from local list")
      local_tumorgenes = data.table::fread(ltg_filename,key='Gene')
      # get rid of rows without a gene name
      local_tumorgenes = subset.data.frame(local_tumorgenes,local_tumorgenes$Gene!='')
      # select only the gene names and the tier (a.k.a. score) 
      local_tumorgenes = local_tumorgenes[,1:2]
      alltumorgenes=unique(c(tg$`Gene Symbol`,local_tumorgenes$Gene))
      alltsg=tg[grep('TSG',tg$`Role in Cancer`),tg$`Gene Symbol`]
      toc()
      
      tic("Getting fusions")
      allfusion=tg[grep('fusion',`Role in Cancer`),`Gene Symbol`]
      allfusionpairs = NULL
      for (i in 1:length(allfusion)) {
        t = trimws(strsplit(tg[allfusion[i], `Translocation Partner`], ', ')[[1]])
        if (length(t) > 0)
          for (j in 1:length(t))
            allfusionpairs = c(allfusionpairs, paste(sort(c(allfusion[i], t[j])), collapse = ' '))
      }
      toc()
      tic("Creating tiers")
      # select tier 1 genes from both
      # then get a list of gene names as the union of the two sets
      tumorgenes_tier1 = subset.data.frame(tg,tg$Tier==1)
      local_tumorgenes_tier1 = subset.data.frame(local_tumorgenes,local_tumorgenes$`Tier 1 and 2 for pediatric cancers final`==1)
      # select the first column (gene names)
      tumorgenes_tier1_genes = tumorgenes_tier1[,1]
      local_tumorgenes_tier1_genes = local_tumorgenes_tier1[,1]
      alltier1 = union(tumorgenes_tier1_genes,local_tumorgenes_tier1_genes)

      # DRY but whatever, do the same for tier 2 TODO make it nice without repetition
      tumorgenes_tier2 = subset.data.frame(tg,tg$Tier==2)
      local_tumorgenes_tier2 = subset.data.frame(local_tumorgenes,local_tumorgenes$`Tier 1 and 2 for pediatric cancers final`==2)
      alltier2 = union(tumorgenes_tier2[,1],local_tumorgenes_tier2[,1])
      toc()
      
      assign("tumorgenes", tg, envir = pfenv)
      assign("allfusion", allfusion, envir = pfenv)
      assign("allfusionpairs", allfusionpairs, envir = pfenv)
      assign("alltier1", alltier1,envir = pfenv)
      assign("alltier2", alltier1,envir = pfenv)
      toc()
    }
}