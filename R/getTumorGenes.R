getTumourGenes <- function(tg_filename, ltg_filename) {
  # get the variable from the package-wide environment
  tumorgenes = get("tumorgenes", envir = pfenv)
  # re-read the file only if it is not defined yet
  if(is.null(tumorgenes)) {
    tic(paste0("Reading tumor genes from file ",tg_filename))
    tumorgenes=data.table::fread(tg_filename,key='Gene Symbol')
    tumorgenes$chr=paste0('chr',str_replace(string=tumorgenes$`Genome Location`,pattern = ':.*',replacement = ''))
    temp=str_replace(string=tumorgenes$`Genome Location`,pattern = '.*:',replacement = '')
    tumorgenes$start=as.numeric(str_replace(string=temp,pattern = '-.*',replacement = ''))
    tumorgenes$end=as.numeric(str_replace(string=temp,pattern = '.*-',replacement = ''))
    tumorgenes$cumstart=NA
    tumorgenes$cumend=NA
    for (i in 1:nrow(chrsz)) {
      ix <- tumorgenes$chr==chrsz$chr[i]
      tumorgenes$cumstart[ix] <- tumorgenes$start[ix]+chrsz$starts[i]
      tumorgenes$cumend[ix] <- tumorgenes$end[ix]+chrsz$starts[i]
    }
    assign("tumorgenes", tumorgenes, envir = pfenv)
    toc()
    
    tic("Merging local genes and selecting TSGs")
    local_tumorgenes = data.table::fread(ltg_filename,key='Gene')
    # get rid of rows without a gene name
    local_tumorgenes = subset.data.frame(local_tumorgenes,local_tumorgenes$Gene!='')
    # select only the gene names and the tier (a.k.a. score)
    local_tumorgenes = local_tumorgenes[,1:2]
    alltumorgenes=unique(c(tumorgenes$`Gene Symbol`,local_tumorgenes$Gene))
    alltsg=tumorgenes[grep('TSG',tumorgenes$`Role in Cancer`),tumorgenes$`Gene Symbol`]

    assign("alltumorgenes", alltumorgenes, envir = pfenv)
    assign("alltsg", alltsg, envir = pfenv)
    toc()
    
    # assign("allfusion", allfusion, envir = pfenv)
    # assign("allfusionpairs", allfusionpairs, envir = pfenv)
    # assign("alltier1", alltier1, envir = pfenv)
    # assign("alltier2", alltier1, envir = pfenv)
  }
}