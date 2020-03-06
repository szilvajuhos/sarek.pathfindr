getInframeHotspots <- function(ih_filename, hotspots_snv) {
  hotspots_inframe <- get("hotspots_inframe", pfenv)
  # this is the local version of near_hotspots
  near_hotspots = NULL
  # if hotspots_inframe is NULL, we assume near_hotspots is also NULL
  if(is.null(hotspots_inframe)) {
    tic("Reading inframe hotspots")
    hotspots_inframe = unique(fread(ih_filename)[,.(Hugo_Symbol,Amino_Acid_Position)])
    hotspots_inframe$start = as.numeric(
      str_replace(
        string = hotspots_inframe$Amino_Acid_Position,
        pattern = '-[0-9]*',
        replacement = ''
      )
    )
    hotspots_inframe$end = as.numeric(
      str_replace(
        string = hotspots_inframe$Amino_Acid_Position,
        pattern = '[0-9]*-',
        replacement = ''
      )
    )
    for (i in 1:nrow(hotspots_inframe))
      near_hotspots = c(near_hotspots,
                        paste(
                          hotspots_inframe$Hugo_Symbol[i],
                          c(
                            seq(hotspots_inframe$start[i] - 2, hotspots_inframe$start[i] + 2),
                            hotspots_inframe$end[i] - 2,
                            hotspots_inframe$end[i] + 2
                          )
                        ))
    for (i in 1:nrow(hotspots_snv))
      near_hotspots = c(near_hotspots,
                        paste(
                          hotspots_snv$Hugo_Symbol[i],
                          seq(hotspots_snv$pos[i] - 2, hotspots_snv$pos[i] +
                                2)
                        ))
    toc()
  }
  assign(x = "near_hotspots", value = near_hotspots, envir = pfenv)
  assign(x = "hotspots_inframe", value = hotspots_inframe, envir = pfenv)
  get("hotspots_inframe",pfenv)
}