scoreStrelka <- function() {
  cat(" ------------ Score Strelka function -------------\n")
  PFconfig<-getEnvVariable('PFconfig')
  strelka_result_files <- strelkaFiles()
  cat('Selected files for Strelka in', getwd(), ':\n\n')
  cat("Strelka Somatic SNVs file:  ", strelka_result_files["strelka_snv_file"], '\n')
  cat("Strelka Somatic indels file:  ", strelka_result_files["strelka_indel_file"], '\n')
  sample = strsplit(basename(strelka_result_files["strelka_snv_file"]), '_somatic')[[1]][1]
  
  strelka_selected <- NULL
  # get data tables from package environment, or read them from files if they are not read yet
  hotspots_snv <- getSNVhotspots(PFconfig$hotspots_snv)
  # defining inframe also defines near_hotspots
  hotspots_inframe <- getInframeHotspots(PFconfig$hotspots_inframe,hotspots_snv)
  near_hotspots <- get("near_hotspots",pfenv)
  
  snptable <- getSNPtable(PFconfig$snptable)
  cosmic_coding <- getCountTable(PFconfig$coding_table,"coding_table")
  cosmic_noncoding <- getCountTable(PFconfig$noncoding_table,"noncoding_table")

  # this function creates more tables, get them via get(tablename,pfenv)
  getTumorGenes(PFconfig$tumorgenes, PFconfig$local_tumorgenes)
  tumorgenes <- get("tumorgenes",pfenv)
  alltumorgenes <- get("alltumorgenes",pfenv)
  # get TSGs
  alltsg <- get("alltsg",pfenv)
  # get tier data
  alltier1 <- get("alltier1",pfenv)
  alltier2 <- get("alltier2",pfenv)

  tic("Strelka indels")
  indels <- scoreStrelkaIndels(strelka_result_files["strelka_indel_file"], sample)
  toc()
  
  tic("Strelka SNVs")
  snvs <- scoreStrelkaSNVs(strelka_result_files["strelka_snv_file"], sample)
  toc()
  
  table_snvs <- snvs[["table_snvs"]]
  table_indels <- indels[["table_indels"]]
  
  # Concatenate the snvs and indels tables for this patient and sample
  table_both = rbind(table_snvs, table_indels)
  table_both$Swegen_count = NA
  table_both$Cosmic_count = NA
  table_both$rank_score = 0
  table_both$rank_terms = ''
  table_both$LOH = ''
  
  # Variant type
  table_both$type = paste0(table_both$ref, '>', table_both$alt)
  table_both$type[nchar(table_both$ref) > nchar(table_both$alt)] = 'del'
  table_both$type[nchar(table_both$ref) < nchar(table_both$alt)] = 'ins'
  table_both$type[nchar(table_both$type) > 3] = 'other'
  table_both$type[table_both$type == 'T>G'] = 'A>C'
  table_both$type[table_both$type == 'T>C'] = 'A>G'
  table_both$type[table_both$type == 'T>A'] = 'A>T'
  table_both$type[table_both$type == 'G>T'] = 'C>A'
  table_both$type[table_both$type == 'G>C'] = 'C>G'
  table_both$type[table_both$type == 'G>A'] = 'C>T'
  
  # concat CSQ lists
  csq <- c(snvs[["snv_csq"]], indels[["indel_csq"]])
  table_both$CSQ = csq
  
  cat("get VEP headers from a vcf object (assumed never to differ between snvs and indels)\n")
  snv_vcf <- snvs[["snv_vcf"]]
  vep_header = strsplit(info(header(snv_vcf))['CSQ', ][, 3], 'Format: ')[[1]][2]
  vep_header = strsplit(vep_header, '\\|')[[1]]
  
  cat("VEP annotation is put in annotation_table\n")
  annotation_table = matrix(
    data = NA,
    nrow = length(unlist(csq)),
    ncol = length(vep_header) + 1
  )
  colnames(annotation_table) = c('ID', vep_header)
  row = 1
  for (i in 1:length(csq)) {
    # for each variant
    # for each VEP annotation:
    for (j in 1:length(csq[[i]])) {
      line = strsplit(csq[[i]][j], '\\|')[[1]]
      annotation_table[row, 1] = table_both$ID[i]
      annotation_table[row, 1 + (1:length(line))] = line
      row = row + 1
    }
  }
  # Annotation table then merged with the "table_both"
  strelka_table = NULL
  strelka_table = rbind(strelka_table, merge(table_both, as.data.table(annotation_table), by =
                                               'ID'))
  setkey(strelka_table, 'sample')
  
  cat("done parsing each sample\n")
  cat("Add Control Freec LOH ... slow. (well, it is not working yet) ")
  # TODO: CF LOH data is not available from this point yet.
  tic("Control-FREEC LOH")
  try({
    if (nrow(strelka_table) > 0 &
        !is.null(freec_loh))
      for (i in 1:nrow(strelka_table)) {
        ix = freec_loh$chr == strelka_table$chr[i] &
          freec_loh$start < strelka_table$start[i] &
          freec_loh$end > strelka_table$end[i]
        ix = ix[!is.na(ix)]
        if (sum(ix) > 0)
          strelka_table$LOH[i] = 'Y'
      }
  }, silent = T)
  toc()
  
  by_pos = snptable[strelka_table$ID, value] # if the variant is as 10:10001801_C/T in database
  by_name1 = snptable[str_extract(strelka_table$Existing_variation, "^rs[0-9]+"), value] # if first rsid
  by_name2 = snptable[str_extract(strelka_table$Existing_variation, "rs[0-9]+$"), value] # if last rsid (often both)
  d = data.table(by_pos, by_name1, by_name2, default = 0)
  d$max = apply(
    X = d,
    MARGIN = 1,
    FUN = max,
    na.rm = T
  )
  strelka_table$Swegen_count = d$max
  
  # TODO: ask Markus about this
  key = paste0(substr(strelka_table$chr, 4, 6),
               ':',
               strelka_table$start,
               '-',
               strelka_table$end)
  # if (project!='BTB')
  key = paste0(strelka_table$chr,
               ':',
               strelka_table$start,
               '-',
               strelka_table$end)
  key = str_replace(key, 'X:', '23:')
  key = str_replace(key, 'Y:', '24:')
  counts = cbind(cosmic_coding[key, value], cosmic_noncoding[key, value], 0)
  max_ = apply(counts, 1, max, na.rm = T)
  strelka_table$Cosmic_count = max_
  
  strelka_table$cumstart = strelka_table$start
  strelka_table$cumend = strelka_table$end
  cat("for each chr get cumulative pos\n")
  for (i in 1:nrow(chrsz)) {
    ix = strelka_table$chr == chrsz$chr[i]
    strelka_table$cumstart[ix] = strelka_table$start[ix] + chrsz$starts[i]
    strelka_table$cumend[ix] = strelka_table$end[ix] + chrsz$starts[i]
  }
  selection = strelka_table[, -c('CSQ')] # not needed after parsing/merging the annotations

  # if there are any in the selection
  if (nrow(selection) > 0) {
    tic("Selecting tiers")
    # Cosmic/local Tier2 :
    ix = selection$SYMBOL %in% alltier2
    selection$rank_score[ix] = 2
    selection$rank_terms[ix] = 'T2_gene'
    # tier 1 priority:
    ix = selection$SYMBOL %in% alltier1
    selection$rank_score[ix] = 2
    selection$rank_terms[ix] = 'T1_gene'
    toc()
    
    tic("Add high impact")
    ix = selection$IMPACT == 'HIGH'
    if (any(ix)) {
      selection$rank_score[ix] = selection$rank_score[ix] + 2
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'high_impact')
    }
    toc()
    tic("Add moderate impact")
    ix = selection$IMPACT == 'MODERATE'
    if (any(ix)) {
      selection$rank_score[ix] = selection$rank_score[ix] + 1
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'moderate_impact')
    }
    toc()
    tic("additional +1 if high impact and TSG")
    ix = selection$IMPACT == 'HIGH' & selection$SYMBOL %in% alltsg
    if (any(ix)) {
      selection$rank_score[ix] = selection$rank_score[ix] + 1
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'high+TSG')
    }
    toc()
    tic("Add clinvar pathogenic")
    ix = grep('pathogenic', selection$CLIN_SIG)
    if (length(ix) > 0) {
      selection$rank_score[ix] = selection$rank_score[ix] + 2
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'clinvar')
    }
    toc()
    tic("Add Polyphen/SIFT damaging/deleterious")
    ix = union(grep('damaging', selection$PolyPhen),
               grep('deleterious', selection$SIFT))
    if (length(ix) > 0) {
      selection$rank_score[ix] = selection$rank_score[ix] + 1
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'polyphen/SIFT')
    }
    toc()
    tic("Add hotspots and near hotspots (within 2 residues of a hotspot)")
    key = paste(
      selection$SYMBOL,
      str_replace(
        string = selection$Protein_position,
        pattern = '/.*',
        replacement = ''
      )
    )
    ix <-
      key %in% paste(hotspots_snv[, Hugo_Symbol], hotspots_snv[, Amino_Acid_Position]) |
      key %in% paste(hotspots_inframe[, Hugo_Symbol], hotspots_inframe[, Amino_Acid_Position])  ## Warning: Exact match used with inframes
    if (any(ix)) {
      selection$rank_score[ix] = selection$rank_score[ix] + 2
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'hotspot')
    }
    ix = !ix &
      key %in% near_hotspots # the near_hotspot excludes those that were a hotspot
    if (any(ix)) {
      selection$rank_score[ix] = selection$rank_score[ix] + 1
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'near_hotspot')
    }
    toc()
    tic("Add cosmic counts")
    ix = which(selection$Cosmic_count > 50)
    if (length(ix) > 0) {
      selection$rank_score[ix] = selection$rank_score[ix] + 2
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'cosmic_>50')
    }
    ix = which(selection$Cosmic_count > 5 &
                 selection$Cosmic_count <= 50)
    if (length(ix) > 0) {
      selection$rank_score[ix] = selection$rank_score[ix] + 1
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'cosmic_>5')
    }
    toc()
    tic("Special case for TERT promoter")
    ix = which(selection$SYMBOL == 'TERT' &
                 selection$Consequence == '5_prime_UTR_variant')
    if (length(ix) > 0) {
      selection$rank_score[ix] = selection$rank_score[ix] + 4
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'TERT_upstr')
    }
    toc()
    tic("Add TF binding variants near (100kb) cancer genes")
    ix = grep('TF', selection$Consequence)
    if (length(ix) > 0) {
      selection$rank_score[ix] = selection$rank_score[ix] + 2
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'TFBS')
      for (i in 1:length(ix)) {
        near = which(
          selection$cumend[ix[i]] > (tumorgenes$cumstart - 100e3) &
            selection$cumstart[ix[i]] < (tumorgenes$cumend + 100e3) &
            !selection$SYMBOL[ix[i]] %in% alltumorgenes
        )
        if (length(near) > 0) {
          genes = paste(tumorgenes$`Gene Symbol`[near], collapse = ',')
          selection$rank_score[ix[i]] = selection$rank_score[ix[i]] + 2
          selection$rank_terms[ix[i]] = paste(selection$rank_terms[ix[i]], paste0('near_', genes))
        }
      }
    }
    ix = which(selection$CANONICAL != 'YES')
    if (length(ix) > 0) {
      selection$rank_score[ix] = selection$rank_score[ix] - 0
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'not_canonical')
    }
    
    ix = which(as.numeric(selection$CADD_PHRED) > 30)
    if (length(ix) > 0) {
      selection$rank_score[ix] = selection$rank_score[ix] + 2
      selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'cadd_>30')
    }
    toc()
    
  }
  firstcols = c(
    'ID',
    'sample',
    'SYMBOL',
    'rank_score',
    'rank_terms',
    'LOH',
    'AFreq',
    'Consequence',
    'IMPACT',
    'CADD_PHRED',
    'SWAF',
    'TOPMED',
    'Swegen_count'
  )
  cols = colnames(selection)
  setcolorder(x = selection, neworder = c(firstcols, cols[!cols %in% firstcols]))
  strelka_selected <-
    selection[order(cumstart, Allele)][order(rank_score, decreasing = T)]
  strelka_selected <- strelka_selected[rank_score > getRankThreshold("strelka_threshold") ]
  strelka_csv_file_name <- paste0(csv_dir, '/', sample, '_tumor.csv')
  cat("Working dir is ",getwd(),"\n")
  cat("Writing results to ", strelka_csv_file_name, "\n")
  fwrite(strelka_selected, file = strelka_csv_file_name)
  cat("Ready.\n")
  strelka_selected
}
