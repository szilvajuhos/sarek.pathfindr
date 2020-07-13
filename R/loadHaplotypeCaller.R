loadHaplotypeCaller <- function(haplotypecaller_files) {
  haplotypecaller_selected <- NULL
  haplotypecaller_ids <- NULL
  reference_genome<-getEnvVariable('reference_genome')
  normal_sample <- strsplit(basename(haplotypecaller_files[1]), '[.]')[[1]][1]
  if (!is.null(haplotypecaller_files)) {
    haplotypecaller_table = NULL
    for (s in 1:length(haplotypecaller_files)) {
      tic("Reading ",haplotypecaller_files[s])
      vcf = readVcf(file = haplotypecaller_files[s], genome = reference_genome)
      sample = strsplit(basename(haplotypecaller_files[s]), '[.]')[[1]][1]
      haplotypecaller_ids[[sample]] = names(vcf)
      
      # if (s>1) # if not the normal, keep only IDs already present (the normal is first)  <--- this is why somatic variants are not shown
      #   vcf=vcf[names(vcf) %in% haplotypecaller_table$ID]
      
      # pre-filtering by SWAF / TOPMED to avoid extreme number of variants
      swaf = as.data.table(info(vcf)$SWAF)
      topmed = as.data.table(info(vcf)$TOPMED)
      values = data.table(swaf$value, topmed$value, 0) # where both are NA, 0 will be kept as AF
      values$max = apply(values, 1, max, na.rm = T)
      keep = unique(swaf$group[which(values$max < 0.01)])  # variant IXs with all alleles <1% in SWAF and TOPMED (swaf$group equals topmed$group)
      vcf = vcf[keep]
      g = geno(vcf)
      inf = (info(vcf))
      rr = rowRanges(vcf)
      
      mutations = as.data.table(ranges(rr))
      mutations$chr <- as.character(seqnames(rr))
      mutations$sample = sample
      mutations = mutations[, .(ID = names, sample, chr, start, end, width)]
      
      mutations$ref = as.character(ref(vcf))
      mutations$alt = as.data.table(alt(vcf))[, .(values = list(value)), by = group][, values]
      mutations$type = paste0(mutations$ref, '>', mutations$alt)
      mutations$type[nchar(mutations$ref) > nchar(mutations$alt)] = 'del'
      mutations$type[nchar(mutations$ref) < nchar(mutations$alt)] = 'ins'
      mutations$type[nchar(mutations$type) > 3] = 'other'
      mutations$type[mutations$type == 'T>G'] = 'A>C'
      mutations$type[mutations$type == 'T>C'] = 'A>G'
      mutations$type[mutations$type == 'T>A'] = 'A>T'
      mutations$type[mutations$type == 'G>T'] = 'C>A'
      mutations$type[mutations$type == 'G>C'] = 'C>G'
      mutations$type[mutations$type == 'G>A'] = 'C>T'
      
      mutations$AFreq = round(sapply(g$AD, "[[", 2) / unlist(g$DP), 2) # not a perfect estimate...
      ad = as.data.table(g$AD)
      colnames(ad) = 'AD'
      mutations = cbind(mutations, ad)

      # add counts from the new sources
      suppressWarnings({
        t = unlist(lapply(
          X = inf$CAF,
          FUN = max,
          na.rm = T
        ))
        t[is.infinite(t)] = 0
        mutations$CAF = t
        t = unlist(lapply(
          X = inf$SWAF,
          FUN = max,
          na.rm = T
        ))
        t[is.infinite(t)] = 0
        mutations$SWAF = t
        t = unlist(lapply(
          X = inf$TOPMED,
          FUN = max,
          na.rm = T
        ))
        t[is.infinite(t)] = 0
        mutations$TOPMED = t
      })
      
      # "info"" annotations also need to be added
      mutations = cbind(mutations, as.data.table(info(vcf))[, .(CSQ)])
      
      # Add Swegen counts
      snptable <- getSNPtable(PFconfig$snptable)
      mutations$Swegen_count = snptable[names(vcf), snptable$value]
      mutations$Swegen_count[is.na(mutations$Swegen_count)] = 0
      # Some have multiple rsIDs

      # Also double check Swegen in case pos in reference but rsID or not properly matched in data
      pos = paste0(mutations$chr,
                   ':',
                   mutations$start,
                   '_',
                   mutations$ref,
                   '/',
                   mutations$alt)
      snps = snptable[c(pos)][!is.na(value)] # extract from reference[those present]
      if (nrow(snps) > 0)
        mutations$Swegen_count[match(snps$name, pos)] = snps$value # put
      # # Some have multiple alt alleles
      # ix=grep(',',mutations$alt)
      # alts=data.table(ix,alt=mutations$alt[ix])
      # alts$counts[ix]=unlist(lapply(alts$pos,function(x) return(max(snptable[x,value],na.rm=T))))
      
      # More thorough filtering
      if (any(mutations$Swegen_count >= 10))
        mutations = mutations[-which(Swegen_count >= 10)]
      
      
      ## Annotate by cosmic (for possibly retaining non PASS hotspots, which is not necessarily smart)
      key = paste0(substr(mutations$chr, 4, 6),
                   ':',
                   mutations$start,
                   '-',
                   mutations$end)
      key = str_replace(key, 'X:', '23:')
      key = str_replace(key, 'Y:', '24:')
      
      #if (project != 'BTB')
      #  key = paste0(mutations$chr, ':', mutations$start, '-', mutations$end)
      cosmic_coding <- getCountTable(PFconfig$coding_table,"coding_table")
      cosmic_noncoding <- getCountTable(PFconfig$noncoding_table,"noncoding_table")
      counts = cbind(cosmic_coding[key, cosmic_coding$value], cosmic_noncoding[key, cosmic_noncoding$value], 0)
      max_ = apply(counts, 1, max, na.rm = T)
      mutations$Cosmic_count = max_
      
      mutations$rank_score = 0
      mutations$rank_terms = ''
      mutations$LOH = ''
      
      # Add Control Freec LOH.
      try({
        freec_loh <- getEnvVariable('freec_loh')
        if (nrow(mutations) > 0 & !is.null(freec_loh)) {
          for (i in 1:nrow(freec_loh)) {
            ix = freec_loh$chr[i] == mutations$chr &
              freec_loh$start[i] < mutations$start &
              freec_loh$end[i] > mutations$end
            ix = ix[!is.na(ix)]
            if (sum(ix) > 0)
              mutations$LOH[ix] = 'Y'
          }
        }
      }, silent = T)
      
      
      vep_header = strsplit(info(header(vcf))['CSQ', ][, 3], 'Format: ')[[1]][2]
      vep_header = strsplit(vep_header, '\\|')[[1]]
      
      # VEP annotation is put in annotation_table
      annotation_table = matrix(
        data = NA,
        nrow = length(unlist(mutations$CSQ)),
        ncol = length(vep_header) + 1
      )
      colnames(annotation_table) = c('ID', vep_header)
      row = 1
      for (i in 1:nrow(mutations)) {
        # for each variant
        # for each VEP annotation:
        for (j in 1:length(mutations$CSQ[[i]])) {
          line = strsplit(mutations$CSQ[[i]][j], '\\|')[[1]]
          annotation_table[row, 1] = mutations$ID[i]
          annotation_table[row, 1 + (1:length(line))] = line
          row = row + 1
        }
      }
      annotation_table = as.data.table(annotation_table)
      
      # Annotation table then concatenated
      haplotypecaller_table = rbind(haplotypecaller_table,
                                    merge(mutations, annotation_table, by = 'ID'))
      toc()
    } # done collecting from vcf files
    
    haplotypecaller_table$cumstart = haplotypecaller_table$start
    haplotypecaller_table$cumend = haplotypecaller_table$end
    # for each chr get cumulative pos
    tic("Get cumulative positions")
    for (i in 1:nrow(chrsz)) {
      ix = haplotypecaller_table$chr == chrsz$chr[i]
      haplotypecaller_table$cumstart[ix] = haplotypecaller_table$start[ix] +
        chrsz$starts[i]
      haplotypecaller_table$cumend[ix] = haplotypecaller_table$end[ix] +
        chrsz$starts[i]
    }
    setkey(haplotypecaller_table, 'sample')
    toc()
    
    # Due to size: remove intron/intergenic variants
    tic("Remove intron/intergenic")
    selection = haplotypecaller_table[!Consequence %in% c('intron_variant', 'intergenic_variant'), -c('CSQ')] # not needed after parsing/merging the annotations
    toc()
    
    tic("Calculate tiers")
    if (nrow(selection) > 0) {
      # Cosmic/local Tier2 :
      alltier1 <- getEnvVariable('alltier1')
      alltier2 <- getEnvVariable('alltier2')
      ix = selection$SYMBOL %in% alltier2
      selection$rank_score[ix] = 2
      selection$rank_terms[ix] = 'T2_gene'
      # tier 1 priority:
      ix = selection$SYMBOL %in% alltier1
      selection$rank_score[ix] = 2
      selection$rank_terms[ix] = 'T1_gene'
      
      # Add high impact
      ix = selection$IMPACT == 'HIGH'
      if (any(ix)) {
        selection$rank_score[ix] = selection$rank_score[ix] + 2
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'high_impact')
      }
      
      # Add moderate impact
      ix = selection$IMPACT == 'MODERATE'
      if (any(ix)) {
        selection$rank_score[ix] = selection$rank_score[ix] + 1
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'moderate_impact')
      }
      
      # additional +1 if high impact and TSG
      alltsg <- getEnvVariable('alltsg')
      ix = selection$IMPACT == 'HIGH' & selection$SYMBOL %in% alltsg
      if (any(ix)) {
        selection$rank_score[ix] = selection$rank_score[ix] + 1
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'high+TSG')
      }
      
      # Add clinvar pathogenic
      ix = grep('pathogenic', selection$CLIN_SIG)
      if (length(ix) > 0) {
        selection$rank_score[ix] = selection$rank_score[ix] + 2
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'clinvar')
      }
      
      # Add Polyphen/SIFT damaging/deleterious
      ix = union(grep('damaging', selection$PolyPhen),
                 grep('deleterious', selection$SIFT))
      if (length(ix) > 0) {
        selection$rank_score[ix] = selection$rank_score[ix] + 1
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'polyphen/SIFT')
      }
      
      # Add hotspots and near hotspots (within 2 residues of a hotspot)
      key = paste(
        selection$SYMBOL,
        str_replace(
          string = selection$Protein_position,
          pattern = '/.*',
          replacement = ''
        )
      )
      hotspots_snv <- getSNVhotspots(PFconfig$hotspots_snv) 
      hotspots_inframe <- getInframeHotspots(PFconfig$hotspots_inframe,hotspots_snv)  
      ix <-
        key %in% paste(hotspots_snv[, Hugo_Symbol], hotspots_snv[, Amino_Acid_Position]) |
        key %in% paste(hotspots_inframe[, Hugo_Symbol], hotspots_inframe[, Amino_Acid_Position])  ## Warning: Exact match used with inframes
      if (any(ix)) {
        selection$rank_score[ix] = selection$rank_score[ix] + 2
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'hotspot')
      }
      near_hotspots <- getEnvVariable('near_hotspots')
      ix = !ix &
        key %in% near_hotspots # the near_hotspot excludes those that were a hotspot
      if (any(ix)) {
        selection$rank_score[ix] = selection$rank_score[ix] + 1
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'near_hotspot')
      }
      
      # Add cosmic counts
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
      
      # Special case for TERT promoter
      ix = which(
        selection$SYMBOL == 'TERT' &
          selection$Consequence %in% c('5_prime_UTR_variant', 'upstream_gene_variant')
      )
      if (length(ix) > 0) {
        selection$rank_score[ix] = selection$rank_score[ix] + 4
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'TERT_upstr')
      }
      
      # Add TF binding variants near (100kb) cancer genes
      # 
      # this assigns more than one variables, so we have to ask for tumorgenes in the next line
      getTumorGenes(PFconfig$tumorgenes, PFconfig$local_tumorgenes)
      tumorgenes <- getEnvVariable('tumorgenes')
      alltumorgenes <- getEnvVariable('alltumorgenes')
      
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
            selection$rank_score[ix[i]] = selection$rank_score[ix[i]] +
              2
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
    }
    toc()
    
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
    
    haplotypecaller_selected <-
      selection[order(cumstart, Allele)][order(rank_score, decreasing = T)]
    haplotypecaller_selected <- haplotypecaller_selected[rank_score > getRankThreshold("haplotypecaller_threshold")]
    hc_csv_name <- paste0(csv_dir, '/', normal_sample, '_ranks.csv')
    fwrite(haplotypecaller_selected, file = hc_csv_name)
    cat("Results written to ",hc_csv_name)
  }
  assign(x='haplotypecaller_ids',value=haplotypecaller_ids,envir=pfenv)
  haplotypecaller_selected
}