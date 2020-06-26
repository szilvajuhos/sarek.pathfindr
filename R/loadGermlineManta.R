loadGermlineManta <- function(manta_files) {
  manta_normal_file <- manta_files["manta_normal_file"]
  swegen_manta_all <- manta_files["swegen_manta_all"]$swegen_manta_all # because it is a list
  manta_normal_table <- NULL
  
  alltier1 <- get("alltier1",pfenv)
  alltier2 <- get("alltier2",pfenv)
  cosmic_fusions <- getCountTable(PFconfig$fusions_table,"fusions_table")
  allfusionpairs <-get("allfusionpairs",pfenv)
  allfusion <-get("allfusion",pfenv)
  
  selected_manta <- NULL
  
  # first collect PASS ids from all samples
  allpass <- NULL
  vcfs <- list()
  cat("Reading Manta normal file ... \n")
  if (length(manta_normal_file) > 0)
    if (!is.na(manta_normal_file))
      # TODO: there is only a single manta normal file up to now, so we can make it simpler
      for (s in 1:length(manta_normal_file)) {
        file_normal <- as.character(manta_normal_file[s])
        vcfs[file_normal] = VariantAnnotation::readVcf(file = file_normal, genome = reference_genome)
        pass = DelayedArray::rowRanges(vcfs[[file_normal]])$FILTER == 'PASS'
        allpass = c(allpass, names(vcfs[[file_normal]])[pass])
      }
  # then collect variants...
  if (length(manta_normal_file[s]) > 0 & !is.na(manta_normal_file)) {
    for (s in 1:length(manta_normal_file)) { # OK if we are going to use more Manta normal files, right now we have only one
      file_normal = as.character(manta_normal_file[s])
      sample = strsplit(basename(file_normal), '[.]')[[1]][1]
      tic(paste("Collecting variants from",file_normal))
      vcf = vcfs[[file_normal]]
      vcf = vcf[names(vcf) %in% allpass]
      
      if (length(vcf) > 0) {
        g = geno(vcf)
        inf = info(vcf)
        rr = DelayedArray::rowRanges(vcf)
        
        tic("Calculate read counts")
        srtable = as.data.table(g$SR)
        colnames(srtable) = paste0('SplitReadSupport_', c('N')) # Verified on a test sample
        prtable = as.data.table(g$PR)
        colnames(prtable) = paste0('PairedReadSupport_', c('N'))
        toc
        
        tic("Calculate allele ratio")
        temp = data.table(
          pr.ref = sapply(prtable$PairedReadSupport_N, "[", 1),
          sr.ref = sapply(srtable$SplitReadSupport_N, "[", 1),
          pr.alt = sapply(prtable$PairedReadSupport_N, "[", 2),
          sr.alt = sapply(srtable$SplitReadSupport_N, "[", 2)
        )
        AFreq = apply(temp[, 3:4], 1, sum, na.rm = T) /
          (apply(temp[, 3:4], 1, sum, na.rm = T) + apply(temp[, 1:2], 1, sum, na.rm =
                                                           T))
        toc()
        tic("make table of variants")
        sv_table = cbind(
          data.table(
            ID = as.character(names(vcf)),
            sample,
            chr = as.character(seqnames(rr))
          ),
          as.data.table(rr)[, -1],
          AFreq,
          prtable,
          srtable,
          data.table(
            Swegen_count = rep(NA, length(vcf)),
            rank_score = 0,
            rank_terms = '',
            LOH = ''
          ),
          as.data.table(inf)
        )
        sv_table$cumstart = sv_table$start
        sv_table$cumend = sv_table$end
        sv_table$altchr = sv_table$chr
        sv_table$altpos = sv_table$END
        sv_table$plot = F
        sv_table$arc = -1
        toc()
        tic("Filter out variants seen 2+ (?) times in reference data")
        ## Key has only chr,start,end
        key = sv_table[, c('chr', 'start', 'end')]
        #if (project == 'BTB')
        #  key$chr = substr(key$chr, 4, 6)
        key$chr = substr(key$chr, 4, 6)
        key$imprecise = '(pr)'
        ## If imprecise, round the pos to 10
        ix = sv_table$IMPRECISE == T
        key$imprecise[ix] = '(impr)'
        key$start[ix] = round(key$start[ix] / 10) * 10
        key$end[ix] = round(key$end[ix] / 10) * 10
        key = paste(key$chr, key$start, key$end, key$imprecise)
        toc()
        tic("put in Swegen count")
        sv_table$Swegen_count = swegen_manta_all[key, swegen_manta_all$value]
        sv_table$Swegen_count[is.na(sv_table$Swegen_count)] = 0
        ## do the filter
        sv_table <- sv_table[Swegen_count < 2]
        toc()
      }
    }
    cat("Loop through all and extract endpoint chr and pos\n")
    tic("Loop through all and extract endpoint chr and pos")
    if (exists('sv_table') & nrow(sv_table) > 0) {
      # loop through all and extract endpoint chr and pos   <------ This one must be remade..
      for (i in 1:nrow(sv_table))
        try({
          # <----  sometimes error here, so try..
          t = strsplit(x = sv_table$ALT[[i]], split = ':')[[1]]
          if (length(t) > 1 & t[1] != "<DUP") {
            tchr = str_extract(t[1], '[0-9,X,Y]*$')
            sv_table$altchr[i] <- paste0('chr', tchr)
            # if (project == 'BTB')
            #   sv_table$altchr[i] <- paste0('chr', tchr)
            # else
            #   sv_table$altchr[i] <- tchr
            tt = str_extract(t[2], '^[0-9]*')
            sv_table$altpos[i] = as.numeric(tt)
          }
        }, silent = T)
      cat("for each chromosome get cumulative pos\n")
      sv_table$altcumpos = sv_table$altpos
      for (i in 1:nrow(chrsz)) {
        ix = sv_table$chr == chrsz$chr[i]
        if (sum(ix) > 0) {
          sv_table$cumstart[ix] = sv_table$start[ix] + chrsz$starts[i]
          sv_table$cumend[ix] = sv_table$end[ix] + chrsz$starts[i]
        }
        ix = sv_table$altchr == chrsz$chr[i]
        if (sum(ix) > 0) {
          sv_table$altcumpos[ix] = sv_table$altpos[ix] + chrsz$starts[i]
        }
      }
      cat("decide how it is to be plotted (not represented elsewhere, up or down arc)\n") 
      for (i in 1:nrow(sv_table)) {
        if (sv_table$chr[i] == sv_table$altchr[i]) {
          # intrachromosomal: plot always, with "positive" arc
          sv_table$plot[i] = T
          sv_table$arc[i] = 1
        } else if (sv_table$altcumpos[i] > sv_table$cumstart[i]) {
          # interchromosomal: plot if mate is to right (else mate will be plotted) with negative arc
          sv_table$plot[i] = T
          sv_table$arc[i] = -1
        }
      }
      
      cat("Add snpEff annotation in the ANN column.\n")
      h = strsplit(info(header(vcf))['ANN', ][, 3], 'annotations: \'')[[1]][2]
      snpEff_header = trimws(strsplit(h, '\\|')[[1]])
      
      # snpEff annotation is put in snpEff_table
      snpEff_table = matrix(
        data = NA,
        nrow = length(unlist(sv_table$ANN)),
        ncol = length(snpEff_header) + 1
      )
      colnames(snpEff_table) = c('ID', snpEff_header)
      row = 1
      for (i in 1:nrow(sv_table)) {
        # for each variant
        # for each VEP annotation:
        for (j in 1:length(sv_table$ANN[[i]]))
          if (length(sv_table$ANN[[i]]) > 0) {
            line = strsplit(sv_table$ANN[[i]][j], '\\|')[[1]]
            snpEff_table[row, 1] = sv_table$ID[i]
            snpEff_table[row, 1 + (1:length(line))] = line
            row = row + 1
          }
      }
      snpEff_table = unique(as.data.table(snpEff_table))[Annotation_Impact ==
                                                           'HIGH']  # <--- only HIGH impact kept for the normal
      
      cat("Filter out annotations of certain types where the ID has >N annotations of that type")
      ids = unique(snpEff_table$ID)
      for (id in ids) {
        common = c(
          'protein_protein_contact',
          'duplication',
          'structural_interaction_variant',
          'inversion',
          'transcript_ablation',
          'feature_ablation',
          'sequence_feature',
          'intergenic_region',
          'downstream_gene_variant',
          'upstream_gene_variant'
        )
        annotations = snpEff_table$Annotation[snpEff_table$ID == id] # the annotations of this variant
        table = table(annotations[annotations %in% common])
        table = table[table > 20] # the common annotations that appear >N times for this variant
        if (length(table) > 0) {
          remove = which(snpEff_table$ID == id &
                           snpEff_table$Annotation %in% names(table))
          snpEff_table = snpEff_table[-remove[-1], ] # saves one to make sure the variant has some annotation
        }
      }
      
      # cat("Add to data (all samples, one table)\n")
      manta_normal_table = rbind(manta_normal_table,
                                 merge(sv_table, snpEff_table, by = 'ID', all = F))
      setkey(manta_normal_table, 'sample')
      
    }  
    toc()
  }
  cat("done parsing each sample\n")
  
  cat("Prepare germline ranking...\n")
  
  # Prepare ranking and (some more) filtering
  if (!is.null(manta_normal_table) & nrow(manta_normal_table) > 0) {
    # ## Make table with most important info for ranking, and report
    # selected <- unique(manta_normal_table[,.(ID,sample,SVTYPE,chr,start,REF,ALT,AFreq,PairedReadSupport_T,SplitReadSupport_T,
    #                        Swegen_count,Rank_score='',Rank_terms='',Gene_Name,Annotation,Annotation_Impact)])
    selection = manta_normal_table[!is.na(Annotation)]
    
    if (nrow(selection) > 0) {
      tic("Increas for known cancer genes")
      # Known cancer genes affect ranking by +1
      for (gene in unique(selection$Gene_Name))
        if (!is.na(gene))
          if (gene != '') {
            ix = selection$Gene_Name == gene
            gene = sort(strsplit(gene, '&')[[1]])
            # cosmic/local Tier2:
            if (any(gene %in% alltier2)) {
              selection$rank_score[ix] = 2
              selection$rank_terms[ix] = 'T2_gene'
            }
            # tier 1 top priority:
            if (any(gene %in% alltier1)) {
              selection$rank_score[ix] = 2
              selection$rank_terms[ix] = 'T1_gene'
            }
          }
      toc()
      
      tic("Add high impact")
      ix = selection$Annotation_Impact == 'HIGH'
      if (any(ix)) {
        selection$rank_score[ix] = selection$rank_score[ix] + 2
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'high_impact')
      }
      toc()
      tic("Add moderate impact")
      ix = selection$Annotation_Impact == 'MODERATE'
      if (any(ix)) {
        selection$rank_score[ix] = selection$rank_score[ix] + 1
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'moderate_impact')
      }
      toc()
      tic("+1 for focal")
      ix = selection$chr == selection$altchr &
        selection$end - selection$start < 3e6 &
        selection$altpos - selection$start < 3e6
      if (any(ix)) {
        selection$rank_score[ix] = selection$rank_score[ix] + 1
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'focal')
      }
      toc()
      tic("cosmic_>xx, fusion_gene or just fusion")
      ix = grep('fusion', selection$Annotation)
      if (length(ix) > 0)
        for (i in ix)
          if (selection$Gene_Name[i] != '') {
            gene = sort(strsplit(selection$Gene_Name[i], '&')[[1]])
            if (paste(gene, collapse = ' ') %in% cosmic_fusions[cosmic_fusions$value > 5, cosmic_fusions$name]) {
              selection$rank_score[i] = selection$rank_score[i] + 3
              selection$rank_terms[i] = paste(selection$rank_terms[i], 'cosmic_>5')
            } else if (paste(gene, collapse = ' ') %in% cosmic_fusions$name) {
              selection$rank_score[i] = selection$rank_score[i] + 2
              selection$rank_terms[i] = paste(selection$rank_terms[i], 'cosmic_>1')
            }
            if (paste(gene, collapse = ' ') %in% allfusionpairs) {
              selection$rank_score[i] = selection$rank_score[i] + 2
              selection$rank_terms[i] = paste(selection$rank_terms[i], 'CGC_fusion')
            } else if (any(gene %in% allfusion)) {
              selection$rank_score[i] = selection$rank_score[i] + 1
              selection$rank_terms[i] = paste(selection$rank_terms[i], 'partial_CGC_fusion')
            } else {
              selection$rank_score[i] = selection$rank_score[i] + 0
              selection$rank_terms[i] = paste(selection$rank_terms[i], 'fusion')
            }
          }
      toc()
      tic("-1 for ablation if long del")
      ix = intersect(
        which(
          selection$SVTYPE == 'DEL' &
            selection$chr == selection$altchr &
            abs(selection$altpos - selection$start) > 3e6
        ),
        grep('ablation', selection$Annotation)
      )
      if (length(ix) > 0) {
        selection$rank_score[ix] = selection$rank_score[ix] - 1
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'long_del')
      }
      toc()
      tic("-1 for duplication if long dup")
      ix = intersect(
        which(
          selection$SVTYPE == 'DUP' &
            selection$chr == selection$altchr &
            abs(selection$altpos - selection$start) > 10e6
        ),
        grep('duplication', selection$Annotation)
      )
      if (length(ix) > 0) {
        selection$rank_score[ix] = selection$rank_score[ix] - 1
        selection$rank_terms[ix] = paste(selection$rank_terms[ix], 'long_dup')
      }
      toc()
      tic("Add Control Freec LOH")
      try({
        if (nrow(selection) > 0 &
            !is.null(freec_loh))
          for (i in 1:nrow(selection)) {
            ix = freec_loh$chr == selection$chr[i] &
              freec_loh$start < selection$start[i] &
              freec_loh$end > selection$end[i]
            ix = ix[!is.na(ix)]
            if (sum(ix) > 0)
              selection$LOH[i] = 'Y'
          }
      }, silent = T)
      toc()
      firstcols = c(
        'ID',
        'sample',
        'Gene_Name',
        'rank_score',
        'rank_terms',
        'LOH',
        'AFreq',
        'Annotation',
        'Annotation_Impact',
        'Swegen_count'
      )
      cols = colnames(selection)
      setcolorder(x = selection, neworder = c(firstcols, cols[!cols %in% firstcols]))
      
      manta_normal_selected <- selection[order(Feature_ID)][order(rank_score, decreasing = T)]
      manta_normal_file <- paste0(csv_dir, '/', sample, '_manta_normal.csv')
      tic(paste("Writing file ",manta_normal_file))
      # sometimes we have CSQ (VEP annotation) sometimes we don't
      ignoredCols = NULL
      if(length(grep('CSQ',colnames(manta_normal_selected),value=TRUE)) > 0) {
        ignoredCols = c('ANN','CSQ','Gene_ID')
      } else {
        cat("No VEP annotations, writing out without CSQ\n")
        ignoredCols = c('ANN','Gene_ID')
      }
      selected_manta <- manta_tumor_selected[,-ignoredCols][rank_score>3]  # original is rank_score > 3
      fwrite(selected_manta,file = manta_normal_file_name)
      toc()
    }
  }
  browser()
  selected_manta
}