input_dir = "."
csv_dir = "."
chrsz = data.table(
  chr = paste0('chr',c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                       "11", "12", "13", "14", "15", "16", "17", "18", 
                       "19", "20", "21", "22", "X", "Y")), 
  label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
            "12", "13", "14", "15", "16", "17", "18", "19", "20", 
            "21", "22", "X", "Y"),
  starts = c(0, 253902921, 500669562, 703689723, 898025604, 1084280925, 
             1259870526, 1424144247, 1574191848, 1717288929, 1855896810, 
             1995944331, 2134165212, 2253485013, 2363996694, 2470839615, 
             2565902256, 2654091777, 2739309138, 2802869979, 2871941700, 
             2923598421, 2979326382, 3140008743),
  length = c(248902921, 241766641, 198020161, 189335881, 181255321, 
             170589601, 159273721, 145047601, 138097081, 133607881, 
             135047521, 133220881, 114319801, 105511681, 101842921, 
             90062641, 83189521, 80217361, 58560841, 64071721, 46656721,
             50727961, 155682361, 56827081)
  )


reference_genome = 'GRCh38'
baf = "the_BAF_file_bugger"
logr = "the_logR_file_bugger"

tic("Shaping tumor genes dataframe")
tumorgenes=data.table::fread('~/reports/reference_data/cancer_gene_census.csv',key='Gene Symbol')
tumorgenes$chr=paste0('chr',str_replace(string=tumorgenes$`Genome Location`,pattern = ':.*',replacement = ''))
temp=str_replace(string=tumorgenes$`Genome Location`,pattern = '.*:',replacement = '')
tumorgenes$start=as.numeric(str_replace(string=temp,pattern = '-.*',replacement = ''))
tumorgenes$end=as.numeric(str_replace(string=temp,pattern = '.*-',replacement = ''))
tumorgenes$cumstart=NA; tumorgenes$cumend=NA
for (i in 1:nrow(chrsz)) {
  ix <- tumorgenes$chr==chrsz$chr[i]
  tumorgenes$cumstart[ix] <- tumorgenes$start[ix]+chrsz$starts[i]
  tumorgenes$cumend[ix] <- tumorgenes$end[ix]+chrsz$starts[i]
}
toc()
tic("Merging locals")
local_tumorgenes = data.table::fread('~/reports/reference_data/2018_gene_list_tere_ref.csv',key='Gene')
# get rid of rows without a gene name
local_tumorgenes = subset.data.frame(local_tumorgenes,local_tumorgenes$Gene!='')
# select only the gene names and the tier (a.k.a. score) 
local_tumorgenes = local_tumorgenes[,1:2]
alltumorgenes=unique(c(tumorgenes$`Gene Symbol`,local_tumorgenes$Gene))
alltsg=tumorgenes[grep('TSG',tumorgenes$`Role in Cancer`),tumorgenes$`Gene Symbol`]
toc()
tic("Getting fusions")
allfusion=tumorgenes[grep('fusion',`Role in Cancer`),`Gene Symbol`]
allfusionpairs=NULL
for (i in 1:length(allfusion)) {
  t=trimws(strsplit(tumorgenes[allfusion[i],`Translocation Partner`],', ')[[1]])
  if (length(t)>0) for (j in 1:length(t))
    allfusionpairs=c(allfusionpairs,paste(sort(c(allfusion[i],t[j])),collapse = ' '))
}
toc()
tic("Creating tiers")
# select tier 1 genes from both
# then get a list of gene names as the union of the two sets
tumorgenes_tier1 = subset.data.frame(tumorgenes,tumorgenes$Tier==1)
local_tumorgenes_tier1 = subset.data.frame(local_tumorgenes,local_tumorgenes$`Tier 1 and 2 for pediatric cancers final`==1)
# select the first column (gene names)
tumorgenes_tier1_genes = tumorgenes_tier1[,1] 
local_tumorgenes_tier1_genes = local_tumorgenes_tier1[,1]
alltier1 = union(tumorgenes_tier1_genes,local_tumorgenes_tier1_genes)

# DRY but whatever, do the same for tier 2 TODO make it nice without repetition
tumorgenes_tier2 = subset.data.frame(tumorgenes,tumorgenes$Tier==2)
local_tumorgenes_tier2 = subset.data.frame(local_tumorgenes,local_tumorgenes$`Tier 1 and 2 for pediatric cancers final`==2)
alltier2 = union(tumorgenes_tier2[,1],local_tumorgenes_tier2[,1])
#alltier1=union(tumorgenes[Tier==1,tumorgenes$`Gene Symbol`],local_tumorgenes[`Tier 1 and 2 for pediatric cancers final`==1,local_tumorgenes$Gene])
#alltier2=union(tumorgenes[Tier==2,tumorgenes$`Gene Symbol`],local_tumorgenes[`Tier 1 and 2 for pediatric cancers final`==2,local_tumorgenes$Gene])
toc()
tic("Reading COSMIC fusions")
cosmic_fusions=fread('~/reports/reference_data/cosmic_fusions_table.csv',key = 'name')
toc()

tic("Reading SNV hotspots")
hotspots_snv=unique(
  fread('~/reports/reference_data/hotspots_v2_snv.csv')[,.(Hugo_Symbol,Amino_Acid_Position)])[-grep('splice',Amino_Acid_Position)]
hotspots_snv$pos <- as.numeric(hotspots_snv$Amino_Acid_Position)
toc()

tic("Reading inframe hotspots")
hotspots_inframe=unique(fread('~/reports/reference_data/hotspots_v2_inframe.csv')[,.(Hugo_Symbol,Amino_Acid_Position)])
hotspots_inframe$start=as.numeric(str_replace(string = hotspots_inframe$Amino_Acid_Position,pattern = '-[0-9]*',replacement = ''))
hotspots_inframe$end=as.numeric(str_replace(string = hotspots_inframe$Amino_Acid_Position,pattern = '[0-9]*-',replacement = ''))
near_hotspots=NULL
for (i in 1:nrow(hotspots_inframe)) 
  near_hotspots=c(near_hotspots,
                  paste(hotspots_inframe$Hugo_Symbol[i],
                        c(seq(hotspots_inframe$start[i]-2,hotspots_inframe$start[i]+2),
                          hotspots_inframe$end[i]-2,hotspots_inframe$end[i]+2)))
for (i in 1:nrow(hotspots_snv))
  near_hotspots=c(near_hotspots,
                  paste(hotspots_snv$Hugo_Symbol[i],
                        seq(hotspots_snv$pos[i]-2,hotspots_snv$pos[i]+2)))
toc()

tic("Loading SweGen SNP table ... ")
# TODO: serialize somehow, for a large table it takes aaages
snptable = fread('~/reports/reference_data/swegen_snp_counts.small.csv',
                 key = 'name')
toc()

tic("Loading COSMIC coding table ... ")
cosmic_coding = fread('~/reports/reference_data/cosmic_coding_table.csv', key = 'name')
toc()

tic("Loading COSMIC non-coding table ... ")
cosmic_noncoding = fread('~/reports/reference_data/cosmic_noncoding_table.csv',
                         key = 'name')
toc()
