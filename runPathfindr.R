# This script goes through pathfindr steps in the right order, and creates a 
# set of CSV files with ranks
# 
# firts, ASCAT and Control-FREEC: we want to use Control-FREEC LOH data later
# 

library(sarek.pathfindr)
library(tictoc) # to measure times
library(config)


browser()

tic("Scoring ASCAT")
score_ascat(PFconfig)
toc()

tic("Control-FREEC scores and LOH")
score_freec(PFconfig)
toc()

tic("Mutect2 scores")
scoreMutect2(PFconfig$mutect2file)
toc()

tic("Strelka SNVs and indels")
scoreStrelka(PFconfig)
toc()
 
tic("Germline GATK HaplotypeCaller calls")
scoreHaplotypeCaller(PFconfig)
toc()

tic("Structural variants by Manta")
scoreManta(PFconfig)
toc()