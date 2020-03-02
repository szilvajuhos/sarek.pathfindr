# control freec files
freec_files <-function(PFconfig) {
  result_files <- dir(path = PFconfig$freec_directory, recursive = F,full.names = T)
  # TODO DRY it out
  file_names_to_remove=unique(c(grep(pattern = '.png',x = result_files)))
  if (length(file_names_to_remove)>0) 
    result_files <- result_files[-file_names_to_remove]
  
  freec_Tratio_file <- result_files[grep(pattern = '[TR].hg38.pileup.gz_ratio.txt',result_files,perl = T)]
  freec_Nratio_file <- result_files[grep(pattern = '[BNT].hg38.pileup.gz_normal_ratio.txt',result_files,perl = T)][1]
  freec_Tbaf_file <- result_files[grep(pattern = '[TR].hg38.pileup.gz_BAF.txt',result_files,perl = T)]
  freec_Nbaf_file <- result_files[grep(pattern = '[BN].hg38.pileup.gz_BAF.txt',result_files,perl = T)][1]
  freec_info_file <- result_files[grep(pattern = 'hg38.pileup.gz_info.txt',result_files,perl = T)]
  freec_cnv_file <- result_files[grep(pattern = "^.*[TR]\\.hg38\\.pileup\\.gz_CNVs$",result_files,perl = T)]
  
  c( "freec_Tratio_file" = freec_Tratio_file,
     "freec_Nratio_file" = freec_Nratio_file,
     "freec_Tbaf_file" = freec_Tbaf_file,
     "freec_Nbaf_file" = freec_Nbaf_file,
     "freec_info_file" = freec_info_file,
     "freec_cnv_file" = freec_cnv_file
     )
}