mantaFiles <- function() {
  # OK, this part has to be at a higher level at init
  result_files <- dir(recursive = T,full.names = T)
  # TODO DRY it out
  file_names_to_remove=unique(c(grep(pattern = '.png',x = result_files),
                                grep(pattern = 'Annotation.old',x = result_files),
                                grep(pattern = '/work/',x = result_files),
                                grep(pattern = 'VEP.CADD/Manta',x = result_files),
                                grep(pattern = 'VEP/haplo',x = result_files),
                                grep(pattern = 'VEP/mute',x = result_files),
                                grep(pattern = 'VEP/Strelk',x = result_files),
                                grep(pattern = 'old.ControlFreec',x = result_files),
                                grep(pattern = 'old.BTB_ControlFreec',x = result_files)
  ))
  if (length(file_names_to_remove)>0) 
    result_files <- result_files[-file_names_to_remove]
  
  # Manta structural variant files
  manta_tumor_file <- grep(pattern = ".*Eff/Manta_.*vs.*somaticSV.vcf.snpEff.ann.vcf$",result_files,value = T)
  manta_normal_file <- grep(pattern = ".*Eff/Manta_.*vs.*diploidSV.vcf.snpEff.ann.vcf$",result_files,value = T)[1]
  # return with a file hash
  c("manta_tumor_file" = manta_tumor_file,
    "manta_normal_file" = manta_normal_file,
    "swegen_manta_all" = '~/reports/reference_data/swegen_sv_counts.csv'
    )
}