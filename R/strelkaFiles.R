strelkaFiles <- function() {
  # OK, this part has to be at a higher level at init
  result_files <- dir(recursive = T,full.names = T)
  # TODO DRY it out
  file_names_to_remove=unique(c(grep(pattern = '.png',x = result_files),
                                grep(pattern = 'Annotation.old',x = result_files),
                                grep(pattern = '/work/',x = result_files),
                                grep(pattern = 'VEP.CADD/Manta',x = result_files),
                                grep(pattern = 'VEP/haplo',x = result_files),
                                grep(pattern = 'VEP/mute',x = result_files),
                                grep(pattern = 'old.ControlFreec',x = result_files),
                                grep(pattern = 'old.BTB_ControlFreec',x = result_files)
  ))
  if (length(file_names_to_remove)>0) 
    result_files <- result_files[-file_names_to_remove]
  
  strelka_snv_file <- grep(pattern = ".*VEP/Strelka_.*_somatic_snvs.*vep.ann.vcf$",result_files,value = T)
  strelka_indel_file <- grep(pattern = ".*VEP/Strelka_.*_somatic_indels.*vep.ann.vcf$",result_files,value = T)

  # return with a file hash
  c("strelka_snv_file" = strelka_snv_file,
    "strelka_indel_file" = strelka_indel_file
  )
}