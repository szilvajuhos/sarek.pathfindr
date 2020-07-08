conda activate sarek.pathfindr
time Rscript <(echo "rmarkdown::render('runPF.Rmd',knit_root_dir='.',output_file='bugger.html')")
