#!/bin/bash

usage() {                                      # Function: Print a help message.
  echo "Usage: $0 [ -s SUBJECT ] [ -a AUTHOR ] [-c CONFIG]" 1>&2 
}

while getopts ":s:a:c:" options; do
  case "${options}" in 
    s) SUBJECT=${OPTARG} ;;
    a) AUTHOR=${OPTARG} ;;
    c) CONFIG=${OPTARG} ;;
    :) echo "Error: -${OPTARG} requires an argument."; usage ; exit 1;;
    \? ) echo "Unknown option: -$OPTARG" >&2; usage; exit 1;;
    *) usage;exit 1;;
  esac
done

if [ -z ${SUBJECT+x} ]; then echo "SUBJECT, AUTHOR and CONFIG parameters are mandatory"; usage; exit 1; fi
if [ -z ${AUTHOR+x} ]; then echo "SUBJECT, AUTHOR and CONFIG parameters are mandatory"; usage; exit 1; fi
if [ -z ${CONFIG+x} ]; then echo "SUBJECT, AUTHOR and CONFIG parameters are mandatory"; usage; exit 1; fi

#time R -e "pftitle='$SUBJECT'; author='$AUTHOR'; rmarkdown::render('CLItest.Rmd',knit_root_dir='.',output_file='$SUBJECT.html',params=list(pfconfig_file='$CONFIG', subject='$SUBJECT'))"
time R -e "pftitle='$SUBJECT'; author= '$AUTHOR'; rmarkdown::render('runPF.Rmd',  knit_root_dir='.',output_file='$SUBJECT.html',params=list(pfconfig_file='$CONFIG', subject='$SUBJECT'))"

