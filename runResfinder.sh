#!/bin/bash
# Run resfinder with multple fasta files contained in a folder
#
#################################################################
# NOTE: be careful with the different paths provided in the code#
#################################################################
#
# Accessing fasta files in input folder
for fastaFile in $(ls ~/Desktop/cge_input)
do
  # Save isolate ID number to make folder for results
  folderName=${fastaFile%%.*}
  mkdir ~/Desktop/cge_output/$folderName
  # run resfinder
  python3 ~/Documents/CGE/resfinder/run_resfinder.py \
    -o ~/Desktop/cge_output/$folderName -s "Escherichia coli" -l 0.6 -t 0.8 \
    --acquired --point -ifa ~/Desktop/cge_input/$fastaFile
  echo "Done $fastaFile"
done

echo "Work is done!"
