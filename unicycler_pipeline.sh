#!/bin/bash
# Accessing folders and subfolders to retrieve information from files
# that contain DNA sequence read and run unicycler
#
# Accessing folders from input 
for folder in $(ls ./Input)
do
    # Accessing subfolders from input
    for file in $(ls ./Input/$folder)
    do
        # Getting the names of arguments for unicycler
        extention=${file##*-}
        if [ "1.fastq" = $extention ]
        then
            forward=$file
        elif [ "2.fastq" = $extention ]
        then
            reverse=$file
        elif [ "L1000.fastq.gz" = $extention ]
        then
            long=$file
        fi
    done
    # Creating directory to save output of unicycler
    isolate=$folder
    mkdir ./Output/$isolate
    # Running unicycler
    echo "Running unicycler with files from folder: $folder"
    unicycler -1 ./Input/$folder/$forward -2 ./Input/$folder/$reverse -l ./Input/$folder/$long -t 32 -o ./Output/$isolate
    echo "Done"
done
