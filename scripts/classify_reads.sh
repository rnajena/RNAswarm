#!/bin/sh

INPUT_DIR="/data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/segemehl/mappings/"

# Classify reads using kraken2
for file in $(ls $INPUT_DIR*/*.fastq); do
    # Get the name of the file without the extension
    filename=$(basename "$file" .fastq)
    # Get the path of the subfolder the file is in
    dirname=$(dirname "$file")
    subdir=$(basename "$dirname")
    # Classify reads using kraken2
    kraken2 --db /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/kraken2/kraken2_db/kraken2_db --threads 24 --fastq-input --output $INPUT_DIR$subdir/$filename.kraken2.out $file
done