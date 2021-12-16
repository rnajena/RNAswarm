#!/bin/sh

INPUT_DIR="/data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/segemehl/mappings/"

# Loop trough samfiles inside a repository
for file in $(ls $INPUT_DIR*/*.sam); do
    # Get the name of the file without the extension
    filename=$(basename "$file" .sam)
    # Get the path of the subfolder the file is in
    dirname=$(dirname "$file")
    subdir=$(basename "$dirname")
    # Generate bam files
    samtools sort -@ 24 $file > $INPUT_DIR$subdir/$filename.bam
    # Generate index files
    samtools index -@ 24 $INPUT_DIR$subdir/$filename.bam
done