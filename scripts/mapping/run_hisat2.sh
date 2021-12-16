#!/bin/sh

GENOME_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/data/dadonaite_2019/genomes
OUTPUT_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/hisat2
READS_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/data/dadonaite_2019/reads

# Generate indexes for genomes of interest
mkdir -p $OUTPUT_DIR/indexes
for file in $GENOME_DIR/*.fasta; do
    echo $file
    hisat2-build $file $OUTPUT_DIR/indexes/$(basename $file .fasta)
done

# Run hisat2 with split read mapping enabled
hisat2 -p 64 -x $OUTPUT_DIR/indexes/pr8_dadonaite -1 $READS_DIR/pr8/SRR7350059.fastq -2 $READS_DIR/pr8/SRR9637504.fastq -S $OUTPUT_DIR/mappings/pr8/SRR7350059.sam