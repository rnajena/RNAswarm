#!/bin/sh

OUTPUT_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/hisat2
READS_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/data/dadonaite_2019/reads

# Generate indexes for genomes of interest
mkdir -p $OUTPUT_DIR/indexes
for file in $GENOME_DIR/*.fasta; do
    echo $file
    hisat2-build $file $OUTPUT_DIR/hisat2_indexes/$(basename $file .fasta)
done

# Run hisat2 with split read mapping enabled
hisat2 -p 64 -x $OUTPUT_DIR/indexes/pr8_dadonaite -U $READS_DIR/pr8/SRR7350059.fastq | samtools sort