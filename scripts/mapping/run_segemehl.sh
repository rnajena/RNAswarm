#!/bin/sh

GENOME_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/data/dadonaite_2019/genomes
OUTPUT_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/segemehl
READS_DIR=/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp

# Generate indexes for genomes of interest
mkdir -p $OUTPUT_DIR/indexes
for file in $GENOME_DIR/*.fasta; do
    echo $file
    segemehl.x -x $OUTPUT_DIR/indexes/$(basename $file .fasta).idx -d $file 2> $OUTPUT_DIR/indexes/$(basename $file .fasta).log
done

# Run segemehl with split read mapping enabled
segemehl.x -i $OUTPUT_DIR/indexes/pr8_dadonaite.idx\
           -d $GENOME_DIR/pr8_dadonaite.fasta\
           -q $READS_DIR/pr8/SRR7350059_trimmed.fastq\
           -t 64 -S $OUTPUT_DIR/mappings/pr8/SRR7350059\
           > $OUTPUT_DIR/mappings/pr8/SRR7350059.sam\
           2> $OUTPUT_DIR/mappings/pr8/SRR7350059.log

segemehl.x -i $OUTPUT_DIR/indexes/pr8_dadonaite.idx\
           -d $GENOME_DIR/pr8_dadonaite.fasta\
           -q $READS_DIR/pr8/SRR9637504_trimmed.fastq\
           -t 64 -S $OUTPUT_DIR/mappings/pr8/SRR9637504\
            > $OUTPUT_DIR/mappings/pr8/SRR9637504.sam\
            2> $OUTPUT_DIR/mappings/pr8/SRR9637504.log

segemehl.x -i $OUTPUT_DIR/indexes/udorn_dadonaite.idx\
           -d $GENOME_DIR/udorn_dadonaite.fasta\
           -q $READS_DIR/udorn/SRR7350060_trimmed.fastq\
           -t 64 -S $OUTPUT_DIR/mappings/udorn/SRR7350060\
           > $OUTPUT_DIR/mappings/udorn/SRR7350060.sam\
           2> $OUTPUT_DIR/mappings/udorn/SRR7350060.log

segemehl.x -i $OUTPUT_DIR/indexes/udorn_dadonaite.idx\
           -d $GENOME_DIR/udorn_dadonaite.fasta\
           -q $READS_DIR/udorn/SRR9637509_trimmed.fastq\
           -t 64 -S $OUTPUT_DIR/mappings/udorn/SRR9637509\
           > $OUTPUT_DIR/mappings/udorn/SRR9637509.sam\
           2> $OUTPUT_DIR/mappings/udorn/SRR9637509.log

segemehl.x -i $OUTPUT_DIR/indexes/wsn_dadonaite.idx\
           -d $GENOME_DIR/wsn_dadonaite.fasta\
           -q $READS_DIR/wsn/SRR6388155_trimmed.fastq\
           -t 64 -S $OUTPUT_DIR/mappings/wsn/SRR6388155\
           > $OUTPUT_DIR/mappings/wsn/SRR6388155.sam\
           2> $OUTPUT_DIR/mappings/wsn/SRR6388155.log

segemehl.x -i $OUTPUT_DIR/indexes/wsn_dadonaite.idx\
           -d $GENOME_DIR/wsn_dadonaite.fasta\
           -q $READS_DIR/wsn/SRR6388157_trimmed.fastq\
           -t 64 -S $OUTPUT_DIR/mappings/wsn/SRR6388157\
           > $OUTPUT_DIR/mappings/wsn/SRR6388157.sam\
           2> $OUTPUT_DIR/mappings/wsn/SRR6388157.log