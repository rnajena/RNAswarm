#!/bin/sh

READS_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/data/dadonaite_2019/reads
OUTPUT_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019
THREADS=72

# Run fastqc on a repository full of fastq files
FASTQ_FILES=$(find $READS_DIR | grep .fastq | xargs ls)

for FASTQ_FILE in $FASTQ_FILES
do
    fastqc -t $THREADS $FASTQ_FILE -o $OUTPUT_DIR &>> $OUTPUT_DIR/fastqc.log
done