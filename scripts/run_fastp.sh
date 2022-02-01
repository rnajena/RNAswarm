#!/bin/sh

READS_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/data/dadonaite_2019/reads
RESUTS_DIR=/home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019
THREADS=72

OUTPUT_DIR="$RESUTS_DIR/01-trimmed_reads"
mkdir $OUTPUT_DIR

FASTQC_DIR="$OUTPUT_DIR/fastqc"
mkdir $FASTQC_DIR

# Run fastp on dadonaite's fastq files
FASTQ_FILES=$(find $READS_DIR | grep .fastq | xargs ls)
for FASTQ_FILE in $FASTQ_FILES
do
    fastp -i $FASTQ_FILE -o $OUTPUT_DIR/$(basename $FASTQ_FILE .fastq)_trimmed.fastq\
          --failed_out $OUTPUT_DIR/$(basename $FASTQ_FILE .fastq)_failed_out.fastq\
          --json $OUTPUT_DIR/$(basename $FASTQ_FILE .fastq).json\
          --html $OUTPUT_DIR/$(basename $FASTQ_FILE .fastq).html\
          &>> $OUTPUT_DIR/$(basename $FASTQ_FILE .fastq).log
    
    fastqc -t $THREADS $OUTPUT_DIR/$(basename $FASTQ_FILE .fastq)_trimmed.fastq\
           -o $FASTQC_DIR &>> $FASTQC_DIR/fastqc.log
done