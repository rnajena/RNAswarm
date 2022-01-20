#!/bin/sh

RESULTS_DIR="/home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019"
READS_DIR="/home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/01-trimmed_reads"
GENOME_DIR="/home/ru27wav/Projects/gl_iav-splash_freiburg/data/dadonaite_2019/genomes"
THREADS=72

mkdir $RESULTS_DIR/02-mappings
INDEXES_DIR="$RESULTS_DIR/02-mappings/segemehl_indexes"
MAPPINGS_DIR="$RESULTS_DIR/02-mappings/segemehl"

mkdir $MAPPINGS_DIR
mkdir $INDEXES_DIR

# Generate indexes for genomes of interest
for GENOME in $GENOME_DIR/*.fasta; do
    segemehl.x -x $INDEXES_DIR/$(basename $GENOME .fasta).idx -d $GENOME 2> $INDEXES_DIR/$(basename $GENOME .fasta).log
done

# Segemehl parameters
accuracy=9
minfragsco=15
minsplicecov=80
minfraglen=15
exclclipping=0

PARAMS="-A $accuracy -U $minfragsco -W $minsplicecov -Z $minfraglen "

# Run segemehl with split read mapping enabled
segemehl.x -i $INDEXES_DIR/pr8_dadonaite.idx\
           -d $GENOME_DIR/MDCK_concat-pr8_dadonaite.fasta\
           -q $READS_DIR/SRR7350059_trimmed.fastq\
           -t $THREADS -S $MAPPINGS_DIR/SRR7350059_concat\
           $PARAMS\
           > $MAPPINGS_DIR/SRR7350059_concat.sam\
           2> $MAPPINGS_DIR/SRR7350059_concat.log

segemehl.x -i $INDEXES_DIR/pr8_dadonaite.idx\
           -d $GENOME_DIR/MDCK_concat-pr8_dadonaite.fasta\
           -q $READS_DIR/SRR9637504_trimmed.fastq\
           -t $THREADS -S $MAPPINGS_DIR/SRR9637504_concat\
           $PARAMS\
           > $MAPPINGS_DIR/SRR9637504_concat.sam\
           2> $MAPPINGS_DIR/SRR9637504_concat.log

segemehl.x -i $INDEXES_DIR/udorn_dadonaite.idx\
           -d $GENOME_DIR/MDCK_concat-udorn_dadonaite.fasta\
           -q $READS_DIR/SRR7350060_trimmed.fastq\
           -t $THREADS -S $MAPPINGS_DIR/SRR7350060_concat\
           $PARAMS\
           > $MAPPINGS_DIR/SRR7350060_concat.sam\
           2> $MAPPINGS_DIR/SRR7350060_concat.log

segemehl.x -i $INDEXES_DIR/udorn_dadonaite.idx\
           -d $GENOME_DIR/MDCK_concat-udorn_dadonaite.fasta\
           -q $READS_DIR/SRR9637509_trimmed.fastq\
           -t $THREADS -S $MAPPINGS_DIR/SRR9637509_concat\
           $PARAMS\
           > $MAPPINGS_DIR/SRR9637509_concat.sam\
           2> $MAPPINGS_DIR/SRR9637509_concat.log

segemehl.x -i $INDEXES_DIR/wsn_dadonaite.idx\
           -d $GENOME_DIR/MDBK_concat-wsn_dadonaite.fasta\
           -q $READS_DIR/SRR6388155_trimmed.fastq\
           -t $THREADS -S $MAPPINGS_DIR/SRR6388155_concat\
           $PARAMS\
           > $MAPPINGS_DIR/SRR6388155_concat.sam\
           2> $MAPPINGS_DIR/SRR6388155_concat.log

segemehl.x -i $INDEXES_DIR/wsn_dadonaite.idx\
           -d $GENOME_DIR/MDBK_concat-wsn_dadonaite.fasta\
           -q $READS_DIR/SRR6388157_trimmed.fastq\
           -t $THREADS -S $MAPPINGS_DIR/SRR6388157_concat\
           $PARAMS\
           > $MAPPINGS_DIR/SRR6388157_concat.sam\
           2> $MAPPINGS_DIR/SRR6388157_concat.log