#!/bin/sh
# Run fastqc on a repository full of fastq files
fastqc -t 32 -o /home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastqc/fastp\
                /home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/*/*.fastq\
            &>> /home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastqc/fastp/fastqc.log