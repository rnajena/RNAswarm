#!/bin/sh
# Run fastp on a dadonaite's fastq files
fastp -i /data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/reads/pr8/SRR7350059.fastq\
      -o /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/pr8/SRR7350059_trimmed.fastq\
      --failed_out /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/pr8/SRR7350059_failed_out.fastq\
      --html /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/pr8/SRR7350059_report.html\
      &>> /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/pr8/SRR7350059_fastp.log

fastp -i /data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/reads/pr8/SRR9637504.fastq\
      -o /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/pr8/SRR9637504_trimmed.fastq\
      --failed_out /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/pr8/SRR9637504_failed_out.fastq\
      --html /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/pr8/SRR9637504_report.html\
      &>> /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/pr8/SRR9637504_fastp.log

fastp -i /data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/reads/udorn/SRR7350060.fastq\
        -o /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/udorn/SRR7350060_trimmed.fastq\
        --failed_out /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/udorn/SRR7350060_failed_out.fastq\
        --html /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/udorn/SRR7350060_report.html\
        &>> /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/udorn/SRR7350060_fastp.log

fastp -i /data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/reads/udorn/SRR9637509.fastq\
        -o /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/udorn/SRR9637509_trimmed.fastq\
        --failed_out /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/udorn/SRR9637509_failed_out.fastq\
        --html /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/udorn/SRR9637509_report.html\
        &>> /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/udorn/SRR9637509_fastp.log

fastp -i /data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/reads/wsn/SRR6388155.fastq\
        -o /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/wsn/SRR6388155_trimmed.fastq\
        --failed_out /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/wsn/SRR6388155_failed_out.fastq\
        --html /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/wsn/SRR6388155_report.html\
        &>> /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/wsn/SRR6388155_fastp.log

fastp -i /data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/reads/wsn/SRR6388157.fastq\
        -o /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/wsn/SRR6388157_trimmed.fastq\
        --failed_out /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/wsn/SRR6388157_failed_out.fastq\
        --html /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/wsn/SRR6388157_report.html\
        &>> /data/dessertlocal/projects/gl_iav-splash_freiburg/results/dadonaite_2019/fastp/wsn/SRR6388157_fastp.log
