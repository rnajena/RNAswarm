#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/dadonaite_2019/reads'
params.trimmed_reads = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/trimmed_reads'

/***********************************************************************
* fasqc REPORT
***********************************************************************/
process fastqcReport {
  label 'preprocessing'

  cpus 8
  time '12h'
  executor 'slurm'
  conda '../envs/preprocessing_qc.yaml'

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), path("${reads.baseName}_fastqc")
  
  publishDir "${params.trimmed_reads}", mode: 'copy'

  script:
  """
  mkdir ${reads.baseName}_fastqc
  fastqc ${reads} -t 8 -o ${reads.baseName}_fastqc
  """
}

/***********************************************************************
* fastp TRIMMING
***********************************************************************/
process fastpTrimming {
  label 'preprocessing'

  cpus 8
  time '12h'
  executor 'slurm'
  conda '../envs/preprocessing_fastp.yaml'

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), path("${reads.baseName}_trimmed.fastq")
  
  publishDir "${params.trimmed_reads}", mode: 'copy'

  script:
  """
  fastp -i ${reads} -o ${reads.baseName}_trimmed.fastq\
        --failed_out ${reads.baseName}_failed_out.fastq\
        --json ${reads.baseName}.json\
        --html ${reads.baseName}.html\
  """
}

/***********************************************************************
* skewer TRIMMING
***********************************************************************/
process skewerTrimming {
  label 'preprocessing'

  cpus 8
  time '12h'
  executor 'slurm'
  conda '../envs/preprocessing_skewer.yaml'

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), path("${reads.baseName}_trimmed.fastq")
  
  publishDir "${params.trimmed_reads}", mode: 'copy'

  script:
  """
  skewer -o ${reads.baseName}_trimmed.fastq\
         ${reads}
  """
}

/************************************************************************
* runs complete preprocessing workflow
************************************************************************/
workflow {

  main:

    reads_ch = Channel
              .fromPath("${params.reads}/*/*.fastq")
              .map{ file -> tuple(file.baseName, file) }.view()
    
    fastpTrimming( reads_ch )

    qc_ch = fastpTrimming.out.concat(reads_ch).view()

    fastqcReport( qc_ch )
}