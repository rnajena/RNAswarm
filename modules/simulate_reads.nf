#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// filepaths
params.reads = '../test_results'
params.genomes = '../test_data'
params.read_len = 150
params.fold_coverage = 72

/***********************************************************************
* ART simulator SIMULATE ILLUMINA READS
***********************************************************************/
process inSilicoSeq_simulate_reads {
  label "simulate_reads"

  cpus 8
  time '12h'
  executor 'slurm'
  conda '../envs/simulate_reads.yaml'

  input:
  tuple val(name), path(genome)

  output:
  tuple val(name), path(genome), path("${name}.fq")

  publishDir "${params.reads}", mode: 'copy'

  script:
  """
  art_illumina -i ${genome} -l ${params.read_len} -f ${params.fold_coverage} -o ${name}
  """
}

workflow {
  genomes_ch = Channel
                .fromPath("${params.genomes}/*.fasta")
                .map{ file -> tuple(file.baseName, file) }.view()

  inSilicoSeq_simulate_reads( genomes_ch )
}