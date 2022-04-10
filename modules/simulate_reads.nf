#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// filepaths
params.reads = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/trimmed_reads'
params.genomes = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/genomes'

/***********************************************************************
* ART simulator SIMULATE READS
***********************************************************************/
process inSilicoSeq_simulate_reads {
  label "simulate_reads"

  pus 8
  time '12h'
  executor 'slurm'
  conda '../envs/simulate_reads.yaml'

  input:
  tuple val(name), path(genome)

  output:
  tuple val(name), path(genome), path("${name}.fq")

  script:
  """
  art_illumina -ss HS25 -sam -i ${genome} -l 150 -f 10 -o ${name}
  """
}


 1) single-end read simulation
       art_illumina -ss HS25 -sam -i reference.fa -l 150 -f 10 -o single_dat

 2) paired-end read simulation
       art_illumina -ss HS25 -sam -i reference.fa -p -l 150 -f 20 -m 200 -s 10 -o paired_dat

 7) generate an extra SAM file with zero-sequencing errors for a paired-end read simulation
       art_illumina -ss HSXn -ef -i reference.fa -p -l 150 -f 20 -m 200 -s 10 -o paired_twosam_dat


workflow {
  genomes_ch = Channel
                .fromPath("${params.genomes}/*.fasta")
                .map{ file -> tuple(file.baseName, file) }.view()

  inSilicoSeq_simulate_reads( genomes_ch )
}