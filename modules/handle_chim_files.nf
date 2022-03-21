#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.chim_files = '/home/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/mappings'
params.genomes = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/genomes'
params.heatmap_dir = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/bwa_heatmaps'

/************************************************************************
* handles .chim files
*************************************************************************/

process handleChimFiles {
  label 'handle_chimeras'
  
  cpus 8
  time '2h'
  executor 'slurm'
  queue 'b_standard,b_fat,s_standard,s_fat'
  conda '../envs/python3.yaml'

  input:
  tuple val(name), path(genome), path(mapping)

  output:
  tuple val(name), path("heatmaps_${mapping.baseName}")
  
  publishDir "${params.heatmap_dir}", mode: 'copy'

  script:
  """
  mkdir heatmaps_${mapping.baseName}
  python /beegfs/ru27wav/Projects/gl_iav-splash_freiburg/src/RNAswarm/bin/handle_chimeras.py ${genome} ${mapping} heatmaps_${mapping.baseName} --bwa_mode
  """
}

/************************************************************************
* runs complete workflow for bwa mappings
*************************************************************************/

workflow {

  main:

      genomes_ch = Channel
                  .fromPath("$params.genomes/*.fasta")
                  .map{ file -> tuple(file.baseName, file) }.view()
 
      mappings_ch = Channel
              .fromPath("$params.chim_files/*.chim")
              .map{ file -> tuple(file.baseName[0..-26], file) }.view()

      handleChimFiles_input_ch = genomes_ch.combine(mappings_ch, by: 0).view()

      handleChimFiles( handleChimFiles_input_ch )
}
