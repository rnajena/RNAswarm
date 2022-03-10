nextflow.enable.dsl=2

params.bwa_mappings = '/home/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/mappings'
params.genomes = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/genomes'
params.heatmap_dir = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/heatmaps'

/************************************************************************
* handles .chim files
*************************************************************************/

process findChimeras {
  label 'handle_bwa_bam_files'

  cpus 8
  time '2h'
  executor 'slurm'
  conda '../envs/python2.yaml'

  input:
  tuple val(name), path(mapping)
  
  output:
  tuple val(name), path("${mapping.baseName}.chim")

  script:
  """
  python /home/ru27wav/Projects/gl_iav-splash_freiburg/src/RNAswarm/bin/find_chimeras_rs.py -i ${mapping} -o ${mapping.baseName}.chim
  """
}

/************************************************************************
* runs complete workflow for bwa mappings
*************************************************************************/

workflow {

  main:

    bwa_mappings_ch = Channel
                     .fromPath("${params.bwa_mappings}/*_bwa.bam")
                     .map{ file -> tuple(file.baseName[0..-27], file) }

    findChimeras( bwa_mappings_ch )
}
