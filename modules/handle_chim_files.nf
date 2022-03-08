nextflow.enable.dsl=2

params.bwa_mappings = ''
params.genomes = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/genomes'
params.heatmap_dir = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/heatmaps'

/************************************************************************
* handles .chim files
*************************************************************************/

process handleBwaBamFiles {
  label 'handle_bwa_bam_files'

  cpus 8
  time '2h'
  executor 'slurm'
  conda '../envs/python2.yaml'

  input:
  tuple val(name), path(mapping)
  
  output:
  tuple val(name), path("heatmaps_${mapping.baseName}")

  script:
  """
  python find_chimeras_rs.py -i ${mapping} -o ${mapping.baseName}.chim
  """
}

/************************************************************************
* runs complete workflow for bwa mappings
*************************************************************************/

workflow {

  main:

      genomes_ch = Channel
                  .fromPath("$params.genomes/*.fasta")
                  .map{ file -> tuple(file.baseName, file) }
 
      mappings_ch = Channel
              .fromPath("$params.mappings/*.trns.txt")
              .map{ file -> tuple(file.baseName[0..-27], file) }

      handleTrnsFiles_input_ch = genomes_ch.combine(mappings_ch, by: 0)

      handleTrnsFiles( handleTrnsFiles_input_ch )
}
