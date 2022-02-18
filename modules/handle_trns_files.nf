nextflow.enable.dsl=2

params.mappings = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/segemehl_mappings'
params.genomes = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/genomes'
params.heatmap_dir = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/heatmaps'

process handleTrnsFiles {
  label 'segemehl'
  
  cpus 8
  executor 'slurm'
  conda '../envs/python3.yaml'

  input:
  tuple val(name), path(genome), path(mapping)

  output:
  tuple val(name), path("heatmaps_${mapping.baseName}")
  
  publishDir "${params.heatmap_dir}"

  script:
  """
  mkdir heatmaps_${mapping.baseName}
  python /home/ru27wav/Projects/gl_iav-splash_freiburg/src/RNAswarm/bin/handle_trns_files.py ${genome} ${mapping} heatmaps_${mapping.baseName}
  """
}

workflow {

    main:

        genomes_ch = Channel
                    .fromPath("$params.genomes/*.fasta")
                    .map{ file -> tuple(file.baseName, file) }
        
        mappings_ch = Channel
                .fromPath("$params.mappings/*.trns.txt")
                .map{ file -> tuple(file.baseName[0..-19], file) }

        handleTrnsFiles_input_ch = genomes_ch.combine(mappings_ch, by: 0)

        handleTrnsFiles( handleTrnsFiles_input_ch )
}