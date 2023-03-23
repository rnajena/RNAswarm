/*************************************************************************
* Plot heatmaps
*************************************************************************/

process plotHeatmaps {
  label 'python3'

  input:
  tuple val(genome_name), path(genome_file), path(trns_files)

  output:
  tuple val(genome_name), path(genome_name)

  publishDir "${params.output}/04-stats_and_plots/heatmaps", mode: 'copy'

  script:
  """
  mkdir ${genome_name}
  plot_heatmaps.py ${trns_files} -g ${genome_file} -o ${genome_name}
  """
}