/***********************************************************************
* heatmaps PLOT
***********************************************************************/
process plotHeatmaps {
  label 'python3'

  input:
  // the input is a genome file and a list of trns files
  tuple path(genome), path(trns_files)

  output:
  path "${genome.baseName}_heatmaps"

  publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

  script:
  """
  plot_heatmaps.py ${trns_files} -g ${genome} -o ${genome.baseName}_heatmaps
  """
}