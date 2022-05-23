/***********************************************************************
* make tables and plots with mapping statistics
***********************************************************************/

process runTableMaker {
  label "simulate_interactions"

  input:
  tuple val(name), path(segemehl_mapping), path(trns_file), path(bwa_mapping), path(chim_file)

  output:
  tuple val(name), path("${name}_summary.csv")

  publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

  script:
  """
  art_templater.py -bs ${segemehl_mapping} -t ${trns_file} -bb ${bwa_mapping} -c ${chim_file} >> ${name}_summary.csv
  """
}

/*************************************************************************
* samtools stats
*************************************************************************/

process getStats {
  label 'mapping_samtools'

  input:
  tuple val(name), path(mappings)

  output:
  tuple val(name), path("${mappings.baseName}.log")

  publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

  script:
  """
  samtools stats -@ ${params.cpus} ${mappings} > ${mappings.baseName}.log
  """
}
