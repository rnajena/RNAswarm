/***********************************************************************
* make tables and plots with mapping statistics
***********************************************************************/

process runTableMaker {
  label "simulate_interactions"

  input:
  tuple val(name)

  output:
  tuple val(name)

  publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

  script:
  """
  art_templater.py -t <trans_file> -c <chim_file> -bs <segemehl_bam_file> -bb <bwa_bam_file> 
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
  samtools flagstats -@ ${params.cpus} ${mappings} > ${mappings.baseName}.log
  """
}
