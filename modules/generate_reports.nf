// /***********************************************************************
// * make tables and plots with mapping statistics
// ***********************************************************************/

// process runTableMaker {
//   label "simulate_interactions"

//   input:
//   tuple val(name)

//   output:
//   tuple val(name)

//   publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

//   script:
//   """
//   art_templater.py -t <trans_file> -c <chim_file> -bs <segemehl_bam_file> -bb <bwa_bam_file> 
//   """
// }

/***********************************************************************
* fasqc REPORT
***********************************************************************/
process fastqcReport {
  label 'preprocessing'

  input:
  tuple val(name), path(reads)

  output:
  path("${reads.baseName}_fastqc")

  publishDir "${params.output}//04-stats_and_plots", mode: 'copy'

  script:
  """
  mkdir ${reads.baseName}_fastqc
  fastqc ${reads} -t 8 -o ${reads.baseName}_fastqc
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
  path("${mappings.baseName}.log")

  publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

  script:
  """
  samtools flagstats -@ ${params.cpus} ${mappings} > ${mappings.baseName}.log
  """
}

/*************************************************************************
* make coverage plots
*************************************************************************/

// process makeCoveragePlots {
//   label 'python3'

//   input:

// }

/*************************************************************************
* make Kraken2 database
*************************************************************************/

process makeKrakenDatabase {
  label 'generate_report_kraken'

  input:
  path(genomes)

  output:
  path('kraken_db')

  publishDir './assets/kraken2', mode: 'copy'

  script:
  """
  kraken2-build --use-ftp --download-taxonomy --db kraken_db
  for file in ${genomes}
  do
      kraken2-build --add-to-library \$file --db kraken_db
  done
  kraken2-build --threads ${params.max_cpus} --build --db kraken_db
  """
}

/*************************************************************************
* run Kraken2
*************************************************************************/

process runKraken {
  label 'generate_report_kraken'

  input:
  tuple val(name), path(reads), path(kraken_db)
  
  output:
  tuple val(name), path("${reads.baseName}_kraken.txt"), path("${reads.baseName}_kraken.out"), path("${reads.baseName}_kraken_classified.txt")

  publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

  script:
  """
  kraken2 --threads ${params.max_cpus}\
          --db ${kraken_db}\
          --report ${reads.baseName}_kraken.txt\
          --classified-out ${reads.baseName}_kraken_classified.txt\
          ${reads} > ${reads.baseName}_kraken.out
  """
}

/*************************************************************************
* Run MultiQC to generate an output HTML report
*************************************************************************/

process runMultiQC {
  label 'generate_report_multiqc'

  input:
  path(logs)

  output:
  path "multiqc_report.html"
  path "multiqc_data"

  publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

  script:
  """
  multiqc -d .
  """
}
