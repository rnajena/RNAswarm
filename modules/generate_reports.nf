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

  output:
  path('kraken_db')

  publishDir './assets/kraken2_db', mode: 'copy'

  script:
  """
  mkdir kraken_db
  cd kraken_db
  wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz
  gunzip k2_standard_20210517.tar.gz
  tar -xvf k2_standard_20210517.tar
  rm k2_standard_20210517.tar
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
  tuple path("${reads.baseName}_kraken.txt"), path("${reads.baseName}_kraken.out"), path("${reads.baseName}_kraken_classified.txt")

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
* make Kraken2 database
*************************************************************************/

process makeSortmernaDatabase {
  label 'generate_report_sortmerna'

  output:
  path('sortmerna_db')

  publishDir './assets/sortmerna_db', mode: 'copy'

  script:
  """
  wget https://github.com/biocore/sortmerna/archive/5458bf3774714275a561fdd4c4e6e7eccf22c769.zip
  unzip 5458bf3774714275a561fdd4c4e6e7eccf22c769.zip
  mv sortmerna-5458bf3774714275a561fdd4c4e6e7eccf22c769/data/rRNA_databases sortmerna_db
  rm -rf sortmerna_db/README.txt sortmerna_db/scripts sortmerna_db/silva_ids_acc_tax.tar.gz
  """
}

/*************************************************************************
* run sortmerna
*************************************************************************/

process runSortmerna {
  label 'generate_report_sortmerna'

  input:
  tuple val(name), path(reads), path(sortmerna_db)

  output:
  tuple val(name), path(reads) //output has to be fixed so that multiQC can work

  publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

  script:
  """
  mkdir workdir
  sortmerna --ref sortmerna_db/rfam-5.8s-database-id98.fasta\
            --ref sortmerna_db/rfam-5s-database-id98.fasta\
            --ref sortmerna_db/silva-arc-16s-id95.fasta\
            --ref sortmerna_db/silva-arc-23s-id98.fasta\
            --ref sortmerna_db/silva-bac-16s-id90.fasta\
            --ref sortmerna_db/silva-bac-23s-id98.fasta\
            --ref sortmerna_db/silva-euk-18s-id95.fasta\
            --ref sortmerna_db/silva-euk-28s-id98.fasta\
            --workdir workdir\
            --reads ${reads}
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
