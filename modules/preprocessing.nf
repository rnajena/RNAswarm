nextflow.enable.dsl=2

params.reads = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/reads'
params.trimmed_reads = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/trimmed_reads'

/***********************************************************************
* fasqc REPORT
***********************************************************************/
process fastqcReport {
  input:
  tuple val(name), path(reads), val(state)

  output:
  tuple val(name), path("${reads.baseName}_fastqc_${state}")
  
  publishDir "${params.trimmed_reads}/${reads.baseName}_fastqc_${state}", mode: 'copy'

  script:
  """
  mkdir ${reads.baseName}_fastqc_${state}"
  fastqc -t ${reads} -o ${reads.baseName}_fastqc_${state}"
  """
}

/***********************************************************************
* fastp TRIMMING
***********************************************************************/
process fastpTrimming {
  input:
  tuple val(name), path(reads), val(state)

  output:
  tuple val(name), path("${reads.baseName}_trimmed.fastq"), val('trimmed')
  
  publishDir "${params.trimmed_reads}", mode: 'copy'

  script:
  """
  fastp -i ${reads} -o ${reads.baseName}_trimmed.fastq\
        --failed_out ${reads.baseName}_failed_out.fastq\
        --json ${reads.baseName}.json\
        --html ${reads.baseName}.html\
  """
}

/************************************************************************
* runs complete preprocessing workflow
************************************************************************/
workflow {

  main:
    reads_ch = Channel
              .fromPath("$params.genomes/*.fasta")
              .map{ file -> tuple(file.baseName, file, 'untrimmed') }

    fastqcReport( reads_ch )
    fastpTrimming( reads_ch )
    fastqcReport( fastpTrimming.out )
}