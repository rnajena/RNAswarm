/***********************************************************************
* fasqc REPORT
***********************************************************************/
process fastqcReport {
  label 'preprocessing'

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), path("${reads.baseName}_fastqc")
  
  publishDir "${params.output}/01-trimmed_reads", mode: 'copy'

  script:
  """
  mkdir ${reads.baseName}_fastqc
  fastqc ${reads} -t 8 -o ${reads.baseName}_fastqc
  """
}

/***********************************************************************
* fastp TRIMMING
***********************************************************************/
process fastpTrimming {
  label 'preprocessing'

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), path("${reads.baseName}_trimmed.fastq")
  
  publishDir "${params.output}/01-trimmed_reads", mode: 'copy'

  script:
  """
  fastp -i ${reads} -o ${reads.baseName}_trimmed.fastq\
        --failed_out ${reads.baseName}_failed_out.fastq\
        --json ${reads.baseName}.json\
        --html ${reads.baseName}.html\
  """
}
