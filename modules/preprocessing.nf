nextflow.enable.dsl=2

params.reads = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/reads'
params.trimmed_reads = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/trimmed_reads'

/***********************************************************************
* fastp TRIMMING
*************************************************************************/
process fastpTrimming {
  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), path("${reads.baseName}_trimmed.fastq")
  
  publishDir "${params.trimmed_reads}", mode: 'copy'

  script:
  """
  fastp -i ${reads} -o ${reads.baseName}_trimmed.fastq\
        --failed_out ${reads.baseName}_failed_out.fastq\
        --json ${reads.baseName}.json\
        --html ${reads.baseName}.html\
    
  fastqc -t ${reads.baseName}_trimmed.fastq\
         -o ${reads.baseName}_fastqc
  """
}

/************************************************************************
* runs complete preprocessing workflow
************************************************************************/
workflow {

  main:
    reads_ch = Channel
              .fromPath("$params.genomes/*.fasta")
              .map{ file -> tuple(file.baseName, file) }

    segemehl( segemehl_input_ch )
}