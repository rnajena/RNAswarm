/***********************************************************************
* fastp TRIMMING
***********************************************************************/
process fastpTrimming {
    label 'fastp'
    publishDir "${params.output}/01-trimmed_reads", mode: 'copy'

    input:
    tuple val(sample_name), path(reads), val(group_name)

    output:
    tuple val(sample_name), path("${reads.baseName}_trimmed.fastq"), val(group_name)

    script:
    """
    fastp -i ${reads} -o ${reads.baseName}_trimmed.fastq\
          --failed_out ${reads.baseName}_failed_out.fastq\
          --json ${reads.baseName}.json\
          --html ${reads.baseName}.html\
    """
}

/***********************************************************************
* concatenate multifasta
***********************************************************************/
process concatenateFasta {
    label 'RNAswarm'

    input:
    tuple val(group_name), path(genome), val(is_genome_concatenated)

    output:
    tuple val(group_name), path("${genome.baseName}_concatenated.fasta"), val(true), path("${genome.baseName}_concatenated.csv")

    script:
    """
    fasta_preprocessor.py -c -i ${genome} -o ${genome.baseName}_concatenated.fasta
    """
}