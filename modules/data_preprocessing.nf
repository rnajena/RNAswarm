/***********************************************************************
* fastp TRIMMING
***********************************************************************/
process fastpTrimming {
    label 'fastp'

    input:
    tuple val(sample_name), path(reads), val(group_name)

    output:
    tuple val(sample_name), path("${reads.baseName}_trimmed.fastq.gz"), val(group_name)

    publishDir "${params.output}/01-trimmed_reads", mode: 'copy'

    script:
    """
    fastp -i ${reads} -o ${reads.baseName}_trimmed.fastq.gz\
          --failed_out ${reads.baseName}_failed_out.fastq.gz\
          --json ${reads.baseName}.json\
          --html ${reads.baseName}.html\
    """
}

/***********************************************************************
* concatenate multifasta
***********************************************************************/
process concatenateFasta {
    label 'python3'

    input:
    tuple val(group_name), path(genome), val(is_genome_concatenated)

    output:
    tuple val(group_name), path("${genome.baseName}_concatenated.fasta"), val(true), path("${genome.baseName}_concatenated.csv")

    script:
    """
    fasta_preprocessor.py -c -i ${genome} -o ${genome.baseName}_concatenated.fasta
    """
}