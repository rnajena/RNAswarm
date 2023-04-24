/*************************************************************************
* import annotation
*************************************************************************/
process annotateArrays {
    label 'python3'

    input:
    tuple val(sample_name), val(genome_name), path(genome), path(sample_arrays)

    output:
    tuple val(sample_name), path(sample_arrays), path("${sample_name}_annotations.tsv")

    publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

    script:
    """
    annotate_interactions.py ${sample_arrays} -g ${genome} -o ${sample_name}_annotations.tsv
    """
}

