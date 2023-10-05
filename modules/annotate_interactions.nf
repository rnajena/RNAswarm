/*************************************************************************
* import annotation
*************************************************************************/
process annotateArrays {
    label 'RNAswarm'

    input:
    tuple val(sample_name), path(genome), path(sample_arrays)

    output:
    tuple val(sample_name), path(sample_arrays), path("${sample_name}_annotations")

    publishDir "${params.output}/06-annotations", mode: 'copy'

    script:
    """
    mkdir ${sample_name}_annotations
    annotate_interactions.py -d ${sample_arrays} -g ${genome} -o ${sample_name}_annotations
    """
}
