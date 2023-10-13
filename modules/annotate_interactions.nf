/*************************************************************************
* annotate arrays de novo
*************************************************************************/
process annotateArrays {
    label 'RNAswarm'

    input:
    tuple val(sample_name), path(genome), path(sample_arrays)

    output:
    tuple val(sample_name), path(genome), path(sample_arrays), path("${sample_name}_annotations/${sample_name}_annotations.csv")

    publishDir "${params.output}/06-annotations", mode: 'copy'

    script:
    """
    mkdir ${sample_name}_annotations
    annotate_interactions.py -d ${sample_arrays} -g ${genome} -o ${sample_name}_annotations
    """
}

/*************************************************************************
* merge annotations
*************************************************************************/
process mergeAnnotations {
    label 'RNAswarm'

    input:
    path(annotations)

    output:
    path("merged_annotations.csv")

    publishDir "${params.output}/06-annotations" , mode: 'copy'

    script:
    """
    echo "id,segment01,start01,end01,segment02,start02,end02" >> merged_annotations.csv
    merge_annotation_tables.py ${annotations} -o merged_annotations.csv
    """
}