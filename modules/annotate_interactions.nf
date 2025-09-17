

process normalizeArrays {
    label 'RNAswarm'

    input:
    tuple val(sample_name), path(genome), path(sample_arrays)

    output:
    tuple val(sample_name), path(genome), path("${sample_arrays}_normalized")

    script:
    """
    mkdir ${sample_arrays}_normalized
    normalize_arrays.py -d ${sample_arrays} -g ${genome} -o ${sample_arrays}_normalized
    """
}


/*************************************************************************
* annotate arrays de novo
*************************************************************************/
process annotateArrays {
    label 'RNAswarm'
    publishDir "${params.output}/06-annotations", mode: 'copy'

    input:
    tuple val(sample_name), path(genome), path(sample_arrays)

    output:
    tuple val(sample_name), path(genome), path(sample_arrays), path("${sample_name}_annotations/${sample_name}_annotations.csv"), path("${sample_name}_annotations/${sample_name}_annotations_gmms.pickle")

    script:
    """
    mkdir ${sample_name}_annotations
    annotate_interactions.py -d ${sample_arrays} -g ${genome} -o ${sample_name}_annotations -m ${params.min_components} -M ${params.max_components} --step_size ${params.step_size} --sigma ${params.sigma}
    """
}

/*************************************************************************
* merge annotations
*************************************************************************/
process mergeAnnotations {
    label 'RNAswarm'
    publishDir "${params.output}/06-annotations" , mode: 'copy'

    input:
    path(annotations)

    output:
    path("merged_annotations.tsv")

    script:
    """
    echo "id\tsegment01\tstart01\tend01\tsegment02\tstart02\tend02" >> merged_annotations.tsv
    merge_annotation_tables.py ${annotations} -o merged_annotations.tsv
    """
}

/*************************************************************************
* deduplicate annotations
*************************************************************************/
process deduplicateAnnotations {
    label 'RNAswarm'
    publishDir "${params.output}/06-annotations" , mode: 'copy'

    input:
    tuple val(group), path(annotation_table), path(count_table)

    output:
    tuple val(group), path("deduplicated_annotations/annotation_table_deduplicated.tsv"), path("deduplicated_annotations/count_table_deduplicated.tsv")

    script:
    """
    mkdir deduplicated_annotations
    deduplicate_annotations.py -a ${annotation_table} -c ${count_table} -o deduplicated_annotations
    """
}