/*************************************************************************
* Fill arrays
*************************************************************************/

process fillArrays {
    label 'python3'

    input:
    tuple val(sample_name), path(trns_files), val(group_name), path(genome)

    output:
    tuple val(sample_name), path(trns_files), val(group_name), path(genome), path("${sample_name}_arrays")

    publishDir "${params.output}/03-arrays", mode: 'copy'

    script:
    """
    mkdir ${sample_name}_arrays
    fill_arrays.py ${trns_files} -g ${genome_file} -o ${sample_name}_arrays
    """
}


/*************************************************************************
* Merge arrays
*************************************************************************/
process mergeArrays {
    label 'python3'

    input:
    tuple val(group_name), path(genome), path(arrays)

    output:
    tuple val(group_name), path(genome), path("${group_name}_merged_arrays")

    publishDir "${params.output}/04-merged-arrays", mode: 'copy'

    script:
    """
    mkdir ${group_name}_merged_arrays
    merge_arrays.py ${arrays} -g ${genome_file} -o ${group_name}_merged_arrays
    """
}
