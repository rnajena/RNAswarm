/*************************************************************************
* Fill arrays
*************************************************************************/

process fillArrays {
    label 'RNAswarm'

    input:
    tuple val(sample_name), path(trns_files), val(group_name), path(genome)

    output:
    tuple val(sample_name), path(trns_files), val(group_name), path(genome), path("${sample_name}_arrays")

    publishDir "${params.output}/03-arrays", mode: 'copy'

    script:
    """
    mkdir ${sample_name}_arrays
    fill_arrays.py ${trns_files} -g ${genome} -o ${sample_name}_arrays
    echo ${sample_name}_arrayss
    """
}

/*************************************************************************
* Fill arrays with ChimericFragments files
*************************************************************************/

process fillArraysCF {
    label 'RNAswarm'

    input:
    tuple val(group_name), path(genome), path(chimeric_fragments)

    output:
    tuple val(group_name), path(genome), path("${group_name}_arrays_CF")

    publishDir "${params.output}/03-arrays", mode: 'copy'

    script:
    """
    mkdir ${group_name}_arrays_CF
    fill_arrays_cf.py ${chimeric_fragments} -g ${genome} -o ${group_name}_arrays_CF
    echo ${group_name}_arrays_CF
    """
}


/*************************************************************************
* Merge arrays
*************************************************************************/
process mergeArrays {
    label 'RNAswarm'

    input:
    tuple val(group_name), path(genome), path(arrays)

    output:
    tuple val(group_name), path(genome), path("${group_name}_merged_arrays")

    publishDir "${params.output}/04-merged-arrays", mode: 'copy'

    script:
    """
    mkdir ${group_name}_merged_arrays
    merge_arrays.py ${arrays} -g ${genome} -o ${group_name}_merged_arrays
    """
}
