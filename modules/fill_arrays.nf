/*************************************************************************
* Fill arrays
*************************************************************************/

process fillArrays {
    label 'python3'

    input:
    tuple val(sample_name), path(trns_files), val(genome_name), path(genome)

    output:
    tuple val(sample_name), path(trns_files), val(genome_name), path(genome), path("${sample_name}_arrays")

    publishDir "${params.output}/03-arrays", mode: 'copy'

    script:
    """
    mkdir ${sample_name}_arrays
    fill_arrays.py ${trns_files} -g ${genome_file} -o ${sample_name}_arrays
    """
}
