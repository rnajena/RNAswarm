/*************************************************************************
* Fill arrays
*************************************************************************/

process fillArrays {
    label 'python3'

    input:
    tuple val(sample_name), path("${reads.baseName}.trns.txt"), path(genome), val(genome.baseName)

    output:
    tuple val(sample_name), path("${sample_name}_arrays"), path(genome), val(genome.baseName)

    publishDir "${params.output}/03-arrays", mode: 'copy'

    script:
    """
    mkdir ${sample_name}_arrays
    fill_aray.py ${trns_files} -g ${genome_file} -o ${sample_name}_arrays
    """
}
