/************************************************************************
* handles .chim files
*************************************************************************/

process handleChimFiles {
    label 'python3'

    input:
    tuple val(name), path(genome), path(mapping)

    output:
    tuple val(name), path("heatmaps_${mapping.baseName}")

    publishDir "${params.output}/03-heatmaps/bwa-mem", mode: 'copy'

    script:
    """
    mkdir heatmaps_${mapping.baseName}
    handle_chimeras.py ${genome} ${mapping} heatmaps_${mapping.baseName} --bwa_mode
    """
}
