/************************************************************************
* handles .trns.txt files
*************************************************************************/

process handleTrnsFiles {
    label 'python3'
  
    input:
    tuple val(name), path(genome), path(trns_file), path(bam_file)

    output:
    tuple val(name), path("heatmaps_${trns_file.baseName}")

    publishDir "${params.output}/03-heatmaps/segmehl", mode: 'copy'

    script:
    """
    mkdir heatmaps_${trns_file.baseName}
    handle_chimeras.py ${genome} ${trns_file} heatmaps_${trns_file.baseName} --segemehl_mode
    """
}


