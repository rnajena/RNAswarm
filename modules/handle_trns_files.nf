/************************************************************************
* handles .trns.txt files
*************************************************************************/

process handleTrnsFiles {
    label 'python3'
  
    input:
    tuple val(name), path(genome), path(trns_file), path(sam_file)

    output:
    tuple val(name), path("${sam_file.baseName}_plots"), path("${sam_file.baseName}_heatmaps.log")

    publishDir "${params.output}/03-heatmaps/segemehl", mode: 'copy'

    script:
    """
    mkdir ${sam_file.baseName}_plots
    handle_chimeras.py -g ${genome} -i ${trns_file} -o ${sam_file.baseName}_plots --segemehl_mode > ${sam_file.baseName}_heatmaps.log
    """
}
