/************************************************************************
* handles .trns.txt files
*************************************************************************/

process handleTrnsFiles {
    label 'python3'
  
    input:
    tuple val(name), path(genome), path(trns_file), path(sam_file)

    output:
    tuple val(name), path("${sam_file.baseName}_heatmaps"), path("${sam_file.baseName}_heatmaps.log")

    publishDir "${params.output}/03-heatmaps/segemehl", mode: 'copy'

    script:
    """
    mkdir ${sam_file.baseName}_heatmaps
    handle_chimeras.py ${genome} ${trns_file} ${sam_file.baseName}_heatmaps --segemehl_mode > ${sam_file.baseName}_heatmaps.log
    echo "mock mokaka change"
    """
}


// process makeConfusionMatrix_trns {
//     label 'python3'

//     input:
//     tuple val(name), path(interaction_reads), path(genomic_reads), path(trns_file) path(bam_file)

//     output:
//     tuple val(name), path("${name}_bwa_confusion_matrix.txt")

//     script:
//     """
//     table_maker.py -i ${interaction_reads} -g ${genomic_reads} -t ${trns_file} -b ${bam_file} > ${name}_bwa_confusion_matrix.txt
//     """
// }
