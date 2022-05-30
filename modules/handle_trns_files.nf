/************************************************************************
* handles .trns.txt files
*************************************************************************/

process handleTrnsFiles {
    label 'python3'
  
    input:
    tuple val(name), path(genome), path(trns_file), path(bam_file)

    output:
    tuple val(name), path("heatmaps_${bam_file.baseName}")

    publishDir "${params.output}/03-heatmaps/segmehl", mode: 'copy'

    script:
    """
    mkdir heatmaps_${bam_file.baseName}
    handle_chimeras.py ${genome} ${trns_file} heatmaps_${bam_file.baseName} --segemehl_mode
    """
}

process makeConfusionMatrix_trns{
    label 'python3'

    input:
    tuple val(name), path(interaction_reads), path(genomic_reads), path(trns_file) path(bam_file)

    output:
    tuple val(name), path("${name}_bwa_confusion_matrix.txt")

    script:
    """
    art_templater.py -i ${interaction_reads} -g ${genomic_reads} -t ${trns_file} -b ${bam_file} > ${name}_bwa_confusion_matrix.txt
    """
}
