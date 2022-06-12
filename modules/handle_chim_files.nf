/************************************************************************
* handles .chim files
*************************************************************************/

process handleChimFiles {
    label 'python3'

    input:
    tuple val(name), path(genome), path(chim_file)

    output:
    tuple val(name), path("heatmaps_${chim_file.baseName}")

    publishDir "${params.output}/03-heatmaps/bwa-mem", mode: 'copy'

    script:
    """
    mkdir heatmaps_${chim_file.baseName}
    handle_chimeras.py ${genome} ${chim_file} heatmaps_${chim_file.baseName} --bwa_mode
    echo "mock change"
    """
}

process makeConfusionMatrix_bwa {
    label 'python3'

    input:
    tuple val(name), path(interaction_reads), path(genomic_reads), path(chim_file) path(bam_file)

    output:
    tuple val(name), path("${name}_bwa_confusion_matrix.txt")

    script:
    """
    art_templater.py -i ${interaction_reads} -g ${genomic_reads} -c ${chim_file} -b ${bam_file} > ${name}_bwa_confusion_matrix.txt
    """
}
