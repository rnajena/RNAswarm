/*************************************************************************
* segemehl INDEX
*************************************************************************/
process segemehlIndex {
    label 'segemehl'

    input:
    tuple val(group_name), path(genome)

    output:
    tuple val(group_name), path(genome), path("${group_name}.idx")

    script:
    """
    segemehl.x -d ${genome} -x ${group_name}.idx
    """
}

/*************************************************************************
* segemehl RUN
*************************************************************************/
process segemehl {
    label 'segemehl'

    input:
    tuple val(sample_name), path(reads), val(group_name), path(genome), path(index)

    output:
    tuple val(sample_name), path("${reads.simpleName}.trns.txt"), path("${reads.simpleName}.sngl.bed"), path("${reads.simpleName}.mult.bed"), path("${reads.simpleName}*_segemehl.sam"), val(group_name), path(genome)

    script:
    """
    segemehl.x -i ${index}\
               -d ${genome}\
               -q ${reads}\
               -S ${reads.baseName}\
               -t ${params.max_cpus}\
               > ${reads.simpleName}_segemehl.sam
    """
}

/*************************************************************************
* segemehl PUBLISH
*************************************************************************/
process segemehlPublish {
    label 'segemehl'
    publishDir "${params.output}/02-mappings/segemehl", mode: 'copy'

    input:
    tuple val(name), path(trns_file), path(sngl_file), path(mult_file), path(sam_file),  val(genome_name), path(genome)

    output:
    tuple val(name), path(trns_file), path(sngl_file), path(mult_file), path(sam_file)

    script:
    """
    ls
    """
}

/*************************************************************************
* samtools CONVERT SAM TO BAM
*************************************************************************/

process convertSAMtoBAM {
    label 'samtools'
    publishDir "${params.output}/02-mappings/${mapper}", mode: 'copy'

    input:
    tuple val(name), path(mappings), val(mapper)

    output:
    tuple val(name), path("${mappings.baseName}.bam")

    script:
    """
    samtools view -@ 8 -S -b ${mappings} > ${mappings.baseName}.bam
    """
}

/************************************************************************
* handles .trns.txt files
*************************************************************************/

process handleTrnsFiles {
    label 'RNAswarm'
    publishDir "${params.output}/03-heatmaps/segemehl", mode: 'copy'
  
    input:
    tuple val(name), path(genome), path(trns_file), path(sam_file)

    output:
    tuple val(name), path("${sam_file.baseName}_plots"), path("${sam_file.baseName}_heatmaps.log")

    script:
    """
    mkdir ${sam_file.baseName}_plots
    handle_chimeras.py -g ${genome} -i ${trns_file} -o ${sam_file.baseName}_plots --segemehl_mode > ${sam_file.baseName}_heatmaps.log
    """
}