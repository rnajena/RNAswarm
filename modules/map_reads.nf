/*************************************************************************
* segemehl INDEX
*************************************************************************/
process segemehlIndex {
    label 'mapping_segemehl'

    input:
    tuple val(name), path(genome)

    output:
    tuple val(name), path(genome), path("${name}.idx")

    script:
    """
    segemehl.x -d ${genome} -x ${name}.idx
    """
}

/*************************************************************************
* segemehl RUN
*************************************************************************/
process segemehl {
    label 'mapping_segemehl'

    input:
    tuple val(sample_name), path(genome), path(index), path(reads)

    output:
    tuple val(sample_name), path(genome), path("${reads.baseName}.trns.txt"), path("${reads.baseName}*_segemehl.sam")

    script:
    """
    segemehl.x -i ${index}\
               -d ${genome}\
               -q ${reads}\
               -S ${reads.baseName}\
               -t ${params.max_cpus}\
               > ${reads.baseName}_segemehl.sam
    """
}

/*************************************************************************
* segemehl PUBLISH
*************************************************************************/
process segemehlPublish {
    label 'mapping_segemehl'

    input:
    tuple val(name), path(trns_file), path(sam_file)

    output:
    tuple val(name), path(trns_file)

    publishDir "${params.output}/02-mappings/segemehl", mode: 'copy'

    script:
    """
    ls
    """
}

/*************************************************************************
* samtools CONVERT SAM TO BAM
*************************************************************************/

process convertSAMtoBAM {
    label 'mapping_samtools'

    input:
    tuple val(name), path(mappings), val(mapper)

    output:
    tuple val(name), path("${mappings.baseName}.bam")

    publishDir "${params.output}/02-mappings/${mapper}", mode: 'copy'

    script:
    """
    samtools view -@ 8 -S -b ${mappings} > ${mappings.baseName}.bam
    """
}

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