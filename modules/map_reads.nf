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
    tuple val(name), path(genome), path(index), path(reads)

    output:
    tuple val(name), path("${reads.baseName}.trns.txt"), path("${reads.baseName}*_segemehl.sam")

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
* bwa-mem INDEX
*************************************************************************/

process bwaIndex {
    label 'mapping_bwa'

    input:
    tuple val(name), path(genome)

    output:
    tuple val(name), path(genome), path("${name}_index")

    script:
    """
    mkdir ${name}_index
    bwa index ${genome} -p ${name}_index/${genome}
    cp ${genome} ${name}_index
    """
    // The cp is a bit hacky, maybe there is a more elegant way of doing this
}

/*************************************************************************
* bwa-mem RUN
*************************************************************************/

process bwaMem {
    label 'mapping_bwa'

    input:
    tuple val(name), path(genome), path(index), path(reads)

    output:
    tuple val(name), path("${reads.baseName}*_bwa.sam")

    script:
    """
    bwa mem -t ${params.max_cpus} -T 20 ${index}/${genome} ${reads} > ${reads.baseName}_bwa.sam
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
* generate .chim files
*************************************************************************/

process findChimeras {
    label 'python2'

    input:
    tuple val(name), path(mapping)
  
    output:
    tuple val(name), path("${mapping.baseName}.chim")
  
    publishDir "${params.output}/02-mappings/bwa-mem", mode: 'copy'
  
    script:
    """
    find_chimeras.py -i ${mapping} -o ${mapping.baseName}.chim
    """
}

/*************************************************************************
* HiSat2 INDEX
*************************************************************************/

process hiSat2Index {
    label 'mapping_hisat2'

    input:
    tuple val(name), path(genome)

    output:
    tuple val(name), path(genome), path("${name}*.ht2")

    script:
    """
    hisat2-build -p ${params.max_cpus} ${genome} ${name}
    """
}

/*************************************************************************
* HiSat2 RUN
*************************************************************************/

process hiSat2 {
    label 'mapping_hisat2'

    input:
    tuple val(name), path(genome), path(index), path(reads)

    output:
    tuple val(name), path("${reads.baseName}*_hisat2.sam")

    script:
    """
    hisat2 -x ${name} -U ${reads} > ${reads.baseName}_hisat2.sam
    """
}
