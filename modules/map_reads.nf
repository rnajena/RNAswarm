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
    tuple val(name), path("${reads.baseName}.trns.txt"), path("${reads.baseName}_segemehl.bam")

    publishDir "${params.output}/02-mappings/segemehl", mode: 'copy'

    script:
    """
    segemehl.x -i ${index}\
               -d ${genome}\
               -q ${reads}\
               -S ${reads.baseName}\
               -A ${params.segemehl_accuracy}\
               -U ${params.segemehl_minfragsco}\
               -W ${params.segemehl_minsplicecov}\
               -Z ${params.segemehl_minfraglen}\
               -t ${params.max_cpus}\
               -b > ${reads.baseName}_segemehl.bam
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
  tuple val(name), path(genome) , path("${name}_index")

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
  tuple val(name), path("${reads.baseName}_bwa.sam")

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
  tuple val(name), path(mappings)

  output:
  tuple val(name), path("${mappings.baseName}.bam")

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
  
  publishDir "${params.output}/02-mappings/chim_files", mode: 'copy'
  
  script:
  """
  find_chimeras_rs.py -i ${mapping} -o ${mapping.baseName}.chim
  """
}
