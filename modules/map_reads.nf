/*************************************************************************
* segemehl INDEX
*************************************************************************/
process segemehlIndex {
  label 'mapping_segemehl'

  input:
  tuple val(name), path(genome), val(is_genome_concatenated)

  output:
  tuple val(name), path(genome), val(is_genome_concatenated), path("${name}.idx")

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
    tuple val(name), path(genome), val(is_genome_concatenated), path(index), path(reads)

    output:
    tuple val(name), path("${reads.baseName}.trns.txt"), path("${reads.baseName}*_segemehl.bam")

    publishDir "${params.output}/02-mappings/segemehl", mode: 'copy'

    script:
    if ( is_genome_concatenated )
      """
      segemehl.x -i ${index}\
                -d ${genome}\
                -q ${reads}\
                -S ${reads.baseName}\
                -t ${params.max_cpus}\
                -b > ${reads.baseName}_concat_segemehl.bam
      """
    else
      """
      segemehl.x -i ${index}\
                -d ${genome}\
                -q ${reads}\
                -S ${reads.baseName}\
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
  tuple val(name), path(genome), val(is_genome_concatenated)

  output:
  tuple val(name), path(genome), val(is_genome_concatenated) , path("${name}_index")

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
  tuple val(name), path(genome), val(is_genome_concatenated), path(index), path(reads)

  output:
  tuple val(name), path("${reads.baseName}*_bwa.sam")

  script:
  if ( is_genome_concatenated )
    """
    bwa mem -t ${params.max_cpus} -T 20 ${index}/${genome} ${reads} > ${reads.baseName}_concat_bwa.sam
    """
  else
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
  tuple val(name), path(genome), val(is_genome_concatenated)

  output:
  tuple val(name), path(genome), val(is_genome_concatenated), path("${name}*.ht2")

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
  tuple val(name), path(genome), val(is_genome_concatenated), path(index), path(reads)

  output:
  tuple val(name), path("${reads.baseName}*_hisat2.sam")

  script:
  if ( is_genome_concatenated )
    """
    hisat2 -x ${name} -U ${reads} > ${reads.baseName}_concat_hisat2.sam
    """
  else
    """
    hisat2 -x ${name} -U ${reads} > ${reads.baseName}_hisat2.sam
    """
}
