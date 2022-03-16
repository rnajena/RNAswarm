nextflow.enable.dsl=2

// filepaths
params.reads = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/trimmed_reads'
params.genomes = '/home/ru27wav/Projects/gl_iav-splash_freiburg/data/dadonaite_2019/genomes'
params.mappings = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/mappings'

// segemehl parameters
params.segemehl_accuracy = 9
params.segemehl_minfragsco = 15
params.segemehl_minsplicecov = 80
params.segemehl_minfraglen = 15
params.segemehl_exclclipping = 0
params.segemehl_threads = 48

// slurm parameters
params.slurm_queue = 'b_standard,b_fat,s_standard,s_fat'

/*************************************************************************
* segemehl INDEX
*************************************************************************/
process segemehlIndex {
  label 'mapping'
  
  cpus 8
  time '12h'
  executor 'slurm'
  conda '../envs/mapping_segemehl.yaml'

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
  label 'mapping'
  
  cpus "${params.segemehl_threads}"
  time '12h'
  executor 'slurm'
  queue 'b_standard,b_fat,s_standard,s_fat'
  conda '../envs/mapping_segemehl.yaml'

  input:
  tuple val(name), path(genome), path(index), path(reads)

  output:
  tuple val(name), path("${reads.baseName}_segemehl.sam"), path("${reads.baseName}.trns.txt") 

  publishDir "${params.mappings}", mode: 'copy'

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
             -t ${params.segemehl_threads}\
             > ${reads.baseName}_segemehl.sam
    """
}

/*************************************************************************
* bwa-mem INDEX
*************************************************************************/

process bwaIndex {
  label 'mapping'
  
  cpus 8
  time '12h'
  executor 'slurm'
  conda '../envs/mapping_bwa.yaml'

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
  label 'mapping'

  cpus "${params.segemehl_threads}"
  time '12h'
  executor 'slurm'
  queue 'b_standard,b_fat,s_standard,s_fat'
  conda '../envs/mapping_bwa.yaml'

  input:
  tuple val(name), path(genome), path(index), path(reads)

  output:
  tuple val(name), path("${reads.baseName}_bwa.sam")

  publishDir "${params.mappings}", mode: 'copy'

  script:
  """
  bwa mem -t ${params.segemehl_threads} -T 20 ${index}/${genome} ${reads} > ${reads.baseName}_bwa.sam
  """
}

/*************************************************************************
* samtools CONVERT SAM TO BAM
*************************************************************************/

process convertSAMtoBAM {
  label 'mapping'

  cpus 8
  time '12h'
  executor 'slurm'
  conda '../envs/mapping_samtools.yaml'

  input:
  tuple val(name), path(mappings)

  output:
  tuple val(name), path("${mappings.baseName}.bam")

  publishDir "${params.mappings}", mode: 'copy'

  script:
  """
  samtools view -@ 8 -S -b ${mappings} > ${mappings.baseName}.bam
  """
}

/************************************************************************
* runs complete mapping workflow
************************************************************************/
workflow {

  main:

    genomes_ch = Channel
                .fromPath("${params.genomes}/*.fasta")
                .map{ file -> tuple(file.baseName, file) }.view()
        
    reads_ch = Channel
              .fromPath("${params.reads}/*.fastq")
              .map{ file -> tuple(file.baseName[0..-20], file) }.view()
        
    segemehlIndex( genomes_ch )
    
    segemehl_input_ch = segemehlIndex.out.combine(reads_ch, by: 0)
    
    segemehl( segemehl_input_ch )
    
    // bwaIndex( genomes_ch )
    
    // bwa_input_ch = bwaIndex.out.combine(reads_ch, by: 0)

    // bwaMem( bwa_input_ch )

    // to_convert_ch = bwaMem.out.concat(segemehl.out)

    // convertSAMtoBAM( to_convert_ch )
}