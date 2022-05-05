/***********************************************************************
* ART simulator SIMULATE ILLUMINA READS
***********************************************************************/

process make_templates {
  label "python3"

  input:
  tuple val(name), path(interactions), path(genome)

  output:
  tuple val(name), path("${name}_interactions.fasta")

  script:
  """
  art_templater.py -i ${interactions} -f ${genome} > ${name}.fasta
  """

}

process simulate_interactions {
  label "simulate_interactions"

  input:
  tuple val(name), path(interactions), path(genome)

  output:
  tuple val(name), path("${name}_concat.fastq")

  publishDir "${params.output}/00-simulated_reads", mode: 'copy'

  script:
  """
  art_illumina -i ${name}_interactions.fasta -l ${params.read_len} -c ${params.rcount_interaction} -o ${name}_interactions.fastq
  art_illumina -i ${genome} -l ${params.read_len} -c ${params.rcount_genome} -o ${name}_genome.fastq
  cat ${name}_genome.fastq ${name}_interactions.fastq > ${name}_concat.fastq
  """
}