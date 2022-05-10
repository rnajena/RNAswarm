/***********************************************************************
* ART simulator SIMULATE ILLUMINA READS
***********************************************************************/

process simulate_interaction_reads {
  label "simulate_interactions"

  input:
  tuple val(name), path(interactions), path(genome)

  output:
  tuple val(name), path("${name}_interactions.fq")

  publishDir "${params.output}/00-simulated_reads", mode: 'copy'

  script:
  """
  art_templater.py -i ${interactions} -f ${genome} > ${name}_interactions.fasta
  art_illumina -sam -i ${name}_interactions.fasta -l ${params.read_len} -c ${params.rcount_interaction} -o ${name}_interactions
  """
}

process simulate_genome_reads {
  label "simulate_interactions"

  input:
  tuple val(name), path(fasta)

  output:
  tuple val(name), path("${fasta.baseName}.fq")

  publishDir "${params.output}/00-simulated_reads", mode: 'copy'

  script:
  """
  art_illumina -sam -i ${fasta} -l ${params.read_len} -c ${params.rcount_genome} -o ${fasta.baseName}
  """
}

process concatenate_reads {
  label "simulate_interactions"

  input:
  tuple val(name),  path(fq_1), path(fq_2)

  output:
  tuple val(name), path("${name}_concat.fq")

  publishDir "${params.output}/00-simulated_reads", mode: 'copy'

  script:
  """
  cat ${fq_1} ${fq_2} > ${name}_concat.fq
  """
}
