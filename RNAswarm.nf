#!/usr/bin/env nextflow

/* 
* RNAswarm: differential RNA-RNA interaction probing pipeline based on RNA proximity ligation data
*
* Authors: 
* - Gabriel Lencioni Lovate <gabriel.lencioni.lovate@uni-jena.de>
* - Kevin Lamkiewicz <kevin.lamkiewicz@uni-jena.de>
*/

nextflow.enable.dsl=2

/************************** 
* MODULES
**************************/

// interaction simulation
include { simulate_interaction_reads; simulate_genome_reads; concatenate_reads } from './modules/interaction_simulator.nf'

workflow simulate_interactions {
    main:
        interaction_tables_ch  = Channel
                                 .fromPath("${params.input}/*.csv")
                                 .map{ file -> tuple(file.baseName, file)}

        genomes_ch = Channel
                    .fromPath("${params.input}/*.fasta")
                    .map{ file -> tuple(file.baseName, file) }
  
        interaction_reads_ch = interaction_tables_ch.combine(genomes_ch, by: 0)

        simulate_interaction_reads( interaction_reads_ch )

        simulate_genome_reads( genomes_ch )

        concat_ch = simulate_genome_reads.out.join(simulate_interaction_reads.out, by: 0).view()

        concatenate_reads( concat_ch )
}

/************************** 
* WORKFLOW ENTRY POINT
**************************/

workflow {
    simulate_interactions()
}
