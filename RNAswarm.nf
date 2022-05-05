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
include { simulate_interactions } from './modules/simulate_interactions.nf'

workflow sim_interactions {
    main:
        interaction_tables_ch  = Channel
                                 .fromPath("${params.input}/*.csv")
                                .map{ file -> tuple(file.baseName, file)}.view()

        genomes_ch = Channel
                    .fromPath("${params.input}/*.fasta")
                     .map{ file -> tuple(file.baseName, file) }.view()
  
        input_ch = interaction_tables_ch.combine(genomes_ch, by: 0).view()

        simulate_interactions( input_ch )
    
    //emit:
        //simulate_interactions.out
}

/************************** 
* WORKFLOW ENTRY POINT
**************************/

workflow {
    // read_simulation
    sim_interactions()
}
