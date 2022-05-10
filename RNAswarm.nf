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
        interaction_tables_ch  = Channel.fromPath("${params.input}/*.csv")
                                        .map{ file -> tuple(file.baseName, file)}

        genomes_ch = Channel.fromPath("${params.input}/*.fasta")
                            .map{ file -> tuple(file.baseName, file) }
  
        interaction_reads_ch = interaction_tables_ch.combine(genomes_ch, by: 0)

        simulate_interaction_reads( interaction_reads_ch )

        simulate_genome_reads( genomes_ch )

        concat_ch = simulate_genome_reads.out.join(simulate_interaction_reads.out, by: 0)

        concatenate_reads( concat_ch )
    emit:
        concatenate_reads.out
}

// preprocessing
include { fastpTrimming; fastqcReport } from './modules/preprocess_reads.nf'

workflow preprocessing {
    take: reads_ch
    main:    
        fastpTrimming( reads_ch )

        qc_ch = fastpTrimming.out.concat(reads_ch)

        fastqcReport( qc_ch )
    emit:
        fastpTrimming.out
}

// mapping with segemehl
include { segemehlIndex; segemehl; convertSAMtoBAM } from './modules/map_reads.nf'
include { makeConfusionMatrix_trns } from './modules/handle_trns_files.nf'

workflow segemehl_mapping {
    take: preprocessed_reads_ch
    main:
        genomes_ch = Channel.fromPath("${params.input}/*.fasta")
                            .map{ file -> tuple(file.baseName, file) }
        
        segemehlIndex( genomes_ch )
    
        segemehl_input_ch = segemehlIndex.out.combine(preprocessed_reads_ch, by: 0)
    
        segemehl( segemehl_input_ch )
    emit:
        segemehl.out
}

// mapping with bwa-mem
include { bwaIndex; bwaMem; findChimeras } from './modules/map_reads.nf'
include { makeConfusionMatrix_bwa } from './modules/handle_chim_files.nf'

workflow bwa_mapping {
    take: preprocessed_reads_ch
    main:
        genomes_ch = Channel.fromPath("${params.input}/*.fasta")
                            .map{ file -> tuple(file.baseName, file) }

        bwaIndex( genomes_ch )

        bwa_input_ch = bwaIndex.out.combine(preprocessed_reads_ch, by: 0)

        bwaMem( bwa_input_ch )

        sam_files_ch = bwaMem.out

        convertSAMtoBAM( sam_files_ch )

        findChimeras( convertSAMtoBAM.out )

        confusion_matrix_ch = 

        makeConfusionMatrix_bwa
    emit:
        findChimeras.out
}


// handle chim files
include { handleChimFiles } from './modules/handle_chim_files.nf'

workflow chim_file_handler {
    take: chim_file_ch
    main:
        genomes_ch = Channel.fromPath("${params.input}/*.fasta")
                            .map{ file -> tuple(file.baseName, file) }

        handleChimFiles_input_ch = genomes_ch.combine(chim_file_ch, by: 0)

        handleChimFiles( handleChimFiles_input_ch )
    emit:
        handleChimFiles.out
}

// handle trns files
include { handleTrnsFiles } from './modules/handle_trns_files.nf'

workflow trns_file_handler {
    take: trns_file_ch
    main:
        genomes_ch = Channel.fromPath("${params.input}/*.fasta")
                            .map{ file -> tuple(file.baseName, file) }

        handleTrnsFiles_input_ch = genomes_ch.combine(trns_file_ch, by: 0)

        handleTrnsFiles( handleTrnsFiles_input_ch )
    emit:
        handleTrnsFiles.out
}

/************************** 
* WORKFLOW ENTRY POINT
**************************/

workflow {
    simulate_interactions()
    preprocessing(simulate_interactions.out)
    // bwa workflow
    bwa_mapping(preprocessing.out)
    chim_file_handler(bwa_mapping.out)
    // segemehl workflow
    segemehl_mapping(preprocessing.out)
    trns_file_handler(segemehl_mapping.out)
}
