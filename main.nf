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
* CHANNELS
**************************/

params.krakenLocalDB = workflow.projectDir + "/assets/kraken2/kraken_db"
params.krakenAssets = workflow.projectDir + "/assets/kraken2/genomes"

/************************** 
* MODULES
**************************/

// preprocessing
include { fastpTrimming} from './modules/preprocessing.nf'
include { fastqcReport } from './modules/generate_reports.nf'

workflow preprocessing {
    take: reads_ch
    main:
        // Trims reads with fastp
        fastpTrimming( reads_ch )
        // Generates fastqc reports
        qc_ch = fastpTrimming.out.concat( reads_ch )
        fastqcReport( qc_ch )
    emit:
        fastpTrimming.out
        fastqcReport.out
}

// mapping with segemehl
include { segemehlIndex; segemehl; segemehlPublish; convertSAMtoBAM } from './modules/map_reads.nf'
include { getStats } from './modules/generate_reports.nf'

workflow segemehl_mapping {
    take:
        preprocessed_reads_ch
        genomes_ch
    main:
        // Indexes the reference genomes
        segemehlIndex( genomes_ch )
        // Maps the reads to the reference genomes
        segemehl_input_ch = segemehlIndex.out.combine( preprocessed_reads_ch.map{ it -> [ it[2], it[0], it[1] ] },  by: 0 ).map{ it -> [ it[3], it[1], it[2], it[4] ] }
        segemehl_input_ch
        segemehl( segemehl_input_ch )
        // Publishes segemehl trns files for inspection
        segemehlPublish( segemehl.out )
        // Converts segemehl's SAM output to BAM file
        convertSAMtoBAM( 
            segemehl.out.map{ it -> [ it[0], it[4], 'segemehl' ] }
            )
        // Runs samtools flagstats on the BAM file
        getStats( segemehl.out.map{ it -> [ it[0], it[4] ] } )
    emit:
        segemehl.out
        convertSAMtoBAM.out
        getStats.out
}

// array filling using numpy
include { fillArrays } from './modules/fill_arrays.nf'

workflow array_filling {
    take:
        segemehl_trns_ch
    main:
        // Fills the arrays
        fillArrays( segemehl_trns_ch )
    emit:
        fillArrays.out
}

// plot heatmaps
include { plotHeatmaps } from './modules/plot_heatmaps.nf'

workflow plot_heatmaps {
    take:
        arrays_ch
    main:
        // Plots the heatmaps
        plotHeatmaps( arrays_ch )
    emit:
        plotHeatmaps.out
}

/************************** 
* WORKFLOW ENTRY POINT
**************************/

include { plotHeatmaps } from './modules/plot_heatmaps.nf'

workflow {
    // parse sample's csv file
    samples_input_ch = Channel
        .fromPath( params.samples, checkIfExists: true )
        .splitCsv()
        .map{ row -> [ "${row[0]}", file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true), "${row[3]}" ] }
    reads_ch = samples_input_ch
        .map{ it -> [ it[0], it[1], it[3] ] }
    genomes_ch = samples_input_ch
        .map{ it -> [ it[3],  it[2] ] }
        .unique()
    // preprocessing workflow
    preprocessing( reads_ch )
    // segemehl workflow
    segemehl_mapping( preprocessing.out[0], genomes_ch )
    // fill arrays with the segemehl output
    array_filling( 
        segemehl_mapping.out[0]
        .map( it -> [ it[0], it[1], it[5], it[6] ] )
    )
    // plot heatmaps using the filled arrays
    plot_heatmaps( array_filling.out )
    
    // Acquire annotation tables
        // If annotations for the arrays are provided, use them

        // If not, we annotate the arrays de novo uding GMMs
    
    // Plot heatmaps with annotations

    // Generate count tables

    // Run differential analysis with DESeq2

    // Generate tables

    // Predict structures

    // Generate circos plots

}

// workflow {
//     // parse sample's csv file
//     Channel.fromPath( params.samples, checkIfExists: true ) \
//         | splitCsv() \
//         | map{ row -> [ "${row[0]}", file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true), "${row[3]}" ] } \
//         | set { samples_input_ch }
//     samples_input_ch.map{ it -> [ it[0], it[1], it[3] ] } \
//         | set { reads_ch }
//     samples_input_ch.map{ it -> [it[3],  it[2] ] } \
//         | unique() \
//         | set { genomes_ch }
//     // preprocessing workflow
//     preprocessing( reads_ch )
//     // segemehl workflow
//     segemehl_mapping( preprocessing.out[0], genomes_ch )
//     // Accumulate trns files mapped to the same genome
//     segemehl_mapping.out[0].groupTuple(by: 1) \
//         | map{ it -> tuple( groupKey(it[1], it[2].size()), it[2] ) } \
//         | set { segemehl_out_grouped_ch }
//     // take a look at the grouped trns files
//     segemehl_out_grouped_ch.view()
// }
