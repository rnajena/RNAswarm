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
include { fastpTrimming} from './modules/data_preprocessing.nf'
include { fastqcReport } from './modules/reports_generation.nf'

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
include { getStats } from './modules/reports_generation.nf'

workflow segemehl_mapping {
    take:
        preprocessed_reads_ch
        genomes_ch
    main:
        // Indexes the reference genomes
        segemehlIndex( genomes_ch )
        // Maps the reads to the reference genomes
        segemehl_input_ch = segemehlIndex.out.combine(            // group name, genome file, idx file
            preprocessed_reads_ch                                 // sample name, reads file, group name 
            .map{ it -> [ it[2], it[0], it[1] ] },  by: 0         // group name, sample name, reads file
            )
            .map{ it -> [ it[3], it[4], it[0], it[1], it[2] ] }   // sample name, reads file, group name, genome file, idx file
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

// handle annotations
include { annotateArrays; } from './modules/annotate_interactions.nf'

workflow annotate_interactions {
    take:
        arrays_ch
    main:
        // Annotates the arrays
        annotateArrays( arrays_ch )
    emit:
        annotateArrays.out
}

// run differential analysis
include { generateCountTables; runDESeq2 } from './modules/differential_analysis.nf'

workflow differential_analysis {
    take:
        arrays_ch
    main:
        // Generates the count tables
        generateCountTables( arrays_ch )
        // Performs the differential analysis
        runDESeq2( generateCountTables.out )
    emit:
        runDESeq2.out
}

// generate circos plots
include { generateCircosPlots } from './modules/generate_circos_plots.nf'

workflow generate_circos_plots {
    take:
        structures_ch
        count_tables_ch
        differential_analysis_results_ch
    main:
        // Generates circos plots
        generateCircosPlots(  )
    emit:
        generateCircosPlots.out
}

/************************** 
* WORKFLOW ENTRY POINT
**************************/

// array filling using numpy
include { fillArrays; mergeArrays } from './modules/handle_arrays.nf'
// plot heatmaps
include { plotHeatmaps } from './modules/data_visualization.nf'

workflow {
    // parse sample's csv file
    samples_input_ch = Channel
        .fromPath( params.samples, checkIfExists: true )
        .splitCsv()
        .map{
            row -> [
                "${row[0]}",                             // sample name
                file("${row[1]}", checkIfExists: true),  // read file
                file("${row[2]}", checkIfExists: true),  // genome file
                "${row[3]}"                              // group name
            ]
        }
    reads_ch = samples_input_ch
        .map{ it -> [ it[0], it[1], it[3] ] }            // sample name, read file, group name
    genomes_ch = samples_input_ch
        .map{ it -> [ it[3],  it[2] ] }                  // group name, genome file
        .unique()

    // preprocessing workflow
    preprocessing( reads_ch )

    // segemehl workflow
    segemehl_mapping( preprocessing.out[0], genomes_ch )

    // fill arrays with the segemehl output
    array_ch = fillArrays(
        segemehl_mapping.out[0]
        .map( it -> [ it[0], it[1], it[5], it[6] ] ) // sample name, trns file, group name, genome
    )

    // plot heatmaps using the filled arrays
    plot_heatmaps( array_ch )

    // accumulate arrays with the same group name
    groupped_arrays_ch = array_ch.groupTuple( by: 2 ) // This can be streamlined by knowing the number of samples in each group beforehand, but should be fine for now
        .map( it -> [ it[2], it[3].unique(), it[4].flatten()] ) // group name, genome, arrays

    // merge arrays with the same group name
    merged_arrays_ch = mergeArrays( groupped_arrays_ch )

    // plot heatmaps using the merged arrays
    plot_heatmaps( merged_arrays_ch ) // I've already called this function, but I'm not sure if it's a problem

    // Check if annotations are present
    if ( params.annotations ) {
        // Create a channel with the annotations
        annotated_arrays_ch = merged_arrays_ch
            .map( it -> [ it[0], it[2]])
            .combine( Channel.fromPath( params.annotations, checkIfExists: true ) )
        annotated_trns_ch = segemehl_mapping.out[0]
            .map( it -> [ it[0], it[1], it[5], it[6] ] ) // sample name, trns file, group name, genome
            .combine( Channel.fromPath( params.annotations, checkIfExists: true ) )
    } else {
        // Annotate interactions de novo
        annotated_arrays_ch = annotate_interactions( array_filling.out ).out
    }

    // Plot the annotations on the heatmaps
    plot_heatmaps( annotated_arrays_ch )

    // Generate count tables
    count_tables_ch = generate_count_tables( annotated_trns_ch )

    // Run differential analysis with DESeq2


    // Generate circos plots
    

    }
