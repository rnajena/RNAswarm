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

// generate circos plots
include { makeCircosTable; runCircos } from './modules/data_visualization.nf'

workflow generate_circos_plots {
    take:
        differential_analysis_results_ch
    main:
        // Generates circos plots
        makeCircosTable( differential_analysis_results_ch )
        runCircos( makeCircosTable.out )
    emit:
        runCircos.out
}

/************************** 
* WORKFLOW ENTRY POINT
**************************/

// array filling using numpy
include { fillArrays; mergeArrays } from './modules/handle_arrays.nf'
// plot heatmaps
include { plotHeatmaps as plotHeatmapsRaw } from './modules/data_visualization.nf'
include { plotHeatmaps as plotHeatmapsMerged } from './modules/data_visualization.nf'
include { plotHeatmapsAnnotated } from './modules/data_visualization.nf'
include { annotateArrays; mergeAnnotations } from './modules/annotate_interactions.nf'
// differential analysis
include { generateCountTables; mergeCountTables; runDESeq2 } from './modules/differential_analysis.nf'
// generate circos plots
include { makeCircosTable; runCircos } from './modules/data_visualization.nf'
workflow {
    // parse sample's csv file
    samples_input_ch = Channel
        .fromPath( params.samples, checkIfExists: true )
        .splitCsv()
        .map{
            row -> [
                "${row[0]}",                                // sample name
                file("${row[1]}", checkIfExists: true),     // read file
                file("${row[2]}", checkIfExists: true),     // genome file
                "${row[3]}"                                 // group name
            ]
        }
    reads_ch = samples_input_ch
        .map{ it -> [ it[0], it[1], it[3] ] }               // sample name, read file, group name
    genomes_ch = samples_input_ch
        .map{ it -> [ it[3], it[2] ] }                      // group name, genome file
        .unique()

    // preprocessing workflow
    preprocessing( reads_ch )

    // segemehl workflow
    segemehl_mapping( preprocessing.out[0], genomes_ch )

    // fill arrays with the segemehl output
    array_ch = fillArrays(
        segemehl_mapping.out[0]
        .map( it -> [ it[0], it[1], it[5], it[6] ] )        // sample name, trns file, group name, genome
    )

    // plot heatmaps using the filled arrays
    plotHeatmapsRaw( 
        array_ch
        .map( it -> [ it[0], it[3], it[4] ] )               // sample name, genome, array
     )

    // accumulate arrays with the same group name
    groupped_arrays_ch = array_ch
        .groupTuple( by: 2 )                                // This can be streamlined by knowing the number of samples in each group beforehand,
                                                            // but should be fine for now
        .map( it -> [ it[2], it[3][0], it[4].flatten()] )   // group name, genome, arrays
    groupped_arrays_ch
    // merge arrays with the same group name
    merged_arrays_ch = mergeArrays( groupped_arrays_ch )

    // plot heatmaps using the merged arrays
    plotHeatmapsMerged( merged_arrays_ch )

    // Check if annotations are present
    if ( params.annotation_table ) {
        // Create a channel with the annotations
        annotated_arrays_ch = merged_arrays_ch
            .combine( Channel.fromPath( params.annotation_table, checkIfExists: true ) )
        annotated_trns_ch = segemehl_mapping.out[0]
            .map( it -> [ it[0], it[1], it[5] ] ) // sample name, trns file, group name
            .combine( Channel.fromPath( params.annotation_table, checkIfExists: true ) )
    } else {
        // Annotate interactions de novo
        annotated_arrays_ch = annotateArrays( 
            merged_arrays_ch 
            )
        // collect annotations from the annotated_arrays_ch channel and merge them
        mergeAnnotations(
            annotated_arrays_ch
                .collect { it[3] }
        ).view()
        //
        annotated_trns_ch = segemehl_mapping.out[0]
            .map( it -> [ it[0], it[1], it[5] ] ) // sample name, trns file, group name
            .combine( mergeAnnotations.out )
    }

    // Plot the annotations on the heatmaps
    plotHeatmapsAnnotated( annotated_arrays_ch )

    // Generate count tables
    count_tables_ch = generateCountTables( annotated_trns_ch.view() )
    merged_count_tables_ch = mergeCountTables(
        count_tables_ch
            .groupTuple( by: 2 )
            .map( it -> [ it[2], it[1] ] ) // group name, count tables
    )

    // Run differential analysis with DESeq2
    differential_analysis_results_ch = runDESeq2(
        samples_input_ch = Channel
            .fromPath( params.comparisons, checkIfExists: true )
            .splitCsv()
            .combine( merged_count_tables_ch, by: 0 )
            .map( it -> [ it[1], it[0], it[2] ] )
            .combine( merged_count_tables_ch, by: 0 )
            .map( it -> [ it[1], it[2], it[0], it[3] ] )
    )

    // Generate circos files and render plots
    circos_tables_ch = makeCircosTable( differential_analysis_results_ch )
    runCircos( circos_tables_ch )
}
