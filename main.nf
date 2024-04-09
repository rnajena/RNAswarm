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
// make alias for mergeCountTables to merge all count tables
include { mergeCountTables as mergeAllCountTables } from './modules/differential_analysis.nf'
// deduplicate annotations
include { deduplicateAnnotations } from './modules/annotate_interactions.nf'
// generate circos plots
include { makeCircosTable_deseq2; makeCircosTable_count_table; runCircos_single; runCircos_comb } from './modules/data_visualization.nf'
workflow {
    if ( params.help ) {
        println("""
    RNAswarm: differential RNA-RNA interaction probing pipeline based on RNA proximity ligation data

    Usage:
            The typical command for running the pipeline is as follows:
            nextflow run gabriellovate/RNAswarm -profile local,apptainer --samples <SAMPLES_CSV_FILE> --comparisons <COMPARISONS_CSV_FILE> --output <OUTDIR>
            nextflow run gabriellovate/RNAswarm -profile slurm,apptainer --samples <SAMPLES_CSV_FILE> --comparisons <COMPARISONS_CSV_FILE> --output <OUTDIR>
            nextflow run gabriellovate/RNAswarm -profile local,apptainer --samples <SAMPLES_CSV_FILE> --comparisons <COMPARISONS_CSV_FILE> --annotation_table <ANNOTATION_TABLE> --output <OUTDIR>
            nextflow run gabriellovate/RNAswarm -profile slurm,apptainer --samples <SAMPLES_CSV_FILE> --comparisons <COMPARISONS_CSV_FILE> --annotation_table <ANNOTATION_TABLE> --output <OUTDIR>

    Mandatory arguments:
            --samples <SAMPLES_CSV_FILE>          CSV file containing the samples to be processed. The file must have the following format:
                                                <SAMPLE_NAME>,<READ_FILE>,<GENOME_FILE>,<GROUP_NAME>
                                                <SAMPLE_NAME>,<READ_FILE>,<GENOME_FILE>,<GROUP_NAME>
                                                ...
                                                <SAMPLE_NAME>,<READ_FILE>,<GENOME_FILE>,<GROUP_NAME>
                                                where:
                                                - <SAMPLE_NAME> is the name of the sample
                                                - <READ_FILE> is the path to the read file
                                                - <GENOME_FILE> is the path to the genome file
                                                - <GROUP_NAME> is the name of the group to which the sample belongs
            --comparisons <COMPARISONS_CSV_FILE>  CSV file containing the comparisons to be performed. The file must have the following format:
                                                <GROUP_NAME_1>,<GROUP_NAME_2>
                                                <GROUP_NAME_1>,<GROUP_NAME_3>
                                                ...
                                                <GROUP_NAME_2>,<GROUP_NAME_3>
                                                where:
                                                - <GROUP_NAME_X> is the name of the group
            --output <OUTDIR>                     Output directory

    Optional arguments:
            --annotation_table <ANNOTATION_TABLE> CSV, TSV or XLSX file containing the annotations to be used.
            --test                                Run in test mode. This will run the pipeline with a small subset of the data
            --help                                Print this help message
            """)
        exit 0
    }
    // parse sample's csv file
    if ( params.test ) {
        println "Running in test mode"
        samples_input_ch = Channel
            .fromPath( params.samples, checkIfExists: true )
            .splitCsv()
            .map{
                row -> [
                    "${row[0]}",                                // sample name
                    file("${baseDir}${row[1]}", checkIfExists: true),     // read file
                    file("${baseDir}${row[2]}", checkIfExists: true),     // genome file
                    "${row[3]}"                                 // group name
                ]
            }
    } else {
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
        )
        //
        annotated_trns_ch = segemehl_mapping.out[0]
            .map( it -> [ it[0], it[1], it[5] ] ) // sample name, trns file, group name
            .combine( mergeAnnotations.out )
    }

    // Plot the annotations on the heatmaps
    plotHeatmapsAnnotated( 
        annotated_arrays_ch.map( it -> [ it[0], it[1], it[2], it[3] ] ) // sample name, genome, array, annotations
    )

    // Generate count tables
    count_tables_ch = generateCountTables( annotated_trns_ch )
    merged_count_tables_ch = mergeCountTables(
        count_tables_ch
            .groupTuple( by: 2 )
            .map( it -> [ it[2], it[1] ] ) // group name, count tables
    )

    // Merge all count tables independently of the group by collecting all count tables
    merged_count_tables_all_ch = mergeAllCountTables(
        count_tables_ch
            .map( it -> [ it[1] ] ) // count tables
            .collect()
            .map( it -> [ "all", it ] ) // group name, count tables
    )

    // Deduplicate annotations
    deduplicate_annotations_input_ch = merged_count_tables_all_ch // group_name, merged_count_table
            .combine( mergeAnnotations.out ) // merged_annotations
            .map( it -> [ it[0], it[2], it[3] ] ) // group name, count table, annotations
    deduplicate_annotations_input_ch.view()
    deduplicateAnnotations( deduplicate_annotations_input_ch )

    // Run differential analysis with DESeq2
    samples_input_ch = Channel
            .fromPath( params.comparisons, checkIfExists: true )
            .splitCsv()
            .combine( merged_count_tables_ch, by: 0 )
            .map( it -> [ it[1], it[0], it[2] ] )
            .combine( merged_count_tables_ch, by: 0 )
            .map( it -> [ it[1], it[2], it[0], it[3] ] )

    runDESeq2( samples_input_ch )

    // Generate circos files
    if ( params.annotation_table ) {
        circos_deseq2_ch = runDESeq2.out
                        .combine( genomes_ch, by: 0 )
                        .map( it -> [ it[1], it[0], it[3], it[2] ] )
                        .combine( genomes_ch, by: 0 )
                        .map( it -> [ it[1], it[2], it[0], it[4], it[3] ] )
                        .combine( Channel.fromPath( params.annotation_table, checkIfExists: true ) )
        circos_count_table_ch = merged_count_tables_ch
                        .combine( genomes_ch, by: 0 )
                        .map( it -> [ it[0], it[2], it[1] ] )
                        .combine( Channel.fromPath( params.annotation_table, checkIfExists: true ))
    } else {
        circos_deseq2_ch = runDESeq2.out
                        .combine( genomes_ch, by: 0 )
                        .map( it -> [ it[1], it[0], it[3], it[2] ] )
                        .combine( genomes_ch, by: 0 )
                        .map( it -> [ it[1], it[2], it[0], it[4], it[3] ] )
                        .combine( deduplicateAnnotations.out )
        circos_count_table_ch = merged_count_tables_ch
                        .combine( genomes_ch, by: 0 )
                        .map( it -> [ it[0], it[2], it[1] ] )
                        .combine( deduplicateAnnotations.out )
    }

    // Create circos tables
    makeCircosTable_deseq2( circos_deseq2_ch )
    makeCircosTable_count_table( circos_count_table_ch )


    // Render circos plots
    runCircos_single( makeCircosTable_count_table.out )
    runCircos_comb( makeCircosTable_deseq2.out )
}