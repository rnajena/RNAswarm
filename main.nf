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
include { segemehlIndex; segemehl; segemehlPublish } from './modules/map_reads.nf'
include { getStats } from './modules/generate_reports.nf'

workflow segemehl_mapping {
    take:
        preprocessed_reads_ch
        genomes_ch
    main:
        // Indexes the reference genomes
        segemehlIndex( genomes_ch )
        // Maps the reads to the reference genomes
        segemehl_input_ch = segemehlIndex.out.combine(preprocessed_reads_ch, by: 0)
        segemehl( segemehl_input_ch ).view()
        // Publishes segemehl trns files for inspection
        segemehlPublish( segemehl.out )
        // Converts segemehl's SAM output to BAM file
        convertSAMtoBAM( 
            segemehl.out.map{ it -> [ it[0], it[2], 'segemehl' ] }
            )
        // Runs samtools flagstats on the BAM file
        getStats( segemehl.out.map{ it -> [ it[0], it[2] ] } )
    emit:
        segemehl.out
        convertSAMtoBAM.out
        getStats.out
}

// handle trns files
include { handleTrnsFiles } from './modules/handle_trns_files.nf'

workflow trns_file_handler {
    take:
        trns_file_ch
        genomes_ch
    main:
        handleTrnsFiles_input_ch = genomes_ch.combine( trns_file_ch, by: 0 )

        handleTrnsFiles( handleTrnsFiles_input_ch )
    emit:
        handleTrnsFiles.out
}

// extract regions from trns files
// include { extractRegions } from './modules/extract_regions.nf'

// make plots from trns files
// include { make_plots } from './modules/make_plots.nf'

// workflow plot_trns_files {
//     take:
//         segemehl
//     main:
        
// }


/************************** 
* WORKFLOW ENTRY POINT
**************************/

include { runMultiQC; makeKrakenDatabase; runKraken; makeSortmernaDatabase; runSortmerna } from './modules/generate_reports.nf'
include { makeDeseq2Table } from './modules/differential_analysis.nf'

workflow {
    // parse sample's csv file
    samples_input_ch = Channel
                    .fromPath( params.samples, checkIfExists: true )
                    .splitCsv()
                    .map { row -> [ "${row[0]}", file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true), "${row[3]}" ] }
    reads_ch = samples_input_ch.map{ it -> [ it[0], it[1] ] }
    genomes_ch = samples_input_ch.map{ it -> [ it[0], it[2] ] }
    // preprocessing workflow
    preprocessing( reads_ch )
    // segemehl workflow
    segemehl_mapping( preprocessing.out[0], genomes_ch )
    trns_file_handler( segemehl_mapping.out[0], genomes_ch )
}
