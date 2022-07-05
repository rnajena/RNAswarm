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
        segemehl( segemehl_input_ch )
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

// mapping with bwa-mem
include { bwaIndex; bwaMem; findChimeras; convertSAMtoBAM } from './modules/map_reads.nf'

workflow bwa_mapping {
    take:
        preprocessed_reads_ch
        preprocessed_genomes_ch
    main:
        bwaIndex( preprocessed_genomes_ch )

        bwa_input_ch = bwaIndex.out.combine(preprocessed_reads_ch, by: 0)

        bwaMem( bwa_input_ch )

        convertSAMtoBAM( 
            bwaMem.out.map{ it -> [ it[0], it[1], 'bwa-mem' ] }
            )

        findChimeras( convertSAMtoBAM.out )

        getStats( convertSAMtoBAM.out )
    emit:
        convertSAMtoBAM.out
        findChimeras.out
        getStats.out
}

// mapping with hisat2
include { hiSat2Index; hiSat2 } from './modules/map_reads.nf'
include { concatenateFasta } from './modules/preprocessing.nf'

workflow hisat2_mapping {
    take:
        preprocessed_reads_ch
        preprocessed_genomes_ch
    main:
        hiSat2Index( preprocessed_genomes_ch )

        hisat2_input_ch = hiSat2Index.out.combine(preprocessed_reads_ch, by: 0)

        hiSat2( hisat2_input_ch )

        convertSAMtoBAM( 
            hiSat2.out.map{ it -> [ it[0], it[1], 'hisat2' ] }
            )
        
        getStats( convertSAMtoBAM.out )
    emit:
        convertSAMtoBAM.out
        getStats.out
}

// handle chim files
include { handleChimFiles } from './modules/handle_chim_files.nf'

workflow chim_file_handler {
    take: chim_file_ch
    main:
        genomes_ch = Channel.fromPath("${params.input}/genomes/*.fasta")
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
        genomes_ch = Channel.fromPath("${params.input}/genomes/*.fasta")
                            .map{ file -> tuple(file.baseName, file) }

        handleTrnsFiles_input_ch = genomes_ch.combine(trns_file_ch, by: 0)

        handleTrnsFiles( handleTrnsFiles_input_ch )
    emit:
        handleTrnsFiles.out
}

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
                    .map { row -> [ "${row[0]}", file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true) ] }
                    // csv table with columns: sample_name, reads, genome
    reads_ch = samples_input_ch.map{ it -> [ it[0], it[1] ] }
    genomes_ch = samples_input_ch.map{ it -> [ it[0], it[2] ] }
    // preprocessing workflow
    preprocessing( reads_ch )
    // segemehl workflow
    segemehl_mapping( preprocessing.out[0], genomes_ch )
    trns_file_handler( segemehl_mapping.out[0] )
    // bwa workflow
    bwa_mapping( preprocessing.out[0], genomes_ch )
    chim_file_handler( bwa_mapping.out[1] )
    // hisat2 workflow
    hisat2_mapping( preprocessing.out[0], genomes_ch )
    // run Kraken2
    makeKrakenDatabase()
    kraken_ch = reads_ch.combine(makeKrakenDatabase.out)
    runKraken( kraken_ch )
    // run sotrmerna
    makeSortmernaDatabase()
    sortmerna_ch = reads_ch.combine(makeSortmernaDatabase.out)
    runSortmerna( sortmerna_ch )
    // generate reports
    logs_ch = bwa_mapping.out[2]
                         .mix( segemehl_mapping.out[2], hisat2_mapping.out[1], preprocessing.out[1], runKraken.out )
                         .collect()
    runMultiQC( logs_ch )
}
