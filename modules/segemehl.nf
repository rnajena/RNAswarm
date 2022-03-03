nextflow.enable.dsl=2

params.reads = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/reads'
params.genomes = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/genomes'
params.mappings = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/segemehl_mappings'

params.segemehl_accuracy = 9
params.segemehl_minfragsco = 15
params.segemehl_minsplicecov = 80
params.segemehl_minfraglen = 15
params.segemehl_exclclipping = 0
params.segemehl_threads = 48

/***********************************************************************
* segemehl INDEX
*************************************************************************/
process segemehlIndex {
    label 'segemehl'
    
    cpus 8
    time '12h'
    executor 'slurm'
    conda '../envs/segemehl.yaml'

    input:
    tuple val(name), path(genome)

    output:
    tuple val(name), path(genome), path("${name}.idx")

    script:
    """
    segemehl.x -d ${genome} -x ${name}.idx
    """
}

/***********************************************************************
* segemehl RUN
*************************************************************************/
process segemehl {
    label 'segemehl'
    
    cpus "${params.segemehl_threads}"
    time '12h'
    executor 'slurm'
    conda '../envs/segemehl.yaml'

    input:
    tuple val(name), path(genome), path(index), path(reads)

    output:
    tuple val(name), path("${reads.baseName}.sam"), path("${reads.baseName}.trns.txt") 

    publishDir "${params.mappings}", mode: 'copy'

    script:
    """
    segemehl.x -i ${index}\
               -d ${genome}\
               -q ${reads}\
               -S ${reads.baseName}\
               -A ${params.segemehl_accuracy}\
               -U ${params.segemehl_minfragsco}\
               -W ${params.segemehl_minsplicecov}\
               -Z ${params.segemehl_minfraglen}\
               -t ${params.segemehl_threads}\
               > ${reads.baseName}.sam
    """
}

/************************************************************************
* runs complete segemehl workflow
************************************************************************/
workflow {

    main:

        genomes_ch = Channel
                    .fromPath("$params.genomes/*.fasta")
                    .map{ file -> tuple(file.baseName, file) }
        
        reads_ch = Channel
                .fromPath("$params.reads/*/*.fastq")
                .map{ file -> tuple(file.baseName[0..-14], file) }
        segemehlIndex(genomes_ch)

        segemehl_input_ch = segemehlIndex.out.combine(reads_ch, by: 0)

        segemehl( segemehl_input_ch )
}