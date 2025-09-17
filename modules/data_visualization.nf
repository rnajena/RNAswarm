/*************************************************************************
* Plot heatmaps
*************************************************************************/

process plotHeatmaps {
    label 'RNAswarm_small'
    publishDir "${params.output}/05-stats_and_plots/heatmaps", mode: 'copy'

    input:
    tuple val(sample_name), path(genome), path(arrays)

    output:
    tuple val(sample_name), path("${sample_name}_heatmaps")

    script:
    """
    mkdir tmp
    mkdir ${sample_name}_heatmaps
    plot_heatmaps.py -d ${arrays} -g ${genome} -o ${sample_name}_heatmaps
    """
}

/*************************************************************************
* Plot heatmaps annotated
*************************************************************************/

process plotHeatmapsAnnotated {
    label 'RNAswarm_small'
    publishDir "${params.output}/05-stats_and_plots/heatmaps_annotated", mode: 'copy'

    input:
    tuple val(sample_name), path(genome), path(arrays), path(annotation_table)

    output:
    tuple val(sample_name), path("${sample_name}_heatmaps")

    script:
    """
    mkdir ${sample_name}_heatmaps
    plot_heatmaps.py -d ${arrays} -g ${genome} -a ${annotation_table} -o ${sample_name}_heatmaps
    """
}

/*************************************************************************
* make circos table for count_table results
*************************************************************************/
process makeCircosTable_count_table {
    label 'RNAswarm_small'

    input:
    tuple val(genome_name), path(genome), path(genome_count_table), val(group_name), path(annotation_table), path(global_count_table)

    output:
    tuple val(genome_name), path("${genome_name}_circos")

    script:
    """
    make_circos_files.py -c ${genome_count_table} -a ${annotation_table} -g ${genome} -o ${genome_name}_circos
    """
}

/*************************************************************************
* make circos table for deseq2 results
*************************************************************************/
process makeCircosTable_deseq2 {
    label 'RNAswarm_small'

    input:
    tuple val(genome_name_01), path(genome_01), val(genome_name_02), path(genome_02), path(results_DESeq2), val(group_name), path(annotation_table), path(global_count_table)

    output:
    tuple val(genome_name_01), val(genome_name_02), path("${genome_name_01}_${genome_name_02}_circos")

    script:
    // It would be important to check if the genomes are of the same size
    """
    make_circos_files.py -d ${results_DESeq2} -a ${annotation_table} -g ${genome_01} -o ${genome_name_01}_${genome_name_02}_circos
    """
}

/*************************************************************************
* run circos for single analysis
*************************************************************************/
process runCircos_single {
    label 'RNAswarm_small'
    publishDir "${params.output}/08-circos_plots", mode: 'copy'

    input:
    tuple val(genome_name), path(circos_dir)

    output:
    tuple val(genome_name), path(circos_dir)

    script:
    """
    cd ${circos_dir}
    circos -conf circos.conf
    """
}

/*************************************************************************
* run circos for comparative analysis
*************************************************************************/
process runCircos_comb {
    label 'RNAswarm_small'
    publishDir "${params.output}/08-circos_plots", mode: 'copy'

    input:
    tuple val(genome_name_01), val(genome_name_02), path(circos_dir)

    output:
    tuple val(genome_name_01), val(genome_name_02), path(circos_dir)

    script:
    """
    cd ${circos_dir}
    circos -conf circos.conf
    """
}
