/*************************************************************************
* Plot heatmaps
*************************************************************************/

process plotHeatmaps {
    label 'RNAswarm'

    input:
    tuple val(sample_name), path(genome), path(arrays)

    output:
    tuple val(sample_name), path("${sample_name}_heatmaps")

    publishDir "${params.output}/04-stats_and_plots/heatmaps", mode: 'copy'

    script:
    """
    mkdir ${sample_name}
    plot_heatmaps.py -d ${arrays} -g ${genome} -o ${sample_name}_heatmaps
    """
}

/*************************************************************************
* Plot heatmaps annotated
*************************************************************************/

process plotHeatmapsAnnotated {
    label 'RNAswarm'

    input:
    tuple val(sample_name), path(genome), path(arrays), path(annotation_table)

    output:
    tuple val(sample_name), path("${sample_name}_heatmaps")

    publishDir "${params.output}/04-stats_and_plots/heatmaps_annotated", mode: 'copy'

    script:
    """
    mkdir ${sample_name}
    plot_heatmaps.py -d ${arrays} -g ${genome} -a ${annotation_table} -o ${sample_name}_heatmaps
    """
}

/*************************************************************************
* make circos table
*************************************************************************/
process makeCircosTable {
    label 'RNAswarm'

    input:
    tuple val(genome_name_01), path(genome_01), val(genome_name_02), path(genome_02), path(results_DESeq2)

    output:
    tuple val(genome_name_01), val(genome_name_02), path("${genome_name_01}_${genome_name_02}_circos_dir")

    script:
    // It would be important to check if the genomes are of the same size
    """
    make_circos_plots.py ${results_DESeq2} -g ${genome_01} -o ${genome_name_01}_${genome_name_02}_circos
    """
}

/*************************************************************************
* run circos
*************************************************************************/
process runCircos {
    label 'RNAswarm'

    input:
    tuple val(genome_name_01), val(genome_name_02), path(circos_dir)

    output:
    tuple val(genome_name_01), val(genome_name_02), path(circos_dir)

    publishDir "", mode: 'copy'

    script:
    """
    cd ${circos_dir}
    circos -conf circos.conf
    """
}
