/*************************************************************************
* Plot heatmaps
*************************************************************************/

process plotHeatmaps {
    label 'RNAswarm'

    input:
    tuple val(sample_name), path(genome), path(arrays)

    output:
    tuple val(sample_name), path("${sample_name}_heatmaps")

    publishDir "${params.output}/05-stats_and_plots/heatmaps", mode: 'copy'

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
    label 'RNAswarm'

    input:
    tuple val(sample_name), path(genome), path(arrays), path(annotation_table)

    output:
    tuple val(sample_name), path("${sample_name}_heatmaps")

    publishDir "${params.output}/05-stats_and_plots/heatmaps_annotated", mode: 'copy'

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
    label 'RNAswarm'

    input:
    tuple val(genome_name), path(genome), path(count_table), path(annotation_table)

    output:
    tuple val(genome_name), path("${genome_name}_circos")

    script:
    """
    make_circos_files.py -c ${count_table} -a ${annotation_table} -g ${genome} -o ${genome_name}_circos
    """
}

/*************************************************************************
* make circos table for count_table results
*************************************************************************/
process makeCircosTable_count_table_30 {
    label 'RNAswarm'

    input:
    tuple val(genome_name), path(genome), path(count_table), path(annotation_table)

    output:
    tuple val(genome_name), path("${genome_name}_circos_30")

    script:
    """
    make_circos_files.py -c ${count_table} -a ${annotation_table} -g ${genome} -o ${genome_name}_circos_30 --number_of_top_hits 30
    """
}

/*************************************************************************
* make circos table for count_table results
*************************************************************************/
process makeCircosTable_count_table_40 {
    label 'RNAswarm'

    input:
    tuple val(genome_name), path(genome), path(count_table), path(annotation_table)

    output:
    tuple val(genome_name), path("${genome_name}_circos_40")

    script:
    """
    make_circos_files.py -c ${count_table} -a ${annotation_table} -g ${genome} -o ${genome_name}_circos_40 --number_of_top_hits 40
    """
}

/*************************************************************************
* make circos table for count_table results
*************************************************************************/
process makeCircosTable_count_table_50 {
    label 'RNAswarm'

    input:
    tuple val(genome_name), path(genome), path(count_table), path(annotation_table)

    output:
    tuple val(genome_name), path("${genome_name}_circos_50")

    script:
    """
    make_circos_files.py -c ${count_table} -a ${annotation_table} -g ${genome} -o ${genome_name}_circos_50 --number_of_top_hits 50
    """
}

/*************************************************************************
* make circos table for deseq2 results
*************************************************************************/
process makeCircosTable_deseq2 {
    label 'RNAswarm'

    input:
    tuple val(genome_name_01), path(genome_01), val(genome_name_02), path(genome_02), path(results_DESeq2), path(annotation_table)

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
    label 'RNAswarm'

    input:
    tuple val(genome_name), path(circos_dir)

    output:
    tuple val(genome_name), path(circos_dir)

    publishDir "${params.output}/08-circos_plots", mode: 'copy'

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
    label 'RNAswarm'

    input:
    tuple val(genome_name_01), val(genome_name_02), path(circos_dir)

    output:
    tuple val(genome_name_01), val(genome_name_02), path(circos_dir)

    publishDir "${params.output}/08-circos_plots", mode: 'copy'

    script:
    """
    cd ${circos_dir}
    circos -conf circos.conf
    """
}
