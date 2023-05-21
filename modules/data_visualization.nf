/*************************************************************************
* Plot heatmaps
*************************************************************************/

process plotHeatmaps {
    label 'python3'

    input:
    tuple val(sample_name), path(trns_files), val(genome_name), path(genome), path(arrays)
    // old input: tuple val(genome_name), path(genome_file), path(trns_files)

    output:
    tuple val(sample_name), path("${sample_name}_heatmaps")

    publishDir "${params.output}/04-stats_and_plots/heatmaps", mode: 'copy'

    script:
    """
    mkdir ${sample_name}
    plot_heatmaps.py -a ${arrays} -g ${genome_file} -o ${sample_name}_heatmaps
    """
}


/*************************************************************************
* make circos table
*************************************************************************/
process makeCircosTable {
    label 'python3'

    input:


    output:


    script:
    """
    
    """
}


/*************************************************************************
* run circos
*************************************************************************/
process runCircos {
    label 'circos'

    input:

    output:

    publishDir "", mode: 'copy'

    script:
    """
    
    """
}
