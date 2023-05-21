/*************************************************************************
# Generate count tables
*************************************************************************/

process generateCountTables {
    label 'python3'

    input:
    // for the input I only need one annotation and one trns file
    tuple val(sample_name), path(trns_files), val(genome_name), path(genome), path(arrays)

    output:
    tuple val(sample_name), path("${sample_name}_heatmaps")

    publishDir "${params.output}/05-stats_and_plots/heatmaps", mode: 'copy'

    script:
    """
    make_counttable.py 
    """
}
