/************************************************************************
* run extract_regions.py script
*************************************************************************/

process handleTrnsFiles {
    label 'python3'
  
    input:
    tuple val(name), path(genome), path(trns_file), path(sam_file)

    output:
    tuple val(name), path("${sam_file.baseName}_plots"), path("${sam_file.baseName}_heatmaps.log")

    publishDir "${params.output}/03-heatmaps/segemehl", mode: 'copy'

    script:
    """
    extract_regions.py <input_file> <input_file>... -g <genome> -o <output_folder> [-m <min_components> -M <max_components> --step_size <step_size>]
    """
}