/*************************************************************************
# Generate count tables
*************************************************************************/
process generateCountTables {
    label 'python3'

    input:
    // for the input I only need one annotation and one trns file
    tuple val(sample_name), path(annotation_table), path(trns_file)

    output:
    tuple val(sample_name), path("${sample_name}_count_table.csv")

    publishDir "${params.output}/06-count_analysis/count_tables", mode: 'copy'

    script:
    """
    make_counttable.py ${trns_file} -a ${annotation_table} -o ${sample_name}_count_table.csv
    """
}


/*************************************************************************
* run DESeq2
*************************************************************************/
process runDESeq2 {
    label 'r'

    input:
    tuple val(genome_name_01), path(genome_01), path(count_table_01), val(genome_name_02), path(genome_02) ,path(count_table_02)

    output:
    tuple val(genome_name_01), path(genome_01), val(genome_name_02), path(genome_02), path("${genome_name_01}_vs_${genome_name_02}_DESeq2.csv")

    publishDir "${params.output}/04-stats_and_plots", mode: 'copy'

    script:
    """
    run_DESeq2.r --count_table1 ${count_table_01} --alias1 ${genome_name_01}\
                 --count_table2 ${count_table_02} --alias2 ${genome_name_02}\
                 --output_file ${genome_name_01}_vs_${genome_name_02}_DESeq2.csv
    """
}



