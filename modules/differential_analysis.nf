/*************************************************************************
# Generate count tables
*************************************************************************/
process generateCountTables {
    label 'RNAswarm'

    input:
    // for the input I only need one annotation and one trns file
    tuple val(sample_name), path(trns_file), val(group_name), path(annotation_table)

    output:
    tuple val(sample_name), path("${sample_name}_count_table.tsv"), val(group_name)

    publishDir "${params.output}/07-count_analysis/count_tables", mode: 'copy'

    script:
    """
    make_counttable.py ${trns_file} -a ${annotation_table} -o ${sample_name}_count_table.tsv
    """
}

/*************************************************************************
# Merge count tables
*************************************************************************/
process mergeCountTables {
    label 'RNAswarm'

    input:
    tuple val(group_name), path(count_tables)
    
    output:
    tuple val(group_name), path("${group_name}_count_table.tsv")

    publishDir "${params.output}/07-count_analysis/count_tables", mode: 'copy'

    script:
    """
    merge_counttable.py ${count_tables} -o ${group_name}_count_table.tsv
    """
}

/*************************************************************************
* run DESeq2
*************************************************************************/
process runDESeq2 {
    label 'RNAswarm'

    input:
    tuple val(group_name_01), path(count_table_01), val(group_name_02), path(count_table_02)

    output:
    tuple val(group_name_01), val(group_name_02), path("${group_name_01}_vs_${group_name_02}_DESeq2.tsv")

    publishDir "${params.output}/07-count_analysis/deseq2", mode: 'copy'

    script:
    """
    run_DESeq2.r --count_table1 ${count_table_01} --alias1 ${group_name_01}\
                 --count_table2 ${count_table_02} --alias2 ${group_name_02}\
                 --output_file ${group_name_01}_vs_${group_name_02}_DESeq2.tsv
    """
}
