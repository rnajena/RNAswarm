process makeDeseq2Table {
    label 'python3'

    input:
    path(files)

    output:
    tuple path('samples.csv'), path('deseq2_input.csv')

    publishDir "${params.output}/05-differential_analysis/", mode: 'copy'

    script:
    """
    make_deseq2_table.py -i samples.csv -s > deseq2_input.csv
    """
}