/*************************************************************************
* annotation table to fasta
**************************************************************************/
process annotationTableToFasta {
    label: "RNAswarm_small"

    input:
    tuple
}


/*************************************************************************
* predict RNA-RNA structures
**************************************************************************/
process runRNAcofold {
    label: "viennarna"

    input:
    tuple 
}