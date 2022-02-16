
params.read = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/reads/SC35M_WT/SC35M_WT_repli01_0120.fastq'
params.genome = '/home/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/genomes/SC35M_WT.fasta'
params.output = '/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/SC35M_WT/'

genome_ch = Channel.fromPath(params.genome)

/************************************************************************
* sgemehl INDEX
************************************************************************/
process segemehlIndex {
    label 'segemehl'

    cpus 8
    memory '64 GB' 
    executor 'slurm'
    conda '../envs/segemehl.yaml'

    input:
    file genome from genome_ch

    output:
    file "SC35M_WT.idx" into genome_idx_ch

    script:
    """
    segemehl.x -x SC35M_WT.idx -d $genome
    """
}
