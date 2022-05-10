#!/usr/bin/env python3

"""table maker.

Usage:
  art_templater.py -i <interaction_reads> -g <genomic_reads> -b <bam>

Options:
  -h --help                              Show this screen.
  -i --interactions=<interaction_reads>  fastq file with interaction reads 
  -f --fasta=<genomic_reads>             fastq file with genomic reads
  -b --bam=<bam_file>                    mapped reads in bam format
"""
from docopt import docopt
import pysam


def check_mapper(bam_file):
    """
    check if the bam file is mapped with bwa-mem or segemehl
    """
    pass

def get_ids(fastq_file):
    """
    get the ids of the reads in the fastq file
    """
    id_list = []
    with open(fastq_file, 'r') as fastq:
        for line in fastq:
            if line.startswith('@'):
                id_list.append(line.split()[0][1:])
    return id_list

def get_yz_flags(bam_file):
    """
    get the YZ:Z: flags of the reads in the bam file
    Args:
        bam_file (str): Path to BAM file.

    Returns:
        dict: Dictionary of flags, with read id as key and YZ:Z: flag as value.
    """
    yz_tags_dict = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            for tag in read.tags:
                if tag[0] == "YZ":
                    yz_tags_dict[read.query_name] = int(tag[1])
    return yz_tags_dict

def get_mapped_reads(bam_file):
    """
    get the ids of the mapped reads in the bam file
    """
    mapped_reads = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_mapped:
                mapped_reads.append(read.query_name)

def get_unmapped_reads(bam_file):
    """
    get the ids of the unmapped reads in the bam file
    """
    unmapped_reads = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                unmapped_reads.append(read.query_name)

def get_chim_ids(chim_file):
    """
    get the ids of the chimeric reads in the chim file
    """
    chim_ids = []
    with open(chim_file, 'r') as chim_file:
        for line in chim_file:
            chim_ids.append(line.line.strip().split("\t")[0])

def make_confusion_matrix_segemehl(ids_genome, ids_interactions, flags_dict):
    """
    make a confusion matrix from the ids and flags
    """
    confusion_matrix = {"true_negative" : 0,
                        "false_positive" : 0,
                        "false_negative" : 0,
                        "true_positive" : 0,
                        "unmapped_genome" : 0,
                        "unmapped_interaction" : 0,
                        "total" : 0}

    for id in ids_genome:
        if id not in flags_dict.keys():
            confusion_matrix["unmapped_genome"] += 1
        elif flags_dict[id] != 8:
                confusion_matrix["true_negative"] += 1
        elif flags_dict[id] == 8:
            confusion_matrix["false_positive"] += 1
    for id in ids_interactions:
        if id not in flags_dict.keys():
            confusion_matrix["unmapped_interaction"] += 1
        elif flags_dict[id] != 8:
            confusion_matrix["false_negative"] += 1
        elif flags_dict[id] == 8:
            confusion_matrix["true_positive"] += 1
    confusion_matrix["total"] = confusion_matrix["true_negative"] + confusion_matrix["false_positive"] + confusion_matrix["false_negative"] + confusion_matrix["true_positive"] + confusion_matrix["unmapped"]
    return confusion_matrix

def make_confusion_matrix_bwa(ids_genome, ids_interactions, chim_flags, bam_file):
    """
    make a confusion matrix from the ids and flags
    """
    confusion_matrix = {"true_negative" : 0,
                        "false_positive" : 0,
                        "false_negative" : 0,
                        "true_positive" : 0,
                        "unmapped_genome" : 0,
                        "unmapped_interaction" : 0,
                        "total" : 0}
                        
    unmapped = get_unmapped_reads(bam_file)
    mapped = get_mapped_reads(bam_file)

    for id in ids_genome:
        if id in unmapped:
            confusion_matrix["unmapped_genome"] += 1
        elif id not in chim_flags and id in mapped:
                confusion_matrix["true_negative"] += 1
        elif id in chim_flags:
            confusion_matrix["false_positive"] += 1
    for id in ids_interactions:
        if id in unmapped:
            confusion_matrix["unmapped_interaction"] += 1
        elif id in chim_flags:
            confusion_matrix["true_positive"] += 1
        elif id not in chim_flags and id in mapped:
            confusion_matrix["false_negative"] += 1
        
    confusion_matrix["total"] = confusion_matrix["true_negative"] + confusion_matrix["false_positive"] + confusion_matrix["false_negative"] + confusion_matrix["true_positive"] + confusion_matrix["unmapped"]
    return confusion_matrix

def main():
    #arguments = docopt(__doc__)
    genome_reads_file = "/mnt/local/ru27wav/Projects/RNAswarm/test/results/0.1_percent_trans/00-simulated_reads/pr8.fq"
    interaction_reads_file = "/mnt/local/ru27wav/Projects/RNAswarm/test/results/0.1_percent_trans/00-simulated_reads/pr8_interactions.fq"
    segemehl_bam_file = "/mnt/local/ru27wav/Projects/RNAswarm/test/results/0.1_percent_trans/02-mappings/segemehl/pr8_concat_trimmed_segemehl.bam"
    genome_ids = get_ids(genome_reads_file)
    interaction_ids = get_ids(interaction_reads_file)
    mapping_flags = get_yz_flags(segemehl_bam_file)
    conf_mx = make_confusion_matrix(genome_ids, interaction_ids, mapping_flags)
    print(conf_mx)

if __name__ == "__main__":
    main()