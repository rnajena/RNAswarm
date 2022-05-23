import matplotlib.pyplot as plt
import itertools
import seaborn as sns
import numpy as np
import itertools
# from matplotlib_venn import venn3

def parse_genome(genome_file):
    """
    """
    genome_dict = {}

    header = ""
    seq = ""

    with open(genome_file) as file:  # read genome file
        for line in file:  # parse genome file
            if line.startswith(">"):  # parse header
            #* if there is a header already, we store the current sequence
            #* to this header.
                if header:
                    genome_dict[header] = seq
                    #* then we flush the sequence
                    seq = ""
                #* and parse the new header
                header = line.strip()[1:]
            else:
                #* if no header is detected, we append the new line
                #* to the current sequence
                seq += line.strip()
        #* after the last line, we have to store
        #* the last sequence as well. Since no new line with ">" occurs, we
        #* do this manually
        genome_dict[header] = seq
    return genome_dict

def make_combination_array(genome_dict):
    """
    Creates a dictionary of numpy array of all possible genome segment combinations.
    Use helper.parse_genome() to create genome_dict.
    """
    combination_arrays = {}
    segments = list(genome_dict.keys())

    # segment_combinations = [
    #     segment_combination
    #     for segment_combination in itertools.combinations_with_replacement(segments, 2)
    # ]
    # * while I usually appreciate the usage of list comprehensions, you can directly transform
    # * the iterator to a list. Actually, we also could just put the iterator in the for loop.
    # * should work as well. Is a tad more memory efficient.
    segment_combinations = list(itertools.combinations_with_replacement(segments, 2))

    for segment_combination in segment_combinations:
        # for segment_combination in itertools.combinations_with_replacement(segments,2): # * this should work as well
        combination_arrays[segment_combination] = np.zeros(
            (
                len(genome_dict[segment_combination[0]]),
                len(genome_dict[segment_combination[1]]),
            )
        )
    return combination_arrays

def plot_heatmap(array, output_dir, filename, combination_0, combination_1):
    heatmap = sns.heatmap(array, square=True, cmap="YlGnBu_r")
    heatmap.set(xlabel=str(combination_1), ylabel=str(combination_0))
    plt.figure()
    heatmap.figure.savefig(f"{output_dir}/{filename}", bbox_inches="tight")
    plt.close("all")

def get_ids(fastq_file):
    """
    get the ids of the reads in the fastq file
    """
    id_list = []
    with open(fastq_file, "r") as fastq:
        for line in fastq:
            if line.startswith("@"):
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


def get_unmapped_reads(bam_file):
    """
    get the ids of the unmapped reads in the bam file
    """
    unmapped_reads = set()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                unmapped_reads.add(read.query_name)
    return unmapped_reads


def get_chim_ids(chim_file):
    """
    get the ids of the chimeric reads in the chim file
    """
    chim_ids = set()
    with open(chim_file, "r") as chim_file:
        for line in chim_file:
            chim_ids.add(line.strip().split("\t")[0])
    return chim_ids


def get_trs_ids(trns_file):
    """
    get the ids of the chimeric reads in the trns file
    """
    trns_ids = set()
    with open(trns_file, "r") as trns_file:
        for line in trns_file:
            trns_ids.add(line.strip().split("\t")[-1])
    return trns_ids