import matplotlib.pyplot as plt
import itertools
import numpy as np
import itertools
from matplotlib.colors import LogNorm
import pandas as pd


def parse_fasta(fasta_file):
    """Parses a fasta file and returns a dictionary of segment names and sequences.

    Parameters
    ----------
    fasta_file : str
        Path to fasta file.

    Returns
    -------
    genome_dict : dict
        Dictionary of segment names and sequences.

    """
    fasta_dict = {}

    header = ""
    seq = ""

    with open(fasta_file) as file:  # read genome file
        for line in file:  # parse genome file
            if line.startswith(">"):  # parse header
                # * if there is a header already, we store the current sequence
                # * to this header.
                if header:
                    fasta_dict[header] = seq
                    # * then we flush the sequence
                    seq = ""
                # * and parse the new header
                header = line.strip()[1:]
            else:
                # * if no header is detected, we append the new line
                # * to the current sequence
                seq += line.strip()
        # * after the last line, we have to store
        # * the last sequence as well. Since no new line with ">" occurs, we
        # * do this manually
        fasta_dict[header] = seq
    return fasta_dict


def make_combination_array(genome_dict, intra_only=False):
    """
    Creates a dictionary of numpy array of all possible genome segment combinations.
    Use helper.parse_genome() to create genome_dict.

    Parameters
    ----------
    genome_dict : dict
        Dictionary of segment names and sequences.

    Returns
    -------
    combination_arrays : dict
        Dictionary of numpy arrays of all relevant genome segment combinations.

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
    if intra_only:
        segment_combinations = list(
            itertools.combinations_with_replacement(segments, 2)
        )
        segment_combinations = [
            segment_combination
            for segment_combination in segment_combinations
            if segment_combination[0] == segment_combination[1]
        ]
    else:
        segment_combinations = list(
            itertools.combinations_with_replacement(segments, 2)
        )
        segment_combinations = [
            segment_combination
            for segment_combination in segment_combinations
            if segment_combination[0] != segment_combination[1]
        ]

    for segment_combination in segment_combinations:
        # for segment_combination in itertools.combinations_with_replacement(segments,2): # * this should work as well
        combination_arrays[segment_combination] = np.zeros(
            (
                len(genome_dict[segment_combination[0]]),
                len(genome_dict[segment_combination[1]]),
            )
        )
    return combination_arrays


def positive_to_negative_strand_point(genome_dict, segment, position):
    """Transpose a point from the positive to the negative strand.

    Args:
        genome_dict (dict): A dictionary containing the genome.
        aSeq (str): The name of the segment.
        a (int): The position of the point on the segment.

    Returns:
        int:
            The position of the point on the negative strand.
    """
    # Get the length of the segment
    aLen = len(genome_dict[segment])
    # Get the position on the negative strand
    return aLen - position + 1


def negative_to_positive_strand(genome_dict, aSeq, cai, caj, bSeq, cbi, cbj):
    """Transpose the region of interest from the negative to the positive strand.

    Args:
        genome_dict (dict): A dictionary containing the genome.
        aSeq (str): The name of the first segment.
        bSeq (str): The name of the second segment.
        cai (int): The start position of the first segment.
        caj (int): The end position of the first segment.
        cbi (int): The start position of the second segment.
        cbj (int): The end position of the second segment.

    Returns:
        tuple: The transposed region of interest.
    """
    # Get the length of the first segment
    aLen = len(genome_dict[aSeq])

    # Get the length of the second segment
    bLen = len(genome_dict[bSeq])

    # Transpose the first segment
    caj_pos = aLen - cai + 1
    cai_pos = aLen - caj

    # Transpose the second segment
    cbj_pos = bLen - cbi + 1
    cbi_pos = bLen - cbj

    return aSeq, cai_pos, caj_pos, bSeq, cbi_pos, cbj_pos


def positive_to_negative_strand(genome_dict, aSeq, cai, caj, bSeq, cbi, cbj):
    """Transpose the region of interest from the positive to the negative strand.

    Args:
        genome_dict (dict): A dictionary containing the genome.
        aSeq (str): The name of the first segment.
        bSeq (str): The name of the second segment.
        cai (int): The start position of the first segment.
        caj (int): The end position of the first segment.
        cbi (int): The start position of the second segment.
        cbj (int): The end position of the second segment.

    Returns:
        tuple: The transposed region of interest.
    """
    # Get the length of the first segment
    aLen = len(genome_dict[aSeq])

    # Get the length of the second segment
    bLen = len(genome_dict[bSeq])

    # Transpose the first segment
    caj_neg = aLen - cai
    cai_neg = aLen + 1 - caj

    # Transpose the second segment
    cbj_neg = bLen - cbi
    cbi_neg = bLen + 1 - cbj

    return aSeq, cai_neg, caj_neg, bSeq, cbi_neg, cbj_neg


def parse_annotation_table(annotation_table):
    """
    Parse annotation table

    Parameters
    ----------
    annotation_table : str
        Path to annotation table as a csv file or a xlsx file

    Returns
    -------
    regions : pandas.DataFrame
        A table containing the annotations for the rectangular regions.
    """
    # read the annotation table
    if annotation_table.lower().endswith(".xlsx"):
        regions = pd.read_excel(annotation_table)
    elif annotation_table.lower().endswith(".csv"):
        # add header
        header = [
            "segment01",
            "start01",
            "end01",
            "segment02",
            "start02",
            "end02",
        ]
        regions = pd.read_csv(annotation_table, names=header)
        # add id column
        regions["id"] = regions.index
    elif annotation_table.lower().endswith(".tsv"):
        regions = pd.read_csv(annotation_table, sep="\t")
    else:
        raise ValueError("Annotation table must be either a csv file or a xlsx file.")
    return regions
