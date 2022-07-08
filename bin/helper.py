import matplotlib.pyplot as plt
import itertools
import seaborn as sns
import numpy as np
import itertools
from matplotlib.colors import LogNorm
# from matplotlib_venn import venn3

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
            #* if there is a header already, we store the current sequence
            #* to this header.
                if header:
                    fasta_dict[header] = seq
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
        fasta_dict[header] = seq
    return fasta_dict

def make_combination_array(genome_dict):
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
    """
    Plots a heatmap of the given array.
    
    Parameters
    ----------
    array : numpy array
        Array to be plotted.
    output_dir : str
        Path to output directory.
    filename : str
        Name of output file.
    combination_0 : str
        First segment of combination.
    combination_1 : str
        Second segment of combination.

    Output
    ------
    None

    """
    heatmap = sns.heatmap(array, square=True, vmin=0, cmap="YlGnBu_r")
    heatmap.set(xlabel=str(combination_1), ylabel=str(combination_0))
    plt.figure()
    heatmap.figure.savefig(f"{output_dir}/{filename}", bbox_inches="tight")
    plt.close("all")

def plot_heatmap_log(array, output_dir, filename, combination_0, combination_1):
    """
    Plots a heatmap of the given array in log scale.

    Parameters
    ----------
    array : numpy array
        Array to be plotted.
    output_dir : str
        Path to output directory.
    filename : str
        Name of output file.
    combination_0 : str
        First segment of combination.
    combination_1 : str
        Second segment of combination.

    Output
    ------
    None

    """
    heatmap = sns.heatmap(array, square=True, vmin=0, cmap="PiYG", norm=LogNorm())
    heatmap.set(xlabel=str(combination_1), ylabel=str(combination_0))
    plt.figure()
    heatmap.figure.savefig(f"{output_dir}/{filename}", bbox_inches="tight")
    plt.close("all")

# def plot_venn_diag(sets_dict):
#     pass