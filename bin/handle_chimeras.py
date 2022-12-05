#!/usr/bin/env python3

"""handle_chimeras.py

Usage:
  handle_chimeras.py -g <genome> -i <input_file> -o <output_folder>

Options:
  -h --help                         Show this screen.
  -g --genome=<genome>              The genome filepath.
  -i --input=<input_file>           The input filepath, either a chim.txt or trns.txt
                                    file, depending on the mode.
  -o --output=<output_folder>       The output folder.

"""


from docopt import docopt
import itertools
import helper
import numpy as np
import matplotlib.pyplot as plt


def __convert_to_int(element):
    """Return the integer value of the element

    Parameters
    ----------
    element : str

    Returns
    -------
    int
    """
    if type(element) == int:
        return element
    elif element.isdigit():
        return int(element)
    else:
        return element


def __check_interaction(currentRow, interaction_arrays):
    """Parses the current row and returns the interaction

    Parameters
    ----------
    currentRow : str
    interaction_arrays : dict

    Returns
    -------
    list
    """
    interaction = []
    currentRow = list(map(__convert_to_int, currentRow))
    if currentRow[1] > currentRow[2]:
        interaction += [currentRow[0]] + [currentRow[2]] + [currentRow[1]]
    else:
        interaction += currentRow[:3]

    if currentRow[4] > currentRow[5]:
        interaction += [currentRow[3]] + [currentRow[5]] + [currentRow[4]]
    else:
        interaction += currentRow[3:]
    # Checks if the cobination is not in the interaction_arrays
    if (interaction[0], interaction[3]) not in interaction_arrays:
        # If not, it reverses the interaction
        interaction = interaction[3:] + interaction[0:3]

    return interaction


def bwaChimera2heatmap(chimFile, interaction_arrays):
    """Parses the chim file and fills the interaction_arrays

    Parameters
    ----------
    chimFile : str
    interaction_arrays : dict

    Returns
    -------
    None
    """
    with open(chimFile) as inputStream:
        for line in inputStream:
            currentRow = line.strip().split("\t")
            currentRow_sanitized = currentRow[1:7]
            interaction = __check_interaction(currentRow_sanitized, interaction_arrays)
            fill_heatmap(interaction, interaction_arrays)


def __extract_start_stop_segemehl(read):
    """Returns the start and stop positions of the read

    Parameters
    ----------
    read : str

    Returns
    -------
    list
    """
    seg = read[0]
    start = int(read[1])
    stop = start + int(read[4])
    return [seg, start, stop]


def segemehlTrans2heatmap(trnsFile, interaction_arrays, intra_combinations=False):
    """Parses the trns file and fills the interaction_arrays

    Parameters
    ----------
    trnsFile : str
    interaction_arrays : dict

    Returns
    -------
    None
    """
    with open(trnsFile) as inputStream:
        for line in inputStream:
            line = line.strip().split()
            firstRead = line[0].split(",")
            secondRead = line[1].split(",")
            currentRow = __extract_start_stop_segemehl(
                firstRead
            ) + __extract_start_stop_segemehl(secondRead)
            interaction = __check_interaction(currentRow, interaction_arrays)
            if intra_combinations:
                fill_heatmap(interaction, interaction_arrays)
            else:
                if interaction[0] != interaction[3]:
                    fill_heatmap(interaction, interaction_arrays)


def fill_heatmap(interaction, interaction_arrays):
    """Fills the interaction_arrays with the interaction

    Parameters
    ----------
    interaction : list
    interaction_arrays : dict

    Returns
    -------
    int
    """
    firstSegment = interaction[0]
    secondSegment = interaction[3]
    interaction_arrays[(firstSegment, secondSegment)][
        interaction[1] : interaction[2], interaction[4] : interaction[5]
    ] += 1
    return 1


def get_histogram_dict(interaction_arrays):
    """Returns the diversity of the interaction_arrays

    Parameters
    ----------
    interaction_arrays : dict

    Returns
    -------
    dict
    """
    histogram_dict = {}
    for combination, interaction_array in interaction_arrays.items():
        highest_point = int(np.nanmax(interaction_array))
        for i in range(0, highest_point):
            if i in histogram_dict.keys():
                histogram_dict[i] = histogram_dict[i] + (interaction_array == i).sum()
            else:
                histogram_dict[i] = (interaction_array == i).sum()

    return histogram_dict


def get_pairwise_arrays(interaction_arrays, genome_dict):
    """Returns the pairwise arrays of the interaction_arrays

    Parameters
    ----------
    interaction_arrays : dict

    Returns
    -------
    dict
    """
    pairwise_arrays = {}

    for segment_combination in itertools.permutations(genome_dict.keys(), 2):
        if segment_combination in interaction_arrays.keys():
            pairwise_arrays[segment_combination] = np.sum(
                interaction_arrays[segment_combination], axis=1
            )
        elif (
            segment_combination[1],
            segment_combination[0],
        ) in interaction_arrays.keys():
            pairwise_arrays[segment_combination] = np.sum(
                interaction_arrays[(segment_combination[1], segment_combination[0])],
                axis=0,
            )
        else:
            raise KeyError(f"{segment_combination} not in the interaction_arrays")
    return pairwise_arrays


def plot_pairwise_arrays(interaction_arrays, genome_dict, foldername):
    """Plots the pairwise arrays of the interaction_arrays

    Parameters
    ----------
    interaction_arrays : dict
    genome_dict : dict
    folderpath : str

    Returns
    -------
    None
    """
    pairwise_arrays = get_pairwise_arrays(interaction_arrays, genome_dict)

    # Plotting the pairwise arrays
    for segment in genome_dict.keys():
        for segment_combination, pairwise_array in pairwise_arrays.items():
            if segment == segment_combination[0]:
                filepath = f"{foldername}/{segment}_pairwise.png"
                plt.plot(pairwise_array, label=segment_combination[1])
        plt.title(f"{segment}")
        plt.legend(loc="best")
        plt.savefig(filepath)
        plt.close("all")


def main():
    """Main function

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # Parse the command line arguments
    arguments = docopt(__doc__)
    genome_file_path = arguments["--genome"]
    readsOfInterest = arguments["--input"]
    output_folder = arguments["--output"]

    # Process input files
    genome_dict = helper.parse_fasta(genome_file_path)
    combination_array = helper.make_combination_array(genome_dict)
    print(
        f"Genome file: {genome_file_path}\n",
        f"File used for interaction parsing: {readsOfInterest}",
    )
    segemehlTrans2heatmap(readsOfInterest, combination_array)


if __name__ == "__main__":
    main()
