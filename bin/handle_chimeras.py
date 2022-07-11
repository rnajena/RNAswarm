#!/usr/bin/env python3

"""handle_chimeras.py

Usage:
  handle_chimeras.py -g <genome> -i <input_file> -o <output_folder> --bwa_mode
  handle_chimeras.py -g <genome> -i <input_file> -o <output_folder> --segemehl_mode

Options:
  -h --help                         Show this screen.
  --segemehl_mode                   Use this mode if you want to use a trns.txt file
                                    outputted by segemehl.
  --bwa_mode                        Use this mode if you want to use a chim.txt file
                                    outputted by bwa.
  -g --genome=<genome>              The genome filepath.
  -i --input=<input_file>           The input filepath, either a chim.txt or trns.txt
                                    file, depending on the mode.
  -o --output=<output_folder>       The output folder.

"""


from docopt import docopt
from os.path import splitext, basename
import sys
import itertools
import helper
import numpy as np
import seaborn as sns
import scipy.ndimage as ndimage
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

    if (interaction[0], interaction[3]) not in interaction_arrays:
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


def __extract_length_segemehl(read):
    """Returns the length of the read

    Parameters
    ----------
    read : list

    Returns
    -------
    int
    """
    return int(read[4])


def segemehlTrans2heatmap(trnsFile, interaction_arrays):
    """Parses the trns file and fills the interaction_arrays

    Parameters
    ----------
    trnsFile : str
    interaction_arrays : dict

    Returns
    -------
    None
    """
    total_mapped_interactions = 0
    sum_interactionArea = 0
    sum_all_matrices = 0
    with open(trnsFile) as inputStream:
        for line in inputStream:
            line = line.strip().split()
            firstRead = line[0].split(",")
            secondRead = line[1].split(",")
            currentRow = __extract_start_stop_segemehl(
                firstRead
            ) + __extract_start_stop_segemehl(secondRead)
            interactionArea = __extract_length_segemehl(
                firstRead
            ) * __extract_length_segemehl(secondRead)
            interaction = __check_interaction(currentRow, interaction_arrays)
            fill_heatmap(interaction, interaction_arrays)
            total_mapped_interactions += 1
            sum_interactionArea += interactionArea
    for array in interaction_arrays.values():
        sum_all_matrices += array.sum()
    print(f"Number of interactions added to the matrices: {total_mapped_interactions}")
    print(f"Interaction area on the trns.txt file: {sum_interactionArea}")
    print(f"Interaction area on the matrices: {sum_all_matrices}")


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


def get_diversity(interaction_arrays):
    """Returns the diversity of the interaction_arrays

    Parameters
    ----------
    interaction_arrays : dict

    Returns
    -------
    dict
    """
    diversity_dict = {}
    for combination, interaction_array in interaction_arrays.items():
        highest_point = int(np.nanmax(interaction_array))
        for i in range(0, highest_point):
            if i in diversity_dict.keys():
                diversity_dict[i] = diversity_dict[i] + (interaction_array == i).sum()
            else:
                diversity_dict[i] = (interaction_array == i).sum()

    return diversity_dict


def detect_peaks(interaction_array):
    """
    Takes an image and detect the peaks using the local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)

    Parameters
    ----------
    interaction_array : ndarray
        Input array.

    Returns
    -------
    boolean ndarray
        A boolean mask of the peaks in the array.

    ndarray
        An array of the same shape as input array. Each peak has an unique value.

    int
        The number of peaks detected.
    """

    # define a 30 by 30 neighborhood
    neighborhood = np.ones((30, 30))

    # apply the local maximum filter; all pixel of maximal value
    # in their neighborhood are set to 1
    local_max = (
        ndimage.maximum_filter(interaction_array, footprint=neighborhood)
        == interaction_array
    )

    # local_max is a mask that contains the peaks we are
    # looking for, but also the background.
    # In order to isolate the peaks we must constraint the peaks to the foreground.

    # we create the mask of the foreground
    foreground = interaction_array > 10

    # we obtain the final mask, containing only peaks,
    # by removing the background from the local_max mask (xor operation)
    detected_peaks = local_max & foreground

    labeled_array, num_features = ndimage.label(detected_peaks)

    detected_peaks = ndimage.maximum_filter(detected_peaks, footprint=neighborhood)

    return detected_peaks, labeled_array, num_features


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
            pairwise_arrays[segment_combination] = np.sum(interaction_arrays[segment_combination], axis=1)
        elif (segment_combination[1], segment_combination[0]) in interaction_arrays.keys():
            pairwise_arrays[segment_combination] = np.sum(interaction_arrays[(segment_combination[1], segment_combination[0])], axis=0)
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
        plt.legend(loc='best')
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
    if arguments["--segemehl_mode"]:
        segemehlTrans2heatmap(readsOfInterest, combination_array)
    elif arguments["--bwa_mode"]:
        bwaChimera2heatmap(readsOfInterest, combination_array)

    # Plotting the pairwise arrays
    plot_pairwise_arrays(combination_array, genome_dict, output_folder) 

    # Creating the diversity plot
    diversity_dict = get_diversity(combination_array)
    x_axis = list(diversity_dict.keys())
    y_axis = list(diversity_dict.values())
    diversity_plot = sns.lineplot(x=x_axis, y=y_axis)
    diversity_plot.set(yscale="log")
    plt.figure()
    diversity_plot.figure.savefig(
        f"{output_folder}/{splitext(splitext(basename(readsOfInterest))[0])[0]}_diversity.png",
        bbox_inches="tight",
    )
    plt.close("all")

    # Plotting heatmaps
    for combination, array in combination_array.items():
        helper.plot_heatmap(
            array,
            output_folder,
            f"{combination[0]}_{combination[1]}",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap_log(
            array,
            output_folder,
            f"{combination[0]}_{combination[1]}_log",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap(
            array > np.power(10, 1.5),
            output_folder,
            f"{combination[0]}_{combination[1]}_cut10to1_5",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap(
            array > np.power(10, 2),
            output_folder,
            f"{combination[0]}_{combination[1]}_cut10to2",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap(
            array > np.power(10, 2.5),
            output_folder,
            f"{combination[0]}_{combination[1]}_cut10to2_5",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap(
            array > np.power(10, 3),
            output_folder,
            f"{combination[0]}_{combination[1]}_cut10to3",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap(
            detect_peaks(array)[0],
            output_folder,
            f"{combination[0]}_{combination[1]}_peaks",
            combination[0],
            combination[1],
        )


if __name__ == "__main__":
    main()
