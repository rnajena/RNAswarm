#!/usr/bin/env python3

import sys
import helper
import numpy as np
import seaborn as sns
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from os.path import splitext, basename

"""

"""

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
        print(highest_point)
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


def main():
    """Main function

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    genome_file_path = sys.argv[1]
    genome_dict = helper.parse_fasta(genome_file_path)
    combination_array = helper.make_combination_array(genome_dict)
    readsOfInterest = sys.argv[2]
    print(
        f"Genome file: {genome_file_path}\n",
        f"File used for interaction parsing: {readsOfInterest}",
    )
    if sys.argv[4] == "--segemehl_mode":
        segemehlTrans2heatmap(readsOfInterest, combination_array)
    elif sys.argv[4] == "--bwa_mode":
        bwaChimera2heatmap(readsOfInterest, combination_array)

    # Creating the diversity plot
    diversity_dict = get_diversity(combination_array)
    x_axis = list(diversity_dict.keys())
    y_axis = list(diversity_dict.values())
    diversity_plot = sns.lineplot(x=x_axis, y=y_axis)
    diversity_plot.set(yscale='log')
    plt.figure()
    diversity_plot.figure.savefig(f"{sys.argv[3]}/{splitext(splitext(basename(readsOfInterest))[0])[0]}_diversity.png", bbox_inches="tight")
    plt.close("all")

    for combination, array in combination_array.items():
        helper.plot_heatmap(
            array,
            sys.argv[3],
            f"{combination[0]}_{combination[1]}",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap_log(
            array,
            sys.argv[3],
            f"{combination[0]}_{combination[1]}_log",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap(
            array > np.power(10,1.5),
            sys.argv[3],
            f"{combination[0]}_{combination[1]}_cut10to1_5",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap(
            array > np.power(10,2),
            sys.argv[3],
            f"{combination[0]}_{combination[1]}_cut10to2",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap(
            array > np.power(10,2.5),
            sys.argv[3],
            f"{combination[0]}_{combination[1]}_cut10to2_5",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap(
            array > np.power(10,3),
            sys.argv[3],
            f"{combination[0]}_{combination[1]}_cut10to3",
            combination[0],
            combination[1],
        )
        helper.plot_heatmap(
            detect_peaks(array)[0],
            sys.argv[3],
            f"{combination[0]}_{combination[1]}_peaks",
            combination[0],
            combination[1],
        )


if __name__ == "__main__":
    main()
