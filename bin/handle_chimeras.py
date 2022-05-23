#!/usr/bin/env python3

import sys
import helper
import numpy as np
import scipy.ndimage as ndimage

def __convert_to_int(element):
    if type(element) == int:
        return element
    elif element.isdigit():
        return int(element)
    else:
        return element


def __check_interaction(currentRow, interaction_arrays):
    """ """
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
    """ """
    with open(chimFile) as inputStream:
        for line in inputStream:
            currentRow = line.strip().split("\t")
            currentRow_sanitized = currentRow[1:7]
            interaction = __check_interaction(currentRow_sanitized, interaction_arrays)
            fill_heatmap(interaction, interaction_arrays)


def __extract_start_stop_segemehl(read):
    """ """
    seg = read[0]
    start = int(read[1])
    stop = start + int(read[4])
    return [seg, start, stop]


def segemehlTrans2heatmap(trnsFile, interaction_arrays):
    """ """
    with open(trnsFile) as inputStream:
        for line in inputStream:
            line = line.strip().split()
            firstRead = line[0].split(",")
            secondRead = line[1].split(",")
            currentRow = __extract_start_stop_segemehl(
                firstRead
            ) + __extract_start_stop_segemehl(secondRead)
            interaction = __check_interaction(currentRow, interaction_arrays)
            fill_heatmap(interaction, interaction_arrays)


def fill_heatmap(interaction, interaction_arrays):
    """ """
    firstSegment = interaction[0]
    secondSegment = interaction[3]
    interaction_arrays[(firstSegment, secondSegment)][
        interaction[1] : interaction[2], interaction[4] : interaction[5]
    ] += 1

def detect_peaks(interaction_array):
    """
    Takes an image and detect the peaks using the local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """

    # define a 30 by 30 neighborhood
    neighborhood = np.ones((30,30))

    # apply the local maximum filter; all pixel of maximal value 
    # in their neighborhood are set to 1
    # local_max = np.isclose(ndimage.maximum_filter(interaction_array, footprint=neighborhood), interaction_array)
    local_max = ndimage.maximum_filter(interaction_array, footprint=neighborhood)==interaction_array

    # local_max is a mask that contains the peaks we are 
    # looking for, but also the background.
    # In order to isolate the peaks we must constraint the peaks to the foreground.

    # we create the mask of the foreground
    foreground = (interaction_array > 10)

    # we obtain the final mask, containing only peaks, 
    # by removing the background from the local_max mask (xor operation)
    detected_peaks = local_max & foreground

    labeled_array, num_features = ndimage.label(detected_peaks)

    detected_peaks = ndimage.maximum_filter(detected_peaks , footprint=neighborhood)

    return detected_peaks, labeled_array, num_features

def main():
    # * see below. I put these lines here as they do not
    # * change, no matter the if result.
    genome_file_path = sys.argv[1]
    genome_dict = helper.parse_genome(genome_file_path)
    combination_array = helper.make_combination_array(genome_dict)
    readsOfInterest = sys.argv[2]

    # ! We can use argparse (or docopt, but thats an extra library)
    # ! to handle the input parsing and optional parameter better
    if sys.argv[4] == "--segemehl_mode":
        # ! This is redundant. It does not matter whether
        # ! argv[4] is segemehl oder bwa, you are reading the genome anyway.
        #! I moved it in front of the if condition
        segemehlTrans2heatmap(readsOfInterest, combination_array)
        # trns_dict = parse_trns_file(readsOfInterest)
        # trns_dict_to_combination_array(interaction_arrays, trns_dict)
    elif sys.argv[4] == "--bwa_mode":
        bwaChimera2heatmap(readsOfInterest, combination_array)
        # chim_dict = parse_chim_file(readsOfInterest)
        # chim_dict_to_combination_array(interaction_arrays, chim_dict)

    # * same as above. the for loop is the same
    # * for both if conditions. So, it can be outside the if clause
    for combination, array in combination_array.items():
        helper.plot_heatmap(
            array,
            sys.argv[3],
            f"{combination[0]}_{combination[1]}",
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

# * general remark / notes:
# * It feels like we are doing way to much parsing here.
# * In the end, we just need to know which reads maps to what segment
# * and where it mapped. The alignment length and everything can be calculated that way.
# * By parsing every field from the trns.splice.txt and chimera.txt respectively,
# * your dictionaries are bloated as hell which, in turn, leads to the spaghetti code.
