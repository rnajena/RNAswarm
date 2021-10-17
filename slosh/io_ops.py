import numpy as np
import os
import re
from scipy import stats

import visualize_data as vd


def regex_from_segments(viral_segments):
    """
    Returns a compiled regex that can be used to search for pairs of viral segment
    abreviations in a given string.

    Precondition:
    The viral segments must be writen as sgmt01_segmt02 in order to be found.
    """
    segments_to_regex = ""
    for viral_segment01 in viral_segments:
        for viral_segment02 in viral_segments:
            segments_to_regex = (
                segments_to_regex + f"{viral_segment01}_{viral_segment02}|"
            )
    compiled_regex = re.compile(segments_to_regex[:-1])
    return compiled_regex


def read_arrays(data_directory, viral_segments):
    """
    Takes a directory containing subdirectories related to diferent replicates and
    returns two dictionaries, one describing the paths to the arrays and another
    containing the arrays for each segment combination.

    Input:

    Return:

    """
    d_repDir2Combinations = {}
    d_combination2Array = {}
    segments_to_regex = regex_from_segments(viral_segments)
    replicateDirs = [
        entry
        for entry in os.listdir(data_directory)
        if os.path.isdir(f"{data_directory}/{entry}")
    ]

    for repDir in replicateDirs:
        allCombinations = os.listdir(f"{data_directory}/{repDir}")
        d_repDir2Combinations[repDir] = allCombinations

        for combination in allCombinations:
            uniqueID = f"{segments_to_regex.search(combination).group()}"
            array = np.load(f"{data_directory}/{repDir}/{combination}")
            if uniqueID in d_combination2Array:
                d_combination2Array[uniqueID].append(array)
            else:
                d_combination2Array[uniqueID] = [array]

    return (d_repDir2Combinations, d_combination2Array)


def calculate_variances(d_arrays):
    """
    Takes a dictionary with numpy arrays and calculates the variance for each unique
    entry.

    Input:
        d_arrays -- {'uniqueName' : [np.arrays]}

    Return:
        d_comb2variance -- {'uniqueName' : np.array}

    """
    d_comb2variance = {}
    for combination, l_countTable in d_arrays.items():
        d_comb2variance[combination] = np.var(l_countTable, axis=0)

    return d_comb2variance


def print_csv_from_dict(dict):
    for comb, characteristics in dict.items():
        for characteristic, value in characteristics.items():
            print("{},{},".format(characteristic, value), end="")


def save_heatmaps(variances, output_dir):
    """
    Takes a dictionary of np.arrays and plots heatmaps for the np.arrays and saves them
     to the output_dir

    Input:

    Return:
    """
    # Retrieves the key (comb) and the value (variance_array)
    for comb, variance_array in variances.items():
        vd.plot_heatmap(variance_array, output_dir, f"{comb}_hist")


def save_histograms(variances, output_dir):
    for comb, variance_array in variances.items():
        vd.plot_histogram(variance_array, f"{output_dir}/{comb}_hist")


def save_characteristics(variances, output_dir):
    d_comb2characteristics = {}

    for comb, variance_array in variances.items():
        statistic, pvalue = stats.normaltest(variance_array.flatten())
        d_comb2characteristics[comb] = {
            "normaltest_statistic": statistic,
            "normaltest_pvalue": pvalue,
            "minimun": np.amin(variance_array),
            "maximun": np.amin(variance_array),
            "range": np.ptp(variance_array),
        }


def format_means_to_table(readcounts, sep=",", output_path=None):
    """Returns teste teste

    Parameters
    ----------
    readcounts: ...
    path: ...

    Returns
    -------
    res: DEseq2_input
    """
    table = ""
    # we should check if all dicts inside readcounts have the same size

    for idx in readcounts[0].keys():
        for sample_id in range(len(readcounts)):
            if sample_id == 0:
                table = f"{table}{idx}{sep}{int(round(readcounts[sample_id][idx]))}"
            elif sample_id == len(readcounts) - 1:
                table = f"{table}{sep}{int(round(readcounts[sample_id][idx]))}\n"
            else:
                table = f"{table}{sep}{int(round(readcounts[1][idx]))}"

    assert sum("\n" in char for char in table) == len(
        readcounts[0].keys()
    ), "Size of table is not consistent with dictionary"

    if output_path:
        with open(output_path, "w") as file:
            file.write(table)

    return table


def format_arrays_to_table():
    """Returns ...

    Parameters
    ----------
    readcounts: ...
    path: ...

    Returns
    -------
    res:
    """
