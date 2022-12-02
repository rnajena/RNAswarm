#!/usr/bin/env python3

"""Make deseq2 table

Usage:
  make_deseq2_table.py -i <samples_csv> -s
  make_deseq2_table.py -i <samples_csv> -b

Options:
  -h --help                 Show this screen.
  -i --input=<samples_csv>  Input folder containing a csv table describing the
                            experiment schema files.
  -s --segemehl_mode        Segemehl mode is used with .txt.trns files as input.
  -b --bwa_mode             Bwa mode is used with .chim files as input.
"""
from docopt import docopt
import helper
import os
import handle_chimeras as hc
from scipy import ndimage
import numpy as np


def parse_samples_csv(csv_filepath):
    samples_dict = {}
    with open(csv_filepath) as inputStream:
        for line in inputStream:
            line = line.strip()
            trns_filepath = line.split(",")[0]
            genome_filepath = line.split(",")[1]
            samples_dict[trns_filepath] = genome_filepath
    return samples_dict


def make_matrices_dict(samples_dict):
    matrices_dict = {}
    for trns_filepath, genome_filepath in samples_dict.items():
        genome_dict = helper.parse_fasta(genome_filepath)
        combination_array = helper.make_combination_array(genome_dict)
        hc.segemehlTrans2heatmap(trns_filepath, combination_array)
        matrices_dict[trns_filepath] = combination_array
    return matrices_dict


def make_labeled_matrices_dict(matrices_dict):
    labeled_matrices_dict = {}
    for trns_filepath, array_dict in matrices_dict.items():
        for combination, array in array_dict.items():
            filtered_array = array > np.power(10, 1.5)
            labeled_array, num_features = ndimage.label(filtered_array)
            labeled_matrices_dict[trns_filepath] = {
                combination: {
                    "labeled_array": labeled_array,
                    "num_features": num_features,
                }
            }
    return labeled_matrices_dict


# def calculate_mean(matrices_dict, labeled_matrices_dict):
#    first_line = True
#    means_dict = {}
#    for lb_trns_filepath, lb_dict in labeled_matrices_dict.items():
#         for trns_filepath, array_dict in matrices_dict.items():
#             for lb_combination, lb_array in lb_dict.items():
#                 for combination, array in array_dict.items():
#                         print(f"{lb_trns_filepath},{lb_combination},{label},{}")
#                     labeled_array = lb_array["labeled_array"]
#                     num_features = lb_array["num_features"]
#                     if num_features > 0:
#                         means_array = ndimage.labeled_comprehension(
#                             array,
#                             labeled_array,
#                             np.arange(1, num_features + 1),
#                             np.mean,
#                             float,
#                             -1,
#                         )
#                         for label in np.arange(1, num_features + 1):
#                             interaction_id = f"{lb_trns_filepath}_{combination[0]}-{combination[1]}_{label}"
#                             means_dict[interaction_id] = {
#                                 trns_filepath: means_array[label - 1]
#                             }

#     return means_dict

def calculate_means(array_dict, labeled_array, num_labels):
    means_dict = {}
    for trns_file, array in array_dict.items():
        means_array = ndimage.labeled_comprehension(
                                    array,
                                    labeled_array,
                                    np.arange(1, num_labels + 1),
                                    np.mean,
                                    float,
                                    -1,
                                )
        means_dict[trns_file] = means_array
    return means_array

def print_deseq2_input_stream(matrices_dict, labeled_matrices_dict):
    

def main():
    arguments = docopt(__doc__)
    samples_dict = parse_samples_csv(arguments["--input"])
    matrices_dict = make_matrices_dict(samples_dict)
    labeled_matrices_dict = make_labeled_matrices_dict(matrices_dict)
    calculate_mean(matrices_dict, labeled_matrices_dict)
    print(means_dict)


if __name__ == "__main__":
    main()
