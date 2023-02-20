#!/usr/bin/env python3

"""parse_peaks.py

This script takes an annotation table, a genome file (fasta format) and one or more trns files, and returns a table with the annotation and the peak base pair for each interaction.
The annotation table must have the following columns: id,segment01,start01,end01,segment02,start02,end02

Usage:
    parse_peaks.py <input_file> <input_file>... -a <annotation_table> -g <genome> -o <output_file>
    parse_peaks.py <input_file> -a <annotation_table> -g <genome> -o <output_file>
    parse_peaks.py -h | --help

Options:
    -h --help                                 Show this screen.
    <input_file>                              The input files to process, 
                                              has to be a trns file generated by segemehl.
    -a --annotation_table=<annotation_table>  The annotation table filepath.
    -g --genome=<genome>                      The genome filepath.
    -o --output=<output_file>                 The output directory.
"""

from docopt import docopt
import os
import numpy as np
import pandas as pd
import helper as hp
import trns_handler as th
import array_handler as ah


def get_peak_cell_from_annotation_table(combination_arrays, annotation):
    """
    Get the peak cell of the interaction matrix from the annotation table.

    Parameters
    ----------
    combination_arrays : dict
        A dictionary of arrays, with the keys being the combination of segments.

    annotation : pandas.DataFrame
        The annotation data frame.

    Returns
    -------
    dict
        A dictionary of peak cells, with the keys being the combination of segments.
    """
    segment_combination = (annotation["segment01"], annotation["segment02"])
    segment_1_start = int(annotation["start01"])
    segment_1_end = int(annotation["end01"])
    segment_2_start = int(annotation["start02"])
    segment_2_end = int(annotation["end02"])

    if (
        segment_combination not in combination_arrays
        and (segment_combination[::-1]) not in combination_arrays
    ):
        raise ValueError("Combination not found")

    if segment_combination in combination_arrays:
        combination_array = combination_arrays[segment_combination]
        region_of_interest = combination_array[
            segment_1_start:segment_1_end, segment_2_start:segment_2_end
        ]
        segment1_peak, segment2_peak = np.unravel_index(
            np.argmax(region_of_interest), region_of_interest.shape
        )
    else:
        segment_combination = segment_combination[::-1]
        combination_array = combination_arrays[segment_combination]
        region_of_interest = combination_array[
            segment_2_start:segment_2_end, segment_1_start:segment_1_end
        ]
        segment2_peak, segment1_peak = np.unravel_index(
            np.argmax(region_of_interest), region_of_interest.shape
        )

    peak_cell_value = np.max(region_of_interest)
    peak_dict = {
        "segment01_peak": segment1_peak + segment_1_start,
        "segment02_peak": segment2_peak + segment_2_start,
        "value_peak": peak_cell_value,
    }
    if (
        peak_cell_value
        != combination_array[peak_dict["segment01_peak"], peak_dict["segment02_peak"]]
    ):
        raise ValueError(
            "The value of the peak cell is not the same as the value of the cell in the peak_cell_coordinates"
        )
    return peak_dict


def main():
    args = docopt(__doc__)
    input_files = args["<input_file>"]
    annotation_table = args["--annotation_table"]
    genome = args["--genome"]

    # Parse annotation table
    # Check if the annotation table is a csv, tsv or excel file
    # if annotation_table.endswith(".csv"):
    #     annotation_table = pd.read_csv(annotation_table, sep=",")
    # elif annotation_table.endswith(".tsv"):
    #     annotation_table = pd.read_csv(annotation_table, sep="\t")
    if annotation_table.endswith(".XLSX") or annotation_table.endswith(".xlsx"):
        annotation_table = pd.read_excel(
            annotation_table,
            index_col=0,
            header=0,
            usecols=lambda x: "Unnamed" not in x,
        )
    else:
        raise ValueError("The annotation table is not an excel file")

    # Parse genome file
    genome_dict = hp.parse_fasta(genome)
    combination_arrays = {}

    for trns_file in input_files:
        # Get the name of the current trns file
        trns_file_name = os.path.basename(trns_file)
        trns_file_name = trns_file_name.split(".")[0]

        # Create and fill combination arrays
        combination_arrays[trns_file_name] = hp.make_combination_array(genome_dict)
        th.segemehlTrans2heatmap(trns_file, combination_arrays[trns_file_name])

    # Merge combination arrays
    merged_combination_arrays = ah.combine_arrays(
        combination_arrays, normalise_array=False
    )

    # Check the peak cell for each annotation
    peak_cell_dict = {}
    for index, row in annotation_table.iterrows():
        peak_cel = get_peak_cell_from_annotation_table(merged_combination_arrays, row)
        peak_cell_dict[index] = peak_cel
        # make peak cell dict into a dataframe
        peak_cell_df = pd.DataFrame.from_dict(peak_cell_dict, orient="index")

    # Merge the peak cell data frame with the annotation table
    merged_df = pd.merge(
        annotation_table, peak_cell_df, left_index=True, right_index=True
    )

    # Check if the output directory exists, if not create it
    if not os.path.exists(args["--output"]):
        os.makedirs(args["--output"])

    # Write the output file
    output_file = os.path.basename(genome)
    output_file = output_file.split(".")[0]
    output_file = args["--output"] + "/" + output_file + "_peak_cells.tsv"
    merged_df.to_csv(output_file, sep="\t")


if __name__ == "__main__":
    main()