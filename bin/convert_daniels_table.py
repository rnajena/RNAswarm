#!/usr/bin/env python3

"""convert_daniels_table.py

This script takes Daniel's table and converts it to a table with the following columns: 
id,segment01,start01,end01,segment02,start02,end02

Usage:
    parse_interactions.py <input_file> -t <interaction_type> -o <output_file>

Options:
    -h --help                                 Show this screen.
    <input_file>                              The input files to process, has to be a 
                                              tsv file with the expected columns.
    -t --interaction_type=<interaction_type>  The interaction type to filter for, it can be set to original, cut_structure or peak_structure
    -o --output=<output_file>                 The output filepath.

"""

from docopt import docopt
import os
import pandas as pd


def convert_to_annotation_table(daniels_table_pd, interaction_type=None):
    """
    Convert Daniel's table (which has to be a pandas dataframe) to a table with the following columns: 
    id,segment01,start01,end01,segment02,start02,end02

    Parameters
    ----------
    daniels_table_pd : pandas dataframe
        The Daniel's table as a pandas dataframe
    interaction_type : str
        The interaction type to filter for, it can be set to original, cut_structure or peak_structure

    Returns
    -------
    annotation_table_pd : pandas dataframe
        The annotation table as a pandas dataframe
    """
    # Create a new dataframe with the columns id,segment01,start01,end01,segment02,start02,end02
    annotation_table_pd = pd.DataFrame(
        columns=["id", "segment01", "start01", "end01", "segment02", "start02", "end02"]
    )

    # Add the columns to the dataframe
    if interaction_type == "original":
        annotation_table_pd["id"] = daniels_table_pd["number"]
        annotation_table_pd["segment01"] = daniels_table_pd["aSeq"]
        annotation_table_pd["start01"] = daniels_table_pd["ai"] - 1
        annotation_table_pd["end01"] = daniels_table_pd["aj"]
        annotation_table_pd["segment02"] = daniels_table_pd["bSeq"]
        annotation_table_pd["start02"] = daniels_table_pd["bi"] - 1
        annotation_table_pd["end02"] = daniels_table_pd["bj"]
    elif interaction_type == "cut_structure":
        annotation_table_pd["id"] = daniels_table_pd["number"]
        annotation_table_pd["segment01"] = daniels_table_pd["aSeq"]
        annotation_table_pd["start01"] = daniels_table_pd["cai"] - 1
        annotation_table_pd["end01"] = daniels_table_pd["caj"]
        annotation_table_pd["segment02"] = daniels_table_pd["bSeq"]
        annotation_table_pd["start02"] = daniels_table_pd["cbi"] - 1
        annotation_table_pd["end02"] = daniels_table_pd["cbj"]
    elif interaction_type == "peak_structure":
        annotation_table_pd["id"] = daniels_table_pd["number"]
        annotation_table_pd["segment01"] = daniels_table_pd["aSeq"]
        annotation_table_pd["start01"] = daniels_table_pd["pai"] - 1
        annotation_table_pd["end01"] = daniels_table_pd["paj"]
        annotation_table_pd["segment02"] = daniels_table_pd["bSeq"]
        annotation_table_pd["start02"] = daniels_table_pd["pbi"] - 1
        annotation_table_pd["end02"] = daniels_table_pd["pbj"]
    else:
        raise ValueError(
            "The interaction type has to be one of the following: original, cut_structure or peak_structure"
        )

    return annotation_table_pd


def main():
    # Parse the arguments
    args = docopt(__doc__)
    input_file = args["<input_file>"]
    output_file = args["--output"]
    interaction_type = args["--interaction_type"]

    # Read the input file
    daniels_table_pd = pd.read_csv(input_file, sep="\t", header=0)

    # Convert the table
    annotation_table_pd = convert_to_annotation_table(
        daniels_table_pd, interaction_type=interaction_type
    )

    # check if the output directory exists, if not create it
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Write the table to file
    annotation_table_pd.to_csv(output_file, sep=",", index=False)


if __name__ == "__main__":
    main()
