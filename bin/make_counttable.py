#!/usr/bin/env python3

"""make_counttable.py

Takes an arbitrary number of trns files, finds and merges interactions, outputing an annotation table.

Usage:
  make_counttable.py <input_file> <input_file>... -a <annotation_table> -o <output_file> 

Options:
    -h --help                                 Show this screen.
    <input_file>                              The input files to process.
    -a --annotation_table=<annotation_table>  The annotation table filepath.
    -o --output=<output_file>                 The output filepath.
  

"""

from docopt import docopt
import pandas as pd
import trns_handler as th


def fill_count_table(interaction, count_table, annotation_table_dict, trns_file):
    """
    Fills the count table with the number of interactions between the given
    regions.

    Parameters
    ----------
    interaction : list
        A list containing the two regions of the interaction and start and end positions.

    count_table : dict
        A dictionary containing the count table.

    annotation_table_dict : dict
        A dictionary containing the annotation table.

    Returns
    -------
    None
    """
    for row in annotation_table_dict.values():
        interaction_name = f"{row['id']}_{row['segment01']}_{row['start01']}_{row['end01']}_{row['segment02']}_{row['start02']}_{row['end02']}"
        if row["segment01"] == interaction[0] and row["segment02"] == interaction[3]:
            if (
                row["start01"] <= interaction[1] <= row["end01"]
                or row["start01"] <= interaction[2] <= row["end01"]
            ) and (
                row["start02"] <= interaction[4] <= row["end02"]
                or row["start02"] <= interaction[5] <= row["end02"]
            ):
                count_table[trns_file][interaction_name] += 1
            else:
                continue
        elif row["segment02"] == interaction[0] and row["segment01"] == interaction[3]:
            if (
                row["start02"] <= interaction[1] <= row["end02"]
                or row["start02"] <= interaction[2] <= row["end02"]
            ) and (
                row["start01"] <= interaction[4] <= row["end01"]
                or row["start01"] <= interaction[5] <= row["end01"]
            ):
                count_table[trns_file][interaction_name] += 1
            else:
                continue


def make_count_table(annotation_table, trns_files):
    """
    Creates a count table from a given annotation table and segemehl trns file.

    Parameters
    ----------
    annotation_table : pandas.DataFrame
        The annotation table containing the regions to count the interactions for.

    trns_files : list
        A list of segemehl trns files.

    Returns
    -------
    count_table : dict
        A dictionary containing the count table.
    """
    annotation_table_dict = annotation_table.to_dict(orient="index")
    count_table = {}
    for trns_file in trns_files:
        count_table[trns_file] = {}
        for row in annotation_table_dict.values():
            count_table[trns_file][
                f"{row['id']}_{row['segment01']}_{row['start01']}_{row['end01']}_{row['segment02']}_{row['start02']}_{row['end02']}"
            ] = 0
    for trns_file in trns_files:
        with open(trns_file) as input_stream:
            for line in input_stream:
                line = line.strip().split()
                firstRead = line[0].split(",")
                secondRead = line[1].split(",")
                currentRow = th.__extract_start_stop_segemehl(
                    firstRead
                ) + th.__extract_start_stop_segemehl(secondRead)
                interaction = th.__check_interaction(currentRow)
                fill_count_table(
                    interaction, count_table, annotation_table_dict, trns_file
                )
    return count_table


def main():
    args = docopt(__doc__)
    # Read annotation table
    annotation_table = pd.read_csv(args["--annotation_table"], sep=",", header=0)
    # Read trns files
    trns_files = args["<input_file>"]
    # Create count table
    count_table = make_count_table(annotation_table, trns_files)
    # Transform count table to pandas.DataFrame
    count_table = pd.DataFrame(count_table)
    # Write count table to file
    count_table.to_csv(args["--output"], sep=",", header=True, index=True)


if __name__ == "__main__":
    main()
