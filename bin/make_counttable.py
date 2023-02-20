#!/usr/bin/env python3

"""make_counttable.py

Takes an arbitrary number of trns files and an annotation table and creates a
count table containing the number of interactions between the regions in the
annotation table.

Usage:
  make_counttable.py <input_file> <input_file>... -a <annotation_table> -o <output_file> [--use_peaks]
  make_counttable.py <input_file> -a <annotation_table> -o <output_file> [--use_peaks]

Options:
  -h --help                                 Show this screen.
  <input_file>                              The input files to process, has to be a trns file generated by segemehl.
  -a --annotation_table=<annotation_table>  The annotation table filepath.
  -o --output=<output_file>
  --use_peaks                               Use the peak regions instead of the full regions.
"""

from docopt import docopt
import pandas as pd
import trns_handler as th


def fill_count_table(interaction, count_table, annotation_table_dict, trns_file, use_peaks=False, window_size=20):
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

    trns_file : str
        The trns file the interaction was found in.

    use_peaks : bool, optional
        Use the peak regions instead of the full regions, by default False

    window_size : int, optional
        The window size to use for the peak regions, by default 20

    Returns
    -------
    None
    """
    # Loop through the annotation table dictionary
    for row in annotation_table_dict.values():
        # Get the interaction id
        interaction_id = int(row['id'])
        # Initialize start and end positions for segments
        start01, end01, start02, end02 = 0, 0, 0, 0
        # Calculate start and end positions for segments
        if use_peaks:
            start01 = row['segment01_peak'] - window_size
            end01 = row['segment01_peak'] + window_size
            start02 = row['segment02_peak'] - window_size
            end02 = row['segment02_peak'] + window_size
        else:
            start01, end01, start02, end02 = row['start01'], row['end01'], row['start02'], row['end02']
        if row["segment01"] == interaction[0] and row["segment02"] == interaction[3]:
            if (
                start01 <= interaction[1] <= end01 or
                start01 <= interaction[2] <= end01 or
                interaction[1] <= start01 <= interaction[2] or
                interaction[1] <= end01 <= interaction[2]
            ) and (
                start02 <= interaction[4] <= end02 or
                start02 <= interaction[5] <= end02 or
                interaction[4] <= start02 <= interaction[5] or
                interaction[4] <= end02 <= interaction[5]
            ):
                count_table[trns_file][interaction_id] += 1
        if row["segment02"] == interaction[0] and row["segment01"] == interaction[3]:
            if (
                start01 <= interaction[4] <= end01 or
                start01 <= interaction[5] <= end01 or
                interaction[4] <= start01 <= interaction[5] or
                interaction[4] <= end01 <= interaction[5]
            ) and (
                start02 <= interaction[1] <= end02 or
                start02 <= interaction[2] <= end02 or
                interaction[1] <= start02 <= interaction[2] or
                interaction[1] <= end02 <= interaction[2]
            ):
                count_table[trns_file][interaction_id] += 1


def make_count_table(annotation_table, trns_files, use_peaks=False):
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
            count_table[trns_file][int(row['id'])] = 0
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
                    interaction, count_table, annotation_table_dict, trns_file, use_peaks=use_peaks
                )
    return count_table


def main():
    args = docopt(__doc__)
    # Read annotation table
    annotation_table = pd.read_csv(args["--annotation_table"], sep="\t", header=0)
    # Read trns files
    trns_files = args["<input_file>"]
    # Create count table
    count_table = make_count_table(annotation_table, trns_files, use_peaks=args["--use_peaks"])
    # Transform count table to pandas.DataFrame
    count_table = pd.DataFrame(count_table)
    # Write count table to file
    count_table.to_csv(args["--output"], sep="\t", header=True, index=True)


if __name__ == "__main__":
    main()