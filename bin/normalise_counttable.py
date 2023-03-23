#!/usr/bin/env python3

"""normalise_counttable.py

Normalises the count table by dividing the counts by the total number mapped
reads.

Usage:
    normalise_counttable.py <samples_table> -c <count_table> --output=<output_repository>
    normalise_counttable.py -h | --help

Options:
    -h --help                               Show this screen.
    <samples_table>                         The samples table filepath.
    -c <count_table>                        The count table filepath.
    --output=<output_repository>            The output repository.

"""

from docopt import docopt
import pandas as pd
import os
import sys

def parse_samples_table(samples_table):
    """
    Parses the samples table and returns a dictionary containing the sample
    names and the corresponding count table filepaths.

    Parameters
    ----------
    samples_table : str
        The samples table filepath.

    Returns
    -------
    dict
        A dictionary containing the sample names and the corresponding count
        table filepaths.
    """
    samples_dict = {}
    with open(samples_table, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            sample_name = line[0]
            library_size = line[1]
            samples_dict[sample_name] = library_size
    return samples_dict

def parse_count_table(count_table):
    """
    Parses the count table and returns a pandas.DataFrame containing the count
    table.

    Parameters
    ----------
    count_table : str
        The count table filepath.

    Returns
    -------
    pandas.DataFrame
        A pandas.DataFrame containing the count table.
    """
    count_table = pd.read_csv(count_table, sep="\t", index_col=0)
    return count_table

def normalise_count_table(count_table, samples_dict, mode="RPM"):
    """
    Normalises the count table by dividing the counts by the total number mapped
    reads.

    Parameters
    ----------
    count_table : pandas.DataFrame
        The count table.

    total_reads : int
        The total number of mapped reads.

    mode : str, optional
        The normalisation mode, by default "RPM"

    Returns
    -------
    pandas.DataFrame
        The normalised count table.
    """
    # For each column get the corresponding library size and divide the counts
    # by the library size.
    for column in count_table.columns:
        # Check if the sample_name is a substring of the column name.
        # if no key is found in the whole dictionary, raise an error.
        for key, value in samples_dict.items():
            if key in column:
                # Get the library size.
                library_size = value
                # Divide the counts by the library size.
                count_table[column] = count_table[column].div(int(library_size))
                break
        else:
            sys.exit(f"ERROR: Sample {column} not found in samples table.")
    # If the normalisation mode is RPM.
    if mode == "RPM":
        # Multiply the counts by 1,000,000.
        count_table = count_table.mul(1e6)
    return count_table

def write_count_table(count_table, normalised_count_table_df, output_repository):
    """
    Writes the normalised count table to a file.

    Parameters
    ----------
    count_table : str
        The count table filepath.

    normalised_count_table_df : pandas.DataFrame
        The normalised count table.

    output_repository : str
        The output repository.

    Returns
    -------
    None
    """
    # Get the count table filename (without the extension)
    count_table_filename = os.path.basename(count_table).split(".")[0]
    # Create the output filepath.
    output_filepath = os.path.join(output_repository, f"{count_table_filename}_normalised.tsv")
    # Write the normalised count table to a file.
    normalised_count_table_df.to_csv(output_filepath, sep="\t")

def main():
    # Parse the command line arguments.
    args = docopt(__doc__)
    # Get the samples table filepath.
    samples_table = args["<samples_table>"]
    # Get the output repository.
    output_repository = args["--output"]
    # Create the output repository if it does not exist.
    if not os.path.exists(output_repository):
        os.makedirs(output_repository)
    # Parse the samples table.
    samples_dict = parse_samples_table(samples_table)
    # Get the count table filepath.
    count_table = args["-c"]
    # Parse the count table.
    count_table_df = parse_count_table(count_table)
    # Normalise the count table.
    normalised_count_table_df = normalise_count_table(count_table_df, samples_dict)
    # Write the normalised count table to a file.
    write_count_table(count_table, normalised_count_table_df, output_repository)

if __name__ == "__main__":
    main()

