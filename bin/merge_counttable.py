#!/usr/bin/env python3

"""merge_counttable.py

Takes an arbitrary number of count tables and merge the columns into one csv

Usage:
  merge_counttable.py <count_table>... -o <output_file>

Options:
  -h --help                                 Show this screen.
  <count_table>                             The count tables to merge.
  -o --output=<output_file>
"""

from docopt import docopt
import pandas as pd


def check_first_column(count_tables):
    """
    Checks if the first column of all count tables are the same.

    Parameters
    ----------
    count_tables : list
        A list containing the count tables.

    Returns
    -------
    bool
        True if the first column of all count tables are the same, False otherwise.
    """
    # Get the first column of the first count table
    first_column = pd.read_csv(count_tables[0], sep="\t", usecols=[0]).iloc[:, 0]

    # Iterate over all count tables
    for count_table in count_tables[1:]:
        # Get the first column of the count table
        first_column_tmp = pd.read_csv(
            count_table, sep="\t", usecols=[0]
        ).iloc[:, 0]

        # Check if the first column of the count table is the same as the first column of the first count table
        if not first_column.equals(first_column_tmp):
            return False

    return True


def main():
    args = docopt(__doc__)
    count_tables = args["<count_table>"]
    output_file = args["--output"]

    # Check if the first column of all count tables are the same
    if not check_first_column(count_tables):
        raise ValueError("The first column of all count tables has to be the same.")
    
    # Read the first count table
    count_table = pd.read_csv(count_tables[0], sep="\t", index_col=0)

    # Iterate over all count tables
    for count_table_tmp in count_tables[1:]:
        # Read the count table
        count_table_tmp = pd.read_csv(count_table_tmp, sep="\t", index_col=0)

        # Merge the count table
        count_table = count_table.join(count_table_tmp)

    # Write the count table
    count_table.to_csv(output_file, sep="\t")

if __name__ == "__main__":
    main()
