#!/usr/bin/env python3

"""merge_annotation_tables.py

Takes an arbitrary number of annotation tables and merge the lines into one csv

Usage:
    merge_annotation_tables.py <annotation_table>... -o <output_file>

Options:
    -h --help                                 Show this screen.
    <annotation_table>                        The annotation tables to merge.
    -o --output=<output_file>
"""

from docopt import docopt
import pandas as pd

def main():
    args = docopt(__doc__)
    annotation_tables = args["<annotation_table>"]
    output_file = args["--output"]

    merged_annotation_table = pd.DataFrame()

    for annotation_table in annotation_tables:
        annotation_table_df = pd.read_csv(annotation_table, header=None)
        merged_annotation_table = pd.concat([merged_annotation_table, annotation_table_df], ignore_index=True)

    merged_annotation_table.to_csv(output_file, sep="\t", mode="a", header=False)

if __name__ == "__main__":
    main()