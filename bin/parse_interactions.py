#!/usr/bin/env python3

"""parse_interactions.py

This script takes an annotation table, a csv file with the following columns:
    - id,
    - segment01,
    - start01,
    - end01,
    - segment02,
    - start02,
    - end02
and a genome file in fasta format and returns a fasta file with the following format:
    >id_segment01
    sequence01
    >id_segment02
    sequence02

Usage:
    parse_interactions.py -a <annotation_table> -g <genome> -o <output_file> [--complement --peaks]

Options:
    -h --help                                 Show this screen.
    -a --annotation_table=<annotation_table>  The annotation table filepath.
    -g --genome=<genome>                      The genome filepath.
    -o --output=<output_file>                 The output filepath.
    --complement                              Return the complement the sequences.
"""

from docopt import docopt
import os
import pandas as pd


def parse_interactions(annotation_table, genome, output_file, complement=False, peaks=False):
    """
    Parse interactions from an annotation table and a genome file.

    Parameters
    ----------
    annotation_table : str
        Path to annotation table.
    genome : str
        Path to genome file.
    output_file : str
        Path to output file.

    Returns
    -------
    None
    """
    annotation_table_df = None
    # Read annotation table which is either an excel, tsv, or csv file
    if annotation_table.endswith(".XLSX"):
        # Ignore any unnamed columns
        annotation_table_df = pd.read_excel(annotation_table, header=0, usecols=lambda x: "Unnamed" not in x)
    elif annotation_table.endswith(".csv"):
        annotation_table_df = pd.read_csv(annotation_table, sep=",", header=0)
    elif annotation_table.endswith(".tsv"):
        annotation_table_df = pd.read_csv(annotation_table, sep="\t", header=0)
    else:
        raise ValueError("Annotation table must be a csv or an excel file.")

    # Read genome
    genome = read_genome(genome)

    # Parse interactions
    with open(output_file, "w") as f:
        for index, row in annotation_table_df.iterrows():
            # Get id
            id = row["id"]

            # Get segment01
            segment01 = row["segment01"]
            start01 = row["start01"]
            end01 = row["end01"]
            if peaks:
                start01 = row["segment01_peak"] - 19
                end01 = row["segment01_peak"] + 20
            sequence01 = genome[segment01][start01:end01]
            if complement:
                sequence01 = sequence01.translate(str.maketrans("ATGC", "TACG"))

            # Get segment02
            segment02 = row["segment02"]
            start02 = row["start02"]
            end02 = row["end02"]
            if peaks:
                start02 = row["segment02_peak"] - 19
                end02 = row["segment02_peak"] + 20
            sequence02 = genome[segment02][start02:end02]
            if complement:
                sequence02 = sequence02.translate(str.maketrans("ATGC", "TACG"))

            is_complement = ""
            if complement:
                is_complement = "_complement"

            # Write to file
            f.write(
                f">{id}_{segment01}_{start01}_{end01}_{is_complement}\n{sequence01}\n>{id}_{segment02}_{start02}_{end02}_{is_complement}\n{sequence02}\n"
            )


def read_genome(genome):
    """
    Read genome from a fasta file.

    Parameters
    ----------
    genome : str
        Path to genome file.

    Returns
    -------
    genome : dict
        Dictionary of sequences.
    """
    genome_dict = {}
    with open(genome, "r") as f:
        for line in f:
            if line.startswith(">"):
                segment = line.strip().split()[0][1:]
                genome_dict[segment] = ""
            else:
                genome_dict[segment] += line.strip()

    return genome_dict


def main():
    # Get arguments
    arguments = docopt(__doc__)

    # Get paths
    annotation_table = os.path.abspath(arguments["--annotation_table"])
    genome = os.path.abspath(arguments["--genome"])
    output_file = os.path.abspath(arguments["--output"])

    # Parse interactions
    parse_interactions(annotation_table, genome, output_file, complement=arguments["--complement"], peaks=arguments["--peaks"])


if __name__ == "__main__":
    main()
