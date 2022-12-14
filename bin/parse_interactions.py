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
    parse_interactions.py -a <annotation_table> -g <genome> -o <output_file>

Options:
    -h --help                                 Show this screen.
    -a --annotation_table=<annotation_table>  The annotation table filepath.
    -g --genome=<genome>                      The genome filepath.
    -o --output=<output_file>                 The output filepath.
"""

from docopt import docopt
import os
import pandas as pd


def parse_interactions(annotation_table, genome, output_file):
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
    # Read annotation table
    annotation_table = pd.read_csv(annotation_table, sep=",", header=0)

    # Read genome
    genome = read_genome(genome)

    # Parse interactions
    with open(output_file, "w") as f:
        for index, row in annotation_table.iterrows():
            # Get id
            id = row["id"]

            # Get segment01
            segment01 = row["segment01"]
            start01 = row["start01"]
            end01 = row["end01"]
            sequence01 = genome[segment01][start01:end01]

            # Get segment02
            segment02 = row["segment02"]
            start02 = row["start02"]
            end02 = row["end02"]
            sequence02 = genome[segment02][start02:end02]

            # Write to file
            f.write(
                f">{id}_{segment01}\n{sequence01}\n>{id}_{segment02}\n{sequence02}\n"
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
    parse_interactions(annotation_table, genome, output_file)


if __name__ == "__main__":
    main()
