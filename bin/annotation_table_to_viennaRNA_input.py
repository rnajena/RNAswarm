#!/usr/bin/env python3

import argparse
import deduplicate_annotations as da
import helper as hp
import pandas as pd
from pathlib import Path


def annotation_to_sequence(annotation, genome_dict, extension=0):
    """Extracts the sequences of the annotation from the genome file.

    Args:
        annotation (pandas.Series): The annotation.
        genomes_dict (dict): A dictionary containing the genome.
        extension (int): The extension to add to the start and end of the annotation.
    Returns:
        sequences (list): The sequences of the annotation.
    """
    sequences = [
        genome_dict[annotation["segment01"]][
            annotation["start01"] + extension : annotation["end01"] + extension
        ],
        genome_dict[annotation["segment02"]][
            annotation["start02"] + extension : annotation["end02"] + extension
        ],
    ]
    return sequences


def write_viennaRNA_input(sequences, output_file):
    """Writes the sequences to an input file for viennaRNA's RNAduplex

    Args:
        sequences (dict): The sequences to write.
        output_file (str): The output file.
    """
    with open(output_file, "w") as file:
        for index, sequence in sequences.items():
            file.write(f">{index}\n")
            file.write(f"{sequence[0]}\n")
            file.write(f"{sequence[1]}\n")


def main():
    parser = argparse.ArgumentParser(description="annotation table to fasta.")
    parser.add_argument("-g", "--genome_file", help="The genome file.")
    parser.add_argument("-a", "--annotation_table", help="The annotation table.")
    parser.add_argument("-o", "--output_folder", help="The output folder.")
    args = parser.parse_args()

    annotation_table_df = da.parse_annotation_table(args.annotation_table)
    genome_dict = hp.parse_fasta(args.genome_file)


    sequences = {}
    # iterate over each line of the annotation_table_df
    for index, row in annotation_table_df.iterrows():
        # extract the sequences of the annotation from the genome file
        interaction = annotation_to_sequence(row, genome_dict)
        sequences[index] = interaction


    sequences_extended = {}
    # iterate over each line of the annotation_table_df
    for index, row in annotation_table_df.iterrows():
        # extract the sequences of the annotation from the genome file
        interaction = annotation_to_sequence(row, genome_dict, extension=25)
        sequences_extended[index] = interaction


    # write the sequences to a fasta file
    output_file = Path(args.output_folder) / Path(f"{Path(args.genome_file).stem}_annotations.fasta")
    output_file_extended = Path(args.output_folder) / Path(f"{Path(args.genome_file).stem}_annotations_extended.fasta")
    write_fasta(sequences, output_file)
    write_fasta(sequences_extended, output_file_extended)


if __name__ == "__main__":
    main()
