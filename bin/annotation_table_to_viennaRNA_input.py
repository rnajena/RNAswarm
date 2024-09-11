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
    # add checks to make sure things are within range
    sequences = []
    extension_start01 = extension
    extension_start02 = extension
    extension_end01 = extension
    extension_end02 = extension
    if annotation["start01"] - extension <= 0:
        extension_start01 = 0
    if annotation["start02"] - extension <= 0:
        extension_start02 = 0
    if annotation["end01"] + extension >= len(genome_dict[annotation["segment01"]]):
        extension_end01 = 0
    if annotation["end02"] + extension >= len(genome_dict[annotation["segment02"]]):
        extension_end02 = 0
    if extension:
        sequences = [
            genome_dict[annotation["segment01"]][
                annotation["start01"]
                - extension_start01 : annotation["end01"]
                + extension_end01
            ],
            genome_dict[annotation["segment02"]][
                annotation["start02"]
                - extension_start02 : annotation["end02"]
                + extension_end02
            ],
        ]
    else:
        sequences = [
            genome_dict[annotation["segment01"]][
                annotation["start01"] : annotation["end01"]
            ],
            genome_dict[annotation["segment02"]][
                annotation["start02"] : annotation["end02"]
            ],
        ]
    return sequences


def create_extended_seqs(
    annotation_table_df, genome_dict, extension_window=[5, 50], extension_step=5
):
    """Creates exteded sequences for structure prediction.

    Args:
        annotation_table_df ()
        genome_dict ()
        extension_window ()
        extension_step ()

    Returns:
        dict ()
    """
    sequences_extended = {}
    # iterate over each line of the annotation_table_df
    for i in range(
        extension_window[0], extension_window[1] + extension_step, extension_step
    ):
        sequences_extended[i] = {}
        for index, row in annotation_table_df.iterrows():
            # extract the sequences of the annotation from the genome file
            interaction = annotation_to_sequence(row, genome_dict, extension=i)
            sequences_extended[i][index] = interaction
    return sequences_extended


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
    parser = argparse.ArgumentParser(description="annotation table to viennaRNA input.")
    parser.add_argument("-g", "--genome_file", help="The genome file.")
    parser.add_argument("-a", "--annotation_table", help="The annotation table.")
    parser.add_argument(
        "-e",
        "--extension_window",
        default=[5, 50],
        help="Interval to extend the annotation (inclusive).",
    )
    parser.add_argument(
        "-s", "--extension_step", default=5, help="Step to extend the annotation"
    )
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

    sequences_extended = create_extended_seqs(
        annotation_table_df, genome_dict, extension_window=[5, 50], extension_step=5
    )

    # write the sequences to a fasta file
    output_file = Path(args.output_folder) / Path(
        f"{Path(args.genome_file).stem}_annotations.fasta"
    )

    write_viennaRNA_input(sequences, output_file)
    for extension, sequences_extended in sequences_extended.items():
        output_file_extended = Path(args.output_folder) / Path(
            f"{Path(args.genome_file).stem}_annotations_extended_{extension}.fasta"
        )
        write_viennaRNA_input(sequences_extended, output_file_extended)


if __name__ == "__main__":
    main()
