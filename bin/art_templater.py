#!/usr/bin/env python3

"""ART templater.

Usage:
  art_templater.py -i <interactions> -f <fasta>

Options:
  -h --help                         Show this screen.
  -i --interactions=<interactions>  Table of interactions to simulate. 
  -f --fasta=<fasta>                Fasta file to use as template.

"""
from docopt import docopt
import numpy as np
import helper


def parse_interactions(interactions):
    """"
    Parse interactions file.

    Args:
        interactions (str): Path to interactions file.

    Returns:
        dict: Dictionary of interactions.
    """
    interaction_dict = {}
    interaction_idx = 0

    with open(interactions) as file:
        for line in file:
            segment_a, start_a, end_a, segment_b, start_b, end_b = line.strip().split(",")
            interaction_dict[interaction_idx] = {
                int(segment_a): (int(start_a), int(end_a)),
                int(segment_b): (int(start_b), int(end_b)),
            }
            interaction_idx += 1
    return interaction_dict


def make_interaction_fasta(fasta_filepath, interaction_dict):
    genome = helper.parse_genome(fasta_filepath)
    interaction_fasta = ""
    for idx, interactions in interaction_dict.items():
        first_subkey = True
        interaction_sequence = ""
        if first_subkey:
            for segment, position in interactions.items():
                interaction_fasta = f"{interaction_fasta}>"
                interaction_sequence = f"{interaction_sequence}{genome[segment][position[0]:position[1]]}"
        elif not first_subkey:
            for segment, position in interactions.items():
                interaction_fasta = f"{interaction_fasta}/n"
                interaction_sequence = f"{interaction_sequence}{genome[segment][position[0]:position[1]]}"
                interaction_fasta = f"{interaction_fasta}{interaction_sequence}/n"
        else:
            raise ValueError("Caught an exception")
    return interaction_fasta


def main():
    arguments = docopt(__doc__)
    interactions = parse_interactions(arguments["--interactions"])
    print(arguments)
    print(interactions)
    print(make_interaction_fasta(arguments["--fasta"], interactions))


if __name__ == "__main__":
    main()
