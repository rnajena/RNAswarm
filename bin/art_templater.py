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
    """
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
            segment_a, start_a, end_a, segment_b, start_b, end_b = line.strip().split(
                ","
            )
            interaction_dict[interaction_idx] = {
                int(segment_a): (int(start_a), int(end_a)),
                int(segment_b): (int(start_b), int(end_b)),
            }
            interaction_idx += 1
    return interaction_dict


def make_interaction_fasta(fasta_filepath, interaction_dict):
    """"
    Make fasta file with interactions.
    
    Args:
        fasta_filepath (str): Path to fasta file.
        interaction_dict (dict): Dictionary of interactions.
        
    Returns:
        str: Fasta file with interactions.
    """
    # This whole function assumes ordered dictionaries, that basically means
    # that it only works on Python 3.6+, which shouldn't be a problem for our
    # purposes.
    genome = [segment for segment in helper.parse_genome(fasta_filepath).values()]
    interaction_fasta = ""
    for idx, interactions in interaction_dict.items():
        first_subkey = True
        interaction_sequence = ""
        for segment, position in interactions.items():
            if first_subkey:
                interaction_fasta = f"{interaction_fasta}>segment_{segment}_start_{position[0]}_end_{position[1]}"
                interaction_sequence = (
                    f"{interaction_sequence}{genome[segment][position[0]:position[1]]}"
                )
                first_subkey = False
            elif not first_subkey:
                interaction_fasta = f"{interaction_fasta}segment_{segment}_start_{position[0]}_end_{position[1]}\n"
                interaction_sequence = (
                    f"{interaction_sequence}{genome[segment][position[0]:position[1]]}"
                )
                interaction_fasta = f"{interaction_fasta}{interaction_sequence}\n"
            else:
                raise ValueError("Caught an exception")
    return interaction_fasta


def main():
    arguments = docopt(__doc__)
    interactions = parse_interactions(arguments["--interactions"])
    interaction_fasta = make_interaction_fasta(arguments["--fasta"], interactions)
    print(interaction_fasta)


if __name__ == "__main__":
    main()


