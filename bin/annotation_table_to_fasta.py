#!/usr/bin/env python3

"""annotation table to fasta.

Usage:
    annotation_table_to_fasta.py -h | --help
    annotation_table_to_fasta.py -a <annotation_table> -g <genome_file>... -o <output_file>

Options:
    -h --help                         Show this screen.
    -a --annotation_table=<annotation_table>     The annotation table.
    -g --genome_file=<genome_file>     The genome file. Can be specified multiple times.
    -o --output_file=<output_file>     The output file.
"""
from docopt import docopt
import deduplicate_annotations as da
from Bio import SeqIO
import pandas as pd


def annotation_to_sequence(annotation, genome_file):
    """Extracts the sequence of the annotation from the genome file.

    Args:
        
    Returns:
        sequence (str): The sequence of the annotation.
    """
    genome = SeqIO.read(genome_file, "fasta")
    sequence = genome.seq[annotation['start']-1:annotation['end']]
    if annotation['strand'] == '-':
        sequence = sequence.reverse_complement()
    return sequence


def main():
    args = docopt(__doc__)
    annotation_table = args['--annotation_table']
    genome_file = args['--genome_file']
    output_file = args['--output_file']
    annotation_table_df = da.parse_annotation_table(annotation_table)


if __name__ == "__main__":
    main()
