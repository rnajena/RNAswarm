#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""merge_arrays.py

Merge multiple numpy arrays into one

Usage:
    merge_arrays.py <array_folder>... -g <genome_file> -o <output>
    merge_arrays.py -h | --help

Options:
    -h --help       Show this screen.
    <array_folder>  Folder containing numpy arrays
    -g <genome_file>    Genome file
    -o <output>     Output file name
"""

from docopt import docopt
import helper as hp
import array_handler as ah

def main():
    args = docopt(__doc__)
    array_folders = args["<array_folder>"]
    genome_file = args["<genome_file>"]
    output = args["<output>"]

    # Process input files
    genome_dict = hp.parse_fasta(genome_file)

    # Create and fill combination arrays
    combination_arrays = {}
    for array_folder in array_folders:
        combination_arrays[array_folder] = hp.make_combination_array(genome_dict, intra_only=True)
        ah.import_combination_arrays(combination_arrays[array_folder], array_folder)

    # Merge combination arrays
    merged_combination_arrays = ah.combine_arrays(combination_arrays)

    # Save merged combination arrays
    ah.save_combination_arrays(merged_combination_arrays, output)

if __name__ == "__main__":
    main()
