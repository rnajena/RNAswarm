#!/usr/bin/env python3

"""fill_arrays.py

Usage:
    fill_arrays.py <trns_file>... -g <genome> [--intra_only] -o <output_folder>

Options:
    -h --help                    Show this screen.
    <trns_file>                  Path to trns files
    -g --genome=<genome>         The genome filepath.
    -o --output=<output_folder>  The output folder.
    --intra_only                 Only plot intra-segment interactions.

"""

from docopt import docopt
import os
import helper as hp
import trns_handler as th
import array_handler as ah


def main():
    args = docopt(__doc__)
    trns_files = args["<trns_file>"]
    genome_file = args["--genome"]
    output_folder = args["--output"]
    # Check if --intra_only is given
    intra_only = False
    if args["--intra_only"]:
        intra_only = True

    # Process input files
    genome_dict = hp.parse_fasta(genome_file)

    # Process trns files
    combination_dicts = {}

    # Create and fill combination dicts
    for trns_file in trns_files:
        trns_file_name = os.path.basename(trns_file)
        combination_dicts[trns_file_name] = hp.make_combination_array(genome_dict, intra_only=intra_only)
        th.segemehlTrans2heatmap(trns_file, combination_dicts[trns_file_name], intra_only=intra_only)

    # Save combination arrays
    for trns_file_name, combination_dict in combination_dicts.items():
        ah.save_combination_arrays(combination_dict, output_folder)


if __name__ == "__main__":
    main()