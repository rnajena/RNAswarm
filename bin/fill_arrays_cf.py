#!/usr/bin/env python3

"""fill_arrays.py

Usage:
    fill_arrays.py <cf_file> -g <genome> [--intra_only] -o <output_folder>

Options:
    -h --help                    Show this screen.
    <cf_file>                    Path to ChimericFragments files
    -g --genome=<genome>         The genome filepath.
    -o --output=<output_folder>  The output folder.
    --intra_only                 Only plot intra-segment interactions.

"""

from docopt import docopt
import os
import helper as hp
import cf_handler as ch
import array_handler as ah


def main():
    args = docopt(__doc__)
    cf_file = args["<cf_file>"]
    genome_file = args["--genome"]
    output_folder = args["--output"]
    # Check if --intra_only is given
    intra_only = False
    if args["--intra_only"]:
        intra_only = True

    # Process input files
    genome_dict = hp.parse_fasta(genome_file)
    
    # Create and fill combination dicts
    combination_dict = hp.make_combination_array(genome_dict, intra_only=intra_only)
    ch.chimericFragments2heatmap(cf_file, combination_dict, intra_only=intra_only)

    # Save combination arrays
    ah.save_combination_arrays(combination_dict, output_folder)


if __name__ == "__main__":
    main()