#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""normalize_arrays.py

    Normalize multiple numpy arrays.

Usage:
    normalize_arrays.py -d <sample_arrays> -g <genome> -o <sample_arrays_normalized> [--max_value=<max_value> --mode=<mode> --do_not_round]
    normalize_arrays.py -h | --help

Options:
    -h --help       Show this screen.
    -d <sample_arrays>  A repository of numpy arrays
    -g <genome_file>    Genome file
    -o <output>     Output file name
    --max_value=<max_value>  Maximum value to normalise to [default: 200000]
    --mode=<mode>   Mode to normalise the array in. Options are "peak_height" and "number_of_data_points" [default: number_of_data_points]
    --round        Whether to round the values to integers [default: True]
"""

from docopt import docopt
import helper as hp
import array_handler as ah

def main():
    args = docopt(__doc__)
    sample_arrays = args["-d"]
    genome_file = args["-g"]
    output_file = args["-o"]
    max_value = int(args["--max_value"])
    mode = args["--mode"]
    do_not_round = args["--do_not_round"]
    if do_not_round in ["True", "true", "1"]:
        round = False
    else:
        round = True

    genome_dict = hp.parse_fasta(genome_file)
    combination_arrays = hp.make_combination_array(genome_dict, intra_only=False)
    combination_arrays_normalized = combination_arrays.copy()
    ah.import_combination_arrays(combination_arrays, sample_arrays)
    for combination, array in combination_arrays.items():
        combination_arrays_normalized[combination] = ah.normalize_array(array, max_value=max_value, mode=mode, round=round)
    ah.save_combination_arrays(combination_arrays_normalized, output_file)

if __name__ == "__main__":
    main()