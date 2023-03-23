#!/usr/bin/env python3

"""get_library_size.py

Usage:
    get_library_size.py -f <bam_file> [--output=<output>]
    get_library_size.py -f <bam_file>... [--output=<output>]
    get_library_size.py -F <bam_file_folder> [--output=<output>]
    get_library_size.py -F <bam_file_folder>... [--output=<output>]
    get_library_size.py -h | --help

Options:
    -h --help           Show this screen.
    -f --file           Bam file
    -F --folder         Bam file folder
    --output=<output>   Output file [default: library_size.txt]

"""

from docopt import docopt
import pysam
import sys
import os
import multiprocessing as mp

def get_library_size(bam_file):
    """Get library size from bam file

    Args:
        bam_file (str): Bam file

    Returns:
        int: Library size
    """
    library_size = pysam.view('-c', '-F', '4', bam_file)
    return library_size

def main():
    """Main function
    """
    args = docopt(__doc__)
    if args['--folder']:
        bam_file_list = []
        for bam_file_folder in args['<bam_file_folder>']:
            bam_file_list.extend([bam_file_folder + '/' + bam_file for bam_file in os.listdir(bam_file_folder) if bam_file.endswith('.bam')])
    else:
        bam_file_list = [args['--file']]
    # Check if --output was given
    if args['--output']:
        output_file = args['--output']
    else:
        # Print to stdout
        output_file = sys.stdout
    # Open output file
    # Check it output folder exists
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))
    # In parallel process each bam file and write to output file
    with mp.Pool() as pool, open(output_file, 'w') as output:
        for bam_file, library_size in zip(bam_file_list, pool.map(get_library_size, bam_file_list)):
            print(f"processing {bam_file}...")
            output.write(f"{bam_file}\t{library_size}")


if __name__ == '__main__':
    main()