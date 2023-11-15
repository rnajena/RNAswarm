#!/usr/bin/env python3

"""trns_parser.py

Usage:
  trns_parser.py -t <trns_file> -f <fastq_file> [-o <output_file>]
  trns_parser.py -h | --help

Options:
  -h --help                         Show this screen.
  -t --trns_file=<trns_file>        trns file to parse.
  -f --fastq_file=<fastq_file>      fastq file to parse.
  -o --output_file=<output_file>    output file to write to.
"""

from docopt import docopt

def load_fastq_to_dict(fastq_file):
    """
    Load the entire fastq file into a dictionary for quick access.

    Args:
        fastq_file (str): Path to fastq file.

    Returns:
        dict: Dictionary with read ids as keys and fastq sequences as values.
    """
    fastq_dict = {}
    with open(fastq_file) as file:
        while True:
            header = file.readline().strip()
            if not header:
                break
            sequence = file.readline().strip()
            _ = file.readline()  # + line
            quality = file.readline().strip()
            read_id = header.split(" ")[0][1:]  # remove '@'
            fastq_dict[read_id] = f"{header}\n{sequence}\n+\n{quality}"
    return fastq_dict


def trns_line_to_read_id(line):
    """
    Convert trns line to read id.

    Args:
        line (str): trns line.

    Returns:
        str: Read id.
    """
    return line.strip().split("\t")[2]


def main():
    arguments = docopt(__doc__)
    trns_file = arguments["--trns_file"]
    fastq_file = arguments["--fastq_file"]
    output_file = arguments["--output_file"]

    fastq_dict = load_fastq_to_dict(fastq_file)
    output_fastq = []

    with open(trns_file) as file:
        for line in file:
            read_id = trns_line_to_read_id(line)
            read = fastq_dict.get(read_id)
            if read:
                output_fastq.append(read)
            else:
                raise ValueError("Read not found")

    if output_file:
        with open(output_file, "w") as output:
            output.write("\n".join(output_fastq))
    else:
        print("\n".join(output_fastq))


if __name__ == "__main__":
    main()
