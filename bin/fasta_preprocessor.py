#!/usr/bin/env python3

"""fasta preprocessor.

Usage:
    fasta_preprocessor.py -c -i <input_multifasta> -o <output_fasta>

Options:
    -h --help                               Show this screen.
    -i --input=<input_fasta>                input multifasta file
    -c --concatenate                        Concatenate multifasta file
    -o --output=<output_fasta>              output fasta file
"""

from fileinput import filename
from posixpath import basename
from docopt import docopt
import helper


def check_multifasta(fasta_file):
    """
    check if the fasta file is a multifasta file
    """
    i = 0
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                i += 1
    if i > 1:
        return True
    elif i == 1:
        return False
    else:
        raise Exception("File is not a fasta file")


def concatenate_multifasta(fasta_file):
    fasta_dict = helper.parse_fasta(fasta_file)
    concatenated_seq = ""
    concatenated_multifasta = ""
    fasta_filename = filename(fasta_file)
    metadata_dict = {}
    for id, seq in fasta_dict.items():
        metadata_dict[id] = {
            "length": len(seq),
            "start": len(concatenated_seq),
            "end": len(concatenated_seq) + len(seq) - 1,
        }
        concatenated_seq += seq
        concatenated_multifasta += f">{fasta_filename}\n{seq}\n"
    return concatenated_multifasta, metadata_dict


def main():
    arguments = docopt(__doc__)
    concatenated_multifasta, metadata_dict = concatenate_multifasta(
        arguments["--input"]
    )
    with open(arguments["--output"], "w") as output:
        output.write(concatenated_multifasta)
    with open(f"{arguments['--output'].basename}.csv", "w") as output:
        output.write(f"id,length,start,end\n")
        for id, metadata in metadata_dict.items():
            output.write(
                f"{id},{metadata['length']},{metadata['start']},{metadata['end']}\n"
            )


if __name__ == "__main__":
    main()
