#!/usr/bin/env python3

"""read_reporter.py: A script to generate a report comparing mappings and kraken2 results.

Usage:
  read_reporter.py -m <mappings> -k <kraken>

Options:
  -h --help                               Show this screen.
  -m --mappings=<mappings>                bam file with mappings.
  -k --kraken=<kraken>                    kraken output file.
"""

from docopt import docopt
import pysam
import helper as hp

def parse_kraken_output(kraken_file):
    """
    Parse the kraken output file.
    Args:
        kraken_file (str): Path to kraken output file.

    Returns:
        dict: Dictionary of kraken taxonomy, with read id as key and taxonomy as value.
    """
    kraken_dict = {}
    with open(kraken_file, "r") as kraken:
        for line in kraken:
            if line.startswith("C"):
                line = line.split()
                kraken_dict[line[1]] = line[2]
    return kraken_dict

def check_mappings(mappings_file, kraken_dict):
    """
    Check the each mapped read in the mappings file to which taxonomy it belongs.
    Args:
        mappings_file (str): Path to mappings file.
        kraken_dict (dict): Dictionary of kraken taxonomy, with read id as key and taxonomy as value.

    Returns:
        dict: Dictionary of kraken taxonomy, with read id as key and taxonomy as value.
    """
    with pysam.AlignmentFile(mappings_file, "rb") as bam:
        summary_dict = {'unmapped': {}, 'mapped': {}}
        for read in bam:
            if read.is_unmapped:
                if read.query_name in kraken_dict:
                    if kraken_dict[read.query_name] in summary_dict['unmapped']:
                        summary_dict['unmapped'][kraken_dict[read.query_name]] += 1
                    else:
                        summary_dict['unmapped'][kraken_dict[read.query_name]] = 1
                elif read.query_name not in kraken_dict:
                    if 'unclassified' in summary_dict['unmapped']:
                        summary_dict['unmapped']['unclassified'] += 1
                    else:
                        summary_dict['unmapped']['unclassified'] = 1
                else:
                    raise ValueError("Caught an exception")
            elif read.is_mapped:
                if read.query_name in kraken_dict:
                    if kraken_dict[read.query_name] in summary_dict['mapped']:
                        summary_dict['mapped'][kraken_dict[read.query_name]] += 1
                    else:
                        summary_dict['mapped'][kraken_dict[read.query_name]] = 1
                elif read.query_name not in kraken_dict:
                    if 'unclassified' in summary_dict['mapped']:
                        summary_dict['mapped']['unclassified'] += 1
                    else:
                        summary_dict['mapped']['unclassified'] = 1
                else:
                    raise ValueError("Caught an exception")
            else:
                raise ValueError("Caught an exception")
    return summary_dict

def main():
    """
    Main function.
    """
    arguments = docopt(__doc__)
    kraken_dict = parse_kraken_output(arguments["--kraken"])
    summary_dict = check_mappings(arguments["--mappings"], kraken_dict)
    print("Unmapped reads:")
    for taxonomy in summary_dict['unmapped']:
        print(taxonomy + ": " + str(summary_dict['unmapped'][taxonomy]))
    print("Mapped reads:")
    for taxonomy in summary_dict['mapped']:
        print(taxonomy + ": " + str(summary_dict['mapped'][taxonomy]))

if __name__ == '__main__':
    main()