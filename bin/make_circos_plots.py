#!/usr/bin/env python3

"""make_circos_plots.py

Takes the results of DESeq2, the annotation table and the genome files and
creates a circos plot for each comparison.

Usage:
  make_circos_plots.py <input_file> -g <genome_file> -o <output_dir>
  make_circos_plots.py -h | --help

Options:
    -h --help                                 Show this screen.
    <input_file>                              The input files to process, has to be a DESeq2 result file.
    -g --genome=<genome_file>                 The genome file filepath.
    -o --output=<output_dir>                  The output directory.
    
    """

from docopt import docopt
import pandas as pd
import helper as hp
import os
import subprocess
import sys


def make_circos_files(
    DESeq2_results,
    genome_dict,
    output_dir,
    number_of_top_hits=20,
):
    """
    Creates a circos plot for the given comparison.
    Sequences with a highly positive log2FoldChange are represented in red,
    sequences with a highly negative log2FoldChange are represented in blue.

    Parameters
    ----------
    DESeq2_results : pandas.DataFrame
        The DESeq2 results.
    genome_dict : dict
        The genome dictionary.
    number_of_top_hits : int, optional
        The number of top hits to plot, by default 20
    output_dir : str
        The output directory.
    Returns
    -------
    None
    """
    # Get the top hits
    top_hits = DESeq2_results.sort_values(by="log2FoldChange", ascending=False).head(
        number_of_top_hits
    )
    # Get the bottom hits
    bottom_hits = DESeq2_results.sort_values(by="log2FoldChange", ascending=True).head(
        number_of_top_hits
    )
    # Combine the top and bottom hits
    hits = pd.concat([top_hits, bottom_hits])
    # Copy the hits data frame to hits_with_positions
    hits_with_positions = hits.copy()
    # Extract hit position in the genome
    # example: 44_SC35M_NP_25_100_SC35M_M_37_150
    # iterate over the hits
    for index, row in hits.iterrows():
        # split the hit id at the underscore
        hit_name = row["id"].split("_")
        # get the start and end positions in segment 1
        start_1 = int(hit_name[3])
        end_1 = int(hit_name[4])
        # get the start and end positions in segment 2
        start_2 = int(hit_name[7])
        end_2 = int(hit_name[8])
        # get the segment names
        segment_1 = f"{hit_name[1]}_{hit_name[2]}"
        segment_2 = f"{hit_name[5]}_{hit_name[6]}"
        # add the segment names and positions to the hits_with_positions data frame
        hits_with_positions.loc[index, "segment_1"] = segment_1
        hits_with_positions.loc[index, "start_1"] = start_1
        hits_with_positions.loc[index, "end_1"] = end_1
        hits_with_positions.loc[index, "segment_2"] = segment_2
        hits_with_positions.loc[index, "start_2"] = start_2
        hits_with_positions.loc[index, "end_2"] = end_2

    # Create a dict of segments
    segments = {}
    # Iterate over the segments in the genome dictionary
    for segment, sequence in genome_dict.items():
        # Get the length of the segment
        length = len(sequence)
        # Add the segment to the segments dict
        segments[segment] = length

    # Create karyotype file
    with open(f"{output_dir}/karyotype.txt", "w") as karyotype_file:
        # Iterate over the segments
        chr_number = 1
        for segment, length in segments.items():
            # Write the segment to the karyotype file with biological positions
            karyotype_file.write(
                f"chr - {segment} {segment} 1 {length} chr{chr_number}\n"
            )
            chr_number += 1

    # Create the hits.txt file
    with open(f"{output_dir}/hits.txt", "w") as hits_file:
        # Iterate over the hits
        for index, row in hits_with_positions.iterrows():
            # Write the hit to the hits.txt file
            hits_file.write(
                f"{row['id']} {row['segment_1']} {int(row['start_1'])} {int(row['end_1'])}\n"
                f"{row['id']} {row['segment_2']} {int(row['start_2'])} {int(row['end_2'])}\n"
            )

    # Create the circos.conf file
    # hits are to be represented as links in the final circos plot
    with open(f"{output_dir}/circos.conf", "w") as circos_file:
        # Write the circos.conf file
        circos_file.write(
            f"karyotype = {output_dir}/karyotype.txt\n"

            f"<ticks>\n"
            f"chromosomes_display_default = yes\n"
            f"</ticks>\n"

            f"<links>\n"
            f"file = {output_dir}/hits.txt\n"
            f"radius = 0.85r\n"
            f"bezier_radius = 0r\n"
            f"color = black\n"
            f"thickness = 2p\n"
            f"</links>\n"

            f"<image>\n"
            f"<<include etc/image.conf>>\n"
            f"</image>\n"
            
            f"<<include etc/housekeeping.conf>>\n"
        )


def main():
    """
    Main function of the script.

    Returns
    -------
    None
    """
    # Parse command line arguments
    args = docopt(__doc__)
    input_file = args["<input_file>"]
    genome_file = args["--genome"]
    output_dir = args["--output"]

    # Check if output directory exists, create it if not
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read input file with pandas, rename unnamed column to id
    DESeq2_results = pd.read_csv(input_file, sep=",")
    DESeq2_results = DESeq2_results.rename(columns={"Unnamed: 0": "id"})

    # Parse genome segment names and sequences
    genome_dict = hp.parse_fasta(genome_file)

    # Create circos plot files
    make_circos_files(DESeq2_results, genome_dict, output_dir)

    # Create circos plot
    run_circos = subprocess.run(
        ["circos", "-conf", f"{output_dir}/circos.conf"], capture_output=True
    )
    # Check if circos ran successfully
    if run_circos.returncode != 0:
        print("Circos failed to run")
        print(run_circos.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
