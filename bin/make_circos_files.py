#!/usr/bin/env python3

"""make_circos_files.py

Takes the results of DESeq2, the annotation table and the genome files and
creates a circos plot for each comparison.

Usage:
  make_circos_files.py <input_file> -a <annotation_table> -g <genome_file> -o <output_dir>
  make_circos_files.py -h | --help

Options:
    -h --help                                 Show this screen.
    <input_file>                              The input files to process, has to be a DESeq2 result file.
    -a --annotation=<annotation_table>        The annotation table.
    -g --genome=<genome_file>                 The genome file filepath.
    -o --output=<output_dir>                  The output directory.
    
    """

from docopt import docopt
import pandas as pd
import helper as hp
import os


def make_circos_files(
    DESeq2_results,
    annotation_table,
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
    # Copy the hits data frame to hits_with_positions
    top_hits_with_positions = top_hits.copy()
    bottom_hits_with_positions = bottom_hits.copy()
    # Create columns for segment_1, start_1, end_1, segment_2, start_2, end_2
    top_hits_with_positions["segment_1"] = ""
    top_hits_with_positions["start_1"] = 0
    top_hits_with_positions["end_1"] = 0
    top_hits_with_positions["segment_2"] = ""
    top_hits_with_positions["start_2"] = 0
    top_hits_with_positions["end_2"] = 0
    bottom_hits_with_positions["segment_1"] = ""
    bottom_hits_with_positions["start_1"] = 0
    bottom_hits_with_positions["end_1"] = 0
    bottom_hits_with_positions["segment_2"] = ""
    bottom_hits_with_positions["start_2"] = 0
    bottom_hits_with_positions["end_2"] = 0
    # Extract hit position in the genome from the annotation table (id,segment01,start01,end01,segment02,start02,end02)
    # and add it to the hits_with_positions data frame
    # iterate over the hits
    for index, row in top_hits_with_positions.iterrows():
        # split the hit id at the underscore
        hit_name = row["id"]
        # get the start and end positions in segment 1
        start_1 = annotation_table.loc[annotation_table["id"] == hit_name, "start01"].iloc[0]
        end_1 = annotation_table.loc[annotation_table["id"] == hit_name, "end01"].iloc[0]
        # get the start and end positions in segment 2
        start_2 = annotation_table.loc[annotation_table["id"] == hit_name, "start02"].iloc[0]
        end_2 = annotation_table.loc[annotation_table["id"] == hit_name, "end02"].iloc[0]
        # get the segment names
        segment_1 = str(annotation_table.loc[annotation_table["id"] == hit_name, "segment01"].iloc[0])
        segment_2 = str(annotation_table.loc[annotation_table["id"] == hit_name, "segment02"].iloc[0])
        # add the segment names and positions to the hits_with_positions data frame
        top_hits_with_positions.loc[index, "segment_1"] = str(segment_1)
        top_hits_with_positions.loc[index, "start_1"] = int(start_1)
        top_hits_with_positions.loc[index, "end_1"] = int(end_1)
        top_hits_with_positions.loc[index, "segment_2"] = str(segment_2)
        top_hits_with_positions.loc[index, "start_2"] = int(start_2)
        top_hits_with_positions.loc[index, "end_2"] = int(end_2)
    for index, row in bottom_hits_with_positions.iterrows():
        # split the hit id at the underscore
        hit_name = row["id"]
        # get the start and end positions in segment 1
        start_1 = annotation_table.loc[annotation_table["id"] == hit_name, "start01"].iloc[0]
        end_1 = annotation_table.loc[annotation_table["id"] == hit_name, "end01"].iloc[0]
        # get the start and end positions in segment 2
        start_2 = annotation_table.loc[annotation_table["id"] == hit_name, "start02"].iloc[0]
        end_2 = annotation_table.loc[annotation_table["id"] == hit_name, "end02"].iloc[0]
        # get the segment names
        segment_1 = str(annotation_table.loc[annotation_table["id"] == hit_name, "segment01"].iloc[0])
        segment_2 = str(annotation_table.loc[annotation_table["id"] == hit_name, "segment02"].iloc[0])
        # add the segment names and positions to the hits_with_positions data frame
        bottom_hits_with_positions.loc[index, "segment_1"] = str(segment_1)
        bottom_hits_with_positions.loc[index, "start_1"] = int(start_1)
        bottom_hits_with_positions.loc[index, "end_1"] = int(end_1)
        bottom_hits_with_positions.loc[index, "segment_2"] = str(segment_2)
        bottom_hits_with_positions.loc[index, "start_2"] = int(start_2)
        bottom_hits_with_positions.loc[index, "end_2"] = int(end_2)

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
                f"chr - {segment} {segment} 1 {length} grey\n"
            )
            chr_number += 1

    # Create the top_hits.txt file
    with open(f"{output_dir}/top_hits.txt", "w") as hits_file:
        # Iterate over the hits
        for index, row in top_hits_with_positions.iterrows():
            # Write the hit to the hits.txt file
            hits_file.write(
                f"{row['id']} {row['segment_1']} {int(row['start_1'])} {int(row['end_1'])}\n"
                f"{row['id']} {row['segment_2']} {int(row['start_2'])} {int(row['end_2'])}\n"
            )

    # Create the bottom_hits.txt file
    with open(f"{output_dir}/bottom_hits.txt", "w") as hits_file:
        # Iterate over the hits
        for index, row in bottom_hits_with_positions.iterrows():
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
            f"<ideogram>\n\n"

            "<spacing>\n"
            "default = 0.005r\n"
            "</spacing>\n\n"

            "radius    = 0.9r\n"
            "thickness = 20p\n"
            "fill      = yes\n\n"

            "</ideogram>\n\n"
            
            f"karyotype = {output_dir}/karyotype.txt\n\n"
            
            "show_ticks          = yes\n"
            "show_tick_labels    = yes\n\n"

            "<ticks>\n\n"

            "radius               = dims(ideogram,radius_outer)\n"
            "label_offset         = 5p\n"
            "orientation          = out\n"
            "label_multiplier     = 1\n"
            "color                = black\n\n"

            "<tick>\n"
            "spacing        = 100u\n"
            "size           = 8p\n"
            "thickness      = 2p\n"
            "color          = black\n"
            "show_label     = yes\n"
            "label_size     = 20p\n"
            "label_offset   = 5p\n"
            "format         = %d\n"
            "</tick>\n\n"

            "</ticks>\n\n"

            "<links>\n\n"

            "<link>\n"
           f"file             = {output_dir}/top_hits.txt\n"
            "radius           = 0.85r\n"
            "bezier_radius    = 0r\n"
            "color            = blue_a3\n"
            "ribbon           = yes\n"
            "flat             = yes\n"
            "</link>\n\n"

            "<link>\n"
           f"file             = {output_dir}/bottom_hits.txt\n"
            "radius           = 0.85r\n"
            "bezier_radius    = 0r\n"
            "color            = red_a3\n"
            "ribbon           = yes\n"
            "flat             = yes\n"
            "</link>\n\n"

            "</links>\n\n"

            "<image>\n"
            "<<include etc/image.conf>>\n"
            "</image>\n\n"

            "<<include etc/colors_fonts_patterns.conf>>\n"
            "<<include etc/housekeeping.conf>>\n"
        )


def main():
    # Parse command line arguments
    args = docopt(__doc__)
    input_file = args["<input_file>"]
    annotation_table = args["--annotation"]
    genome_file = args["--genome"]
    output_dir = args["--output"]

    # Check if output directory exists, create it if not
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read input file with pandas, rename unnamed column to id
    DESeq2_results = pd.read_csv(input_file, sep="\t")
    DESeq2_results = DESeq2_results.rename(columns={"Unnamed: 0": "id"})

    # Read annotation table
    annotation_table = pd.read_csv(annotation_table, sep="\t")

    # Parse genome segment names and sequences
    genome_dict = hp.parse_fasta(genome_file)

    # Create circos plot files
    make_circos_files(DESeq2_results, annotation_table, genome_dict, output_dir)


if __name__ == "__main__":
    main()
