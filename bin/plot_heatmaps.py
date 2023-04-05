#!/usr/bin/env python3

"""plot_heatmaps.py

Usage:
    plot_heatmaps.py <trns_file> <trns_file>... -g <genome> [-a <annotation_table> --intra_only] -o <output_folder>
    plot_heatmaps.py <trns_file> -g <genome> [-a <annotation_table> --intra_only] -o <output_folder>


Options:
    -h --help                                 Show this screen.
    <trns_file>                               Path to trns files
    -g --genome=<genome>                      The genome filepath.
    -o --output=<output_folder>               The output folder.
    -a --annotation_table=<annotation_table>  The annotation table filepath.
    --intra_only                              Only plot intra-segment interactions.

"""

from docopt import docopt
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import helper as hp
import trns_handler as th
import array_handler as ah


def plot_heatmaps(
    merged_combination_arrays, plots_folder, colour_palette="PiYG", regions=None
):
    """
    Plot heatmaps for manual fitting of GMMs

    Parameters
    ----------
    merged_combination_arrays : dict
        Dictionary of combination arrays
    plots_folder : str
        Path to folder where plots should be saved
    colour_palette : str
        Colour palette to use for plotting

    Returns
    
    -------
    None
    """
    for combination in merged_combination_arrays:
        # Plot raw heatmap
        plot_heatmap(
            merged_combination_arrays[combination],
            plots_folder,
            colour_palette,
            combination,
            colorbar_label="Read counts",
            regions=regions,
        )

        # Plot log10 transformed data
        plot_heatmap(
            np.log10(merged_combination_arrays[combination] + 1),
            plots_folder,
            colour_palette,
            combination,
            colorbar_label="log10(read counts + 1",
            suffix="_log10",
            regions=regions,
        )


def plot_heatmap(
    combination_array,
    plots_folder,
    colour_palette,
    combination,
    regions=None,
    colorbar_label="Read counts",
    suffix="",
):
    """
    Plot heatmap from a given combination.

    Parameters
    ----------
    merged_combination_arrays : dict
        Dictionary of combination arrays
    plots_folder : str
        Path to folder where plots should be saved
    colour_palette : str
        Colour palette to use for plotting
    combination : tuple
        Tuple of combination
    regions : pandas.DataFrame
        A table containing the annotations for the rectangular regions.
    colorbar_label : str
        Label for colorbar

    Returns
    -------
    None
    """
    ax = plt.gca()
    plt.imshow(combination_array, cmap=colour_palette)
    plt.xticks(
        np.arange(0, combination_array.shape[1], 25),
        np.arange(0, combination_array.shape[1], 25),
        rotation=90,
    )
    plt.yticks(
        np.arange(0, combination_array.shape[0], 25),
        np.arange(0, combination_array.shape[0], 25),
    )
    plt.tick_params(
        axis="both", which="major", labelsize=3, labeltop=True, labelright=True
    )
    plt.grid(which="major", axis="x", linestyle="-", linewidth="0.05", color="grey")
    plt.grid(which="major", axis="y", linestyle="-", linewidth="0.05", color="grey")
    plt.colorbar(label=colorbar_label)
    plt.xlabel(f"{combination[1]}")
    plt.ylabel(f"{combination[0]}")
    if regions is not None:
        suffix = f"{suffix}_annotated"
        for region in regions.itertuples():
            # check if region is in combination
            if (
                region.segment01 == combination[0]
                and region.segment02 == combination[1]
            ):
                ax.add_patch(
                    Rectangle(
                        (region.start02, region.start01),
                        region.end02 - region.start02,
                        region.end01 - region.start01,
                        linewidth=0.5,
                        edgecolor="silver",
                        facecolor="none",
                    )
                )
                # label the region
                ax.text(
                    int((region.start02 + region.end02) / 2),
                    int((region.start01 + region.end01) / 2),
                    region.id,
                    color="silver",
                    fontsize=3,
                    verticalalignment="center",
                    horizontalalignment="center",
                )
            elif (
                region.segment01 == combination[1]
                and region.segment02 == combination[0]
            ):
                ax.add_patch(
                    Rectangle(
                        (region.start01, region.start02),
                        region.end01 - region.start01,
                        region.end02 - region.start02,
                        linewidth=0.5,
                        edgecolor="silver",
                        facecolor="none",
                    )
                )
                # label the region
                ax.text(
                    int((region.start01 + region.end01) / 2),
                    int((region.start02 + region.end02) / 2),
                    region.id,
                    fontsize=3,
                    color="silver",
                    verticalalignment="center",
                    horizontalalignment="center",
                )
    plt.savefig(
        os.path.join(plots_folder, f"{combination[0]}_{combination[1]}{suffix}.pdf"),
        format="pdf",
    )
    plt.close()


def parse_annotation_table(annotation_table):
    """
    Parse annotation table

    Parameters
    ----------
    annotation_table : str
        Path to annotation table as a csv file

    Returns
    -------
    regions : pandas.DataFrame
        A table containing the annotations for the rectangular regions.
    """
    # First line is a header and the separators are tabs
    regions = pd.read_csv(annotation_table, sep="\t", header=0)
    # generate the peak pandas.DataFrame
    peaks = pd.DataFrame(
        {
            "id": regions["id"],
            "segment01": regions["segment01"],
            "01type": regions["01type"],
            "start01": regions["segment01_peak"] - 20,
            "end01": regions["segment01_peak"] + 20,
            "segment02": regions["segment02"],
            "02type": regions["02type"],
            "start02": regions["segment02_peak"] - 20,
            "end02": regions["segment02_peak"] + 20,
        }
    )
    return regions, peaks


def main():
    """
    Main function

    Returns
    -------
    None
    """
    args = docopt(__doc__)
    trns_files = args["<trns_file>"]
    genome_file = args["--genome"]
    output_folder = args["--output"]
    # Check if --intra_only is given
    intra_only = False
    if args["--intra_only"]:
        intra_only = True
    # check if --annotation_table is given
    if args["--annotation_table"]:
        annotation_table = args["--annotation_table"]
    else:
        annotation_table = None

    # Check if output folder exists, if not create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Process input files
    genome_dict = hp.parse_fasta(genome_file)
    combination_arrays = {}

    # Check if trns_files is a list or a single file
    if isinstance(trns_files, str):
        trns_file_name = os.path.basename(trns_file)
        trns_file_name = trns_file_name.split(".")[0]

        # Create and fill combination arrays
        combination_arrays[trns_file_name] = hp.make_combination_array(genome_dict, intra_only=intra_only)
        th.segemehlTrans2heatmap(trns_file, combination_arrays[trns_file_name], intra_only=intra_only)
        merged_combination_arrays = combination_arrays

    elif isinstance(trns_files, list):
        for trns_file in trns_files:
            # Get the name of the current trns file
            trns_file_name = os.path.basename(trns_file)
            trns_file_name = trns_file_name.split(".")[0]

            # Create and fill combination arrays
            combination_arrays[trns_file_name] = hp.make_combination_array(genome_dict, intra_only=intra_only )
            th.segemehlTrans2heatmap(trns_file, combination_arrays[trns_file_name], intra_only=intra_only)

        # Merge combination arrays
        merged_combination_arrays = ah.combine_arrays(
            combination_arrays, normalise_array=False
        )
    
    # Define colour palettes
    colour_palette = "gist_stern"

    # Create subfolder for colour palette if it does not exist
    raw_folder = os.path.join(output_folder, "raw")
    peak_folder = os.path.join(raw_folder, "peaks")
    if not os.path.exists(raw_folder):
        os.makedirs(raw_folder)
    if not os.path.exists(peak_folder):
        os.makedirs(peak_folder)
    if annotation_table is not None:
        regions, peaks = parse_annotation_table(annotation_table)
        plot_heatmaps(
            merged_combination_arrays,
            raw_folder,
            colour_palette=colour_palette,
            regions=regions,
        )
        plot_heatmaps(
            merged_combination_arrays,
            peak_folder,
            colour_palette=colour_palette,
            regions=peaks,
        )
    else:
        plot_heatmaps(
            merged_combination_arrays, output_folder, colour_palette=colour_palette
        )


if __name__ == "__main__":
    main()
