#!/usr/bin/env python3

"""plot_heatmaps.py

Usage:
    plot_heatmaps.py <trns_file> <trns_file>... -g <genome> -o <output_folder>
    plot_heatmaps.py <trns_file> <trns_file>... -g <genome> -a <annotation_table> -o <output_folder>


Options:
    -h --help                                 Show this screen.
    <trns_file>                               Path to trns files
    -g --genome=<genome>                      The genome filepath.
    -o --output=<output_folder>               The output folder.
    -a --annotation_table=<annotation_table>  The annotation table filepath.

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

        # Plot log2 transformed data
        plot_heatmap(
            np.log2(merged_combination_arrays[combination] + 1),
            plots_folder,
            colour_palette,
            combination,
            colorbar_label="log2(read counts + 1)",
            suffix="_log2",
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
                    fontsize=5,
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
                    fontsize=5,
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
    # First line is header and separator is a comma
    regions = pd.read_csv(annotation_table, sep=",", header=0)
    return regions


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
    annotation_table = args["--annotation_table"]
    colour_palettes = [
        "gist_stern",
        "terrain",
    ]

    # Process input files
    genome_dict = hp.parse_fasta(genome_file)
    combination_arrays = {}

    for trns_file in trns_files:
        # Get the name of the current trns file
        trns_file_name = os.path.basename(trns_file)
        trns_file_name = trns_file_name.split(".")[0]

        # Create and fill combination arrays
        combination_arrays[trns_file_name] = hp.make_combination_array(genome_dict)
        th.segemehlTrans2heatmap(trns_file, combination_arrays[trns_file_name])

    # Normalise arrays and create density arrays for GMM fitting
    merged_combination_arrays = ah.combine_arrays(
        combination_arrays, normalise_array=False
    )

    for colour_palette in colour_palettes:
        # Create subfolder for colour palette
        palette_folder = os.path.join(output_folder, colour_palette)
        os.mkdir(palette_folder)
        if annotation_table is not None:
            regions = parse_annotation_table(annotation_table)
            plot_heatmaps(
                merged_combination_arrays,
                palette_folder,
                colour_palette=colour_palette,
                regions=regions,
            )
        else:
            plot_heatmaps(
                merged_combination_arrays, palette_folder, colour_palette=colour_palette
            )


if __name__ == "__main__":
    main()
