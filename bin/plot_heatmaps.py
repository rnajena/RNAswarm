#!/usr/bin/env python3

"""plot_heatmaps.py

Usage:
    plot_heatmaps.py <trns_file> <trns_file>... -g <genome> -o <output_folder>

Options:
    -h --help                             Show this screen.
    <trns_file>                           Path to trns files
    -g --genome=<genome>                  The genome filepath.
    -o --output=<output_folder>           The output folder.

"""

from docopt import docopt
import os
import numpy as np
import matplotlib.pyplot as plt
import helper as hp
import trns_handler as th
import array_handler as ah


def plot_heatmaps(merged_combination_arrays, plots_folder, colour_palette="PiYG"):
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
        plot_heatmap(merged_combination_arrays, plots_folder, colour_palette, combination)

        # Plot log2 transformed data
        plt.imshow(
            np.log2(merged_combination_arrays[combination] + 1), cmap=colour_palette
        )
        plt.xticks(
            np.arange(0, merged_combination_arrays[combination].shape[1], 25),
            np.arange(0, merged_combination_arrays[combination].shape[1], 25),
            rotation=90,
        )
        plt.yticks(
            np.arange(0, merged_combination_arrays[combination].shape[0], 25),
            np.arange(0, merged_combination_arrays[combination].shape[0], 25),
        )
        plt.tick_params(
            axis="both", which="major", labelsize=3, labeltop=True, labelright=True
        )
        plt.grid(which="major", axis="x", linestyle="-", linewidth="0.05", color="grey")
        plt.grid(which="major", axis="y", linestyle="-", linewidth="0.05", color="grey")
        plt.colorbar(label="log2(read counts + 1)")
        plt.xlabel(f"{combination[1]}")
        plt.ylabel(f"{combination[0]}")
        plt.savefig(
            os.path.join(plots_folder, f"{combination[0]}_{combination[1]}_log2.pdf"),
            format="pdf",
        )
        plt.close()

        # Plot log10 transformed data
        plt.imshow(
            np.log10(merged_combination_arrays[combination] + 1), cmap=colour_palette
        )
        plt.xticks(
            np.arange(0, merged_combination_arrays[combination].shape[1], 25),
            np.arange(0, merged_combination_arrays[combination].shape[1], 25),
            rotation=90,
        )
        plt.yticks(
            np.arange(0, merged_combination_arrays[combination].shape[0], 25),
            np.arange(0, merged_combination_arrays[combination].shape[0], 25),
        )
        plt.tick_params(
            axis="both", which="major", labelsize=3, labeltop=True, labelright=True
        )
        plt.grid(which="major", axis="x", linestyle="-", linewidth="0.05", color="grey")
        plt.grid(which="major", axis="y", linestyle="-", linewidth="0.05", color="grey")
        plt.colorbar(label="log10(read counts + 1)")
        plt.xlabel(f"{combination[1]}")
        plt.ylabel(f"{combination[0]}")
        plt.savefig(
            os.path.join(plots_folder, f"{combination[0]}_{combination[1]}_log10.pdf"),
            format="pdf",
        )
        plt.close()

def plot_heatmap(combination_arrays, plots_folder, colour_palette, combination):
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

    Returns
    -------
    None
    """
    plt.imshow(combination_arrays[combination], cmap=colour_palette)
    plt.xticks(
            np.arange(0, combination_arrays[combination].shape[1], 25),
            np.arange(0, combination_arrays[combination].shape[1], 25),
            rotation=90,
        )
    plt.yticks(
            np.arange(0, combination_arrays[combination].shape[0], 25),
            np.arange(0, combination_arrays[combination].shape[0], 25),
        )
    plt.tick_params(
            axis="both", which="major", labelsize=3, labeltop=True, labelright=True
        )
    plt.grid(which="major", axis="x", linestyle="-", linewidth="0.05", color="grey")
    plt.grid(which="major", axis="y", linestyle="-", linewidth="0.05", color="grey")
    plt.colorbar(label="read counts")
    plt.xlabel(f"{combination[1]}")
    plt.ylabel(f"{combination[0]}")
    plt.savefig(
            os.path.join(plots_folder, f"{combination[0]}_{combination[1]}.pdf"),
            format="pdf",
        )
    plt.close()


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
        plot_heatmaps(
            merged_combination_arrays, palette_folder, colour_palette=colour_palette
        )

if __name__ == "__main__":
    main()