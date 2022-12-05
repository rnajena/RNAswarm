#!/usr/bin/env python3


import extract_regions as er
import numpy as np
import matplotlib.pyplot as plt
import helper as hp
import handle_chimeras as hc
import os


# Load all files with trns.txt extension to a dictionary
trns_files = [
    "/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/20221202/02-mappings/segemehl/7wtintakt_S7_R1_001_trimmed.trns.txt",
    "/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/20221202/02-mappings/segemehl/8-wt_S8_R1_001_trimmed.trns.txt",
    "/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/20221202/02-mappings/segemehl/H-4_S14_R1_001_trimmed.trns.txt",
    "/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/20221202/02-mappings/segemehl/SC35M_WTWT_repli01_0120_trimmed.trns.txt",
    "/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/20221202/02-mappings/segemehl/SC35M_WTWT_repli01_1021_trimmed.trns.txt",
    "/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/20221202/02-mappings/segemehl/SC35M_WTWT_repli02_1120_trimmed.trns.txt",
    "/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/20221202/02-mappings/segemehl/SC35M_WTWT_repli03_1120_trimmed.trns.txt",
]
# Load all files with fasta extension to a dictionary
genome_file_path = "/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/data/schwemmle_group/genomes/SC35M_WTWT.fasta"


# Output folder
plots_folder = "/beegfs/ru27wav/Projects/gl_iav-splash_freiburg/results/schwemmle_group/plots_20221203"


# Process input files
genome_dict = hp.parse_fasta(genome_file_path)
combination_arrays = {}

for trns_file in trns_files:
    print(trns_file)
    # Get the name of the current trns file
    trns_file_name = os.path.basename(trns_file)
    print(trns_file_name)
    trns_file_name = trns_file_name.split(".")[0]
    print(trns_file_name)

    # Create and fill combination arrays
    combination_arrays[trns_file_name] = hp.make_combination_array(genome_dict)
    hc.segemehlTrans2heatmap(trns_file, combination_arrays[trns_file_name])


# Normalise arrays and create density arrays for GMM fitting
merged_combination_arrays = er.combine_arrays(combination_arrays, normalise_array=False)


def plot_heatmaps_for_manual_fitting(
    merged_combination_arrays, plots_folder, colour_palette="PiYG"
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
        print(combination)
        plt.imshow(merged_combination_arrays[combination], cmap=colour_palette)
        plt.xticks(
            np.arange(0, merged_combination_arrays[combination].shape[1], 25),
            np.arange(0, merged_combination_arrays[combination].shape[1], 25),
            rotation=90,
        )
        plt.yticks(
            np.arange(0, merged_combination_arrays[combination].shape[0], 25),
            np.arange(0, merged_combination_arrays[combination].shape[0], 25),
        )
        plt.tick_params(axis="both", which="major", labelsize=3, labeltop=True, labelright=True)
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
        plt.tick_params(axis="both", which="major", labelsize=3, labeltop=True, labelright=True)
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
        plt.tick_params(axis="both", which="major", labelsize=3, labeltop=True, labelright=True)
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



colour_palettes = [
    "gist_stern",
    "terrain",
]
for colour_palette in colour_palettes:
    # Create subfolder for colour palette
    palette_folder = os.path.join(plots_folder, colour_palette)
    os.mkdir(palette_folder)
    plot_heatmaps_for_manual_fitting(
        merged_combination_arrays, palette_folder, colour_palette=colour_palette
    )
