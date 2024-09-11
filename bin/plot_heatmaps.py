#!/usr/bin/env python3

"""plot_heatmaps.py

Usage:
    plot_heatmaps.py -t <trns_file>... -g <genome> [-a <annotation_table> --intra_only] -o <output_folder>
    plot_heatmaps.py -d <array_dir> -g <genome> [-a <annotation_table> --intra_only] -o <output_folder>

Options:
    -h --help                                 Show this screen.
    -t --trns_file=<trns_file>...             The trns files (space-separated).
    -d --array_dir=<array_dir>...             The array directories (space-separated).
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
import helper as hp


def plot_heatmaps(
    combination_array, plots_folder, color_palette="PiYG", regions=None
):
    """
    Plot heatmaps for manual fitting of GMMs

    Parameters
    ----------
    combination_array : dict
        Dictionary of combination arrays
    plots_folder : str
        Path to folder where plots should be saved
    color_palette : str
        color palette to use for plotting

    Returns

    -------
    None
    """
    for combination in combination_array.keys():
        # Plot raw heatmap
        plot_heatmap(
            combination_array[combination],
            plots_folder,
            color_palette,
            combination,
            colorbar_label="read count",
            regions=regions,
        )

        # Plot log10 transformed data
        plot_heatmap(
            np.log10(combination_array[combination] + 1),
            plots_folder,
            color_palette,
            combination,
            colorbar_label="log10(read count + 1)",
            suffix="_log10",
            regions=regions,
        )


def plot_heatmap(
    combination_array: np.ndarray,
    plots_folder: str,
    color_palette: str,
    combination: tuple,
    regions: pd.DataFrame = None,
    colorbar_label: str = "read count",
    suffix: str = "",
) -> None:
    """
    Plot heatmap from a given combination.

    Parameters
    ----------
    combination_array : np.ndarray
        Array representing the heatmap data
    plots_folder : str
        Path to folder where plots should be saved
    color_palette : str
        color palette to use for plotting
    combination : tuple
        Tuple of combination
    regions : pandas.DataFrame, optional
        A table containing the annotations for the rectangular regions, by default None
    colorbar_label : str, optional
        Label for colorbar, by default "Read counts"
    suffix : str, optional
        Suffix for the output file, by default ""

    Returns
    -------
    None
    """
    ax = plt.gca()
    plt.imshow(combination_array, cmap=color_palette)
    set_ticks_and_grid(combination_array, ax)
    plt.colorbar(label=colorbar_label)
    plt.xlabel(f"{combination[1]}")
    plt.ylabel(f"{combination[0]}")
    if regions is not None:
        suffix = f"{suffix}_annotated"
        annotate_regions(ax, combination, regions)
    else:
        suffix = f"{suffix}_unannotated"
    save_plot(plots_folder, combination, suffix)

    plt.close()


def set_ticks_and_grid(combination_array: np.ndarray, ax: plt.Axes) -> None:
    """
    Set ticks, grid and tick parameters for the given axis.

    Parameters
    ----------
    combination_array : np.ndarray
        Array representing the heatmap data
    ax : plt.Axes
        Axis to set ticks and grid

    Returns
    -------
    None
    """
    ax.set_xticks(np.arange(0, combination_array.shape[1], 25))
    ax.set_xticklabels(np.arange(0, combination_array.shape[1], 25), rotation=90)
    ax.set_yticks(np.arange(0, combination_array.shape[0], 25))
    ax.set_yticklabels(np.arange(0, combination_array.shape[0], 25))
    ax.tick_params(
        axis="both", which="major", labelsize=3, labeltop=True, labelright=True
    )
    ax.grid(which="major", axis="x", linestyle="-", linewidth="0.05", color="grey")
    ax.grid(which="major", axis="y", linestyle="-", linewidth="0.05", color="grey")


def annotate_regions(ax: plt.Axes, combination: tuple, regions: pd.DataFrame) -> None:
    """
    Annotate the regions on the given axis based on the combination.

    Parameters
    ----------
    ax : plt.Axes
        Axis to annotate regions on
    combination : tuple
        Tuple of combination
    regions : pandas.DataFrame
        A table containing the annotations for the rectangular regions

    Returns
    -------
    None
    """
    for region in regions.itertuples():
        if region.segment01 == combination[0] and region.segment02 == combination[1]:
            ax.add_patch(
                create_rectangle(
                    region.start01, region.start02, region.end01, region.end02
                )
            )
            ax.text(
                get_center(region.start02, region.end02),
                get_center(region.start01, region.end01),
                region.id,
                color="black",
                fontsize=4,
                verticalalignment="center",
                horizontalalignment="center",
            )
        elif region.segment01 == combination[1] and region.segment02 == combination[0]:
            ax.add_patch(
                create_rectangle(
                    region.start02, region.start01, region.end02, region.end01
                )
            )
            ax.text(
                get_center(region.start01, region.end01),
                get_center(region.start02, region.end02),
                region.id,
                color="black",
                fontsize=4,
                verticalalignment="center",
                horizontalalignment="center",
            )


def create_rectangle(start1: int, start2: int, end1: int, end2: int) -> Rectangle:
    """
    Create a Rectangle object with given start and end points.

    Parameters
    ----------
    start1 : int
        Start point of the first axis
    start2 : int
        Start point of the second axis
    end1 : int
        End point of the first axis
    end2 : int
        End point of the second axis

    Returns
    -------
    Rectangle
        Rectangle object with the given start and end points
    """
    return Rectangle(
        (start2, start1),
        end2 - start2,
        end1 - start1,
        linewidth=1,
        edgecolor="black",
        facecolor="none",
    )


def get_center(start: int, end: int) -> int:
    """
    Calculate the center of two given points.

    Parameters
    ----------
    start : int
        Start point
    end : int
        End point

    Returns
    -------
    int
        Center of the two given points
    """
    return int((start + end) / 2)


def save_plot(plots_folder: str, combination: tuple, suffix: str) -> None:
    """
    Save the plot in the specified folder with the given suffix.

    Parameters
    ----------
    plots_folder : str
        Path to folder where plots should be saved
    combination : tuple
        Tuple of combination
    suffix : str
        Suffix for the output file

    Returns
    -------
    None
    """
    plt.savefig(
        os.path.join(plots_folder, f"{combination[0]}_{combination[1]}{suffix}.pdf"),
        format="pdf",
    )
    plt.savefig(
        os.path.join(plots_folder, f"{combination[0]}_{combination[1]}{suffix}.svg"),
        format="svg",
    )



def prepare_arrays(
    array_dir=None, intra_only=True, genome_dict=None
):
    """
    Prepare arrays for plotting and merge them.

    Parameters
    ----------
    array_dir : str
        Path to a single array folder or a list of paths to multiple array folders
    intra_only : bool
        If True, only intra-chromosomal interactions are considered
    genome_dict : dict
        Dictionary of genome segments

    Returns
    -------
    merged_combination_arrays : dict
        Dictionary of combination arrays
    """
    combination_array = {}
    # Create a dictionary to store the combination arrays
    if array_dir is not None:
        # Get the name of the current array folder
        array_dir_name = os.path.basename(array_dir)
        array_dir_name = array_dir_name.split(".")[0]

        # Create and fill combination arrays
        combination_array = hp.make_combination_array(
            genome_dict, intra_only=intra_only
        )
        ah.import_combination_arrays(combination_array, array_dir)
    return combination_array


def main():
    """
    Main function

    Returns
    -------
    None
    """
    args = docopt(__doc__)


    # Get other arguments
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

    if args.get("--array_dir"):
        array_dir = args.get("--array_dir")
        combination_array = prepare_arrays(
            array_dir=array_dir,
            intra_only=intra_only,
            genome_dict=genome_dict,
        )

    # Define color palettes
    color_palette = "Greens"

    # Plot heatmaps
    if annotation_table is not None:
        regions = hp.parse_annotation_table(annotation_table)
        plot_heatmaps(
            combination_array,
            output_folder,
            color_palette=color_palette,
            regions=regions,
        )
    else:
        plot_heatmaps(
            combination_array, output_folder, color_palette=color_palette
        )


if __name__ == "__main__":
    main()
