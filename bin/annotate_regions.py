#!/usr/bin/env python3

"""interaction_finder.py

Usage:
  interaction_finder.py -g <genome> -i <input_files> -o <output_folder>
  interaction_finder.py -g <genome> -i <input_files> -o <output_folder> -m <min_components> -M <max_components>
  interaction_finder.py -g <genome> -i <input_files> -o <output_folder> --ignore_intra
  interaction_finder.py -g <genome> -i <input_files> -o <output_folder> -m <min_components> -M <max_components> --ignore_intra

Options:
  -h --help                             Show this screen.
  -g --genome=<genome>                  The genome filepath.
  -i --input=<input_files>              The input trns.txt files.
  -o --output=<output_folder>           The output folder.
  -m --min_components=<min_components>  The minimum number of components to use
                                        for the Gaussian Mixture Model [default: 28].
  -M --max_components=<max_components>  The maximum number of components to use 
                                        for the Gaussian Mixture Model [default: 30].
  --ignore_intra                        Ignore intra-segment interactions.

"""

from docopt import docopt
import math
import helper as hp
import handle_chimeras as hc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.mixture as mix
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle


def normalize_array(array, max_value=2000000, mode="number_of_data_points"):
    """
    Normalize an array to a range of 0 to max_value (default 2000000).

    Parameters
    ----------
    array : array-like
        The array to normalize.

    max_value : int
        The maximum value to normalize to.

    mode : str
        The mode to normalize the array in. Options are "peak_height" and "total_sum".

    Returns
    -------
    array-like
        The normalized array.
    """
    if mode == "peak_height":
        if np.max(array) < 60:
            return array
        else:
            return np.around((array / np.max(array)) * max_value)
    elif mode == "number_of_data_points":
        return np.rint((array / np.sum(array)) * max_value)
    else:
        raise ValueError("Invalid mode")

def make_polled_matrix()

def convert_to_density_array(interaction_matrix):
    """
    Convert an array to a density array. For each cell on the interaction matrix,
    append the position of the cell to the array the number of times equal to the
    value of the cell.

    Parameters
    ----------
    interaction_matrix : array-like
        The interaction matrix to convert to a density array.

    Returns
    -------
    array-like
        The density array.
    """
    density_list = []
    # iterate over each i,j coordinate in the array
    for (y, x), value in np.ndenumerate(interaction_matrix):
        for i in range(int(value)):
            density_list.append((x, y))
    return np.array(density_list)


def plot_bic_scores(gmm_dict, output_folder=None):
    """
    Plot the BIC scores for each number of components.

    Parameters
    ----------
    gmm_dict : dict
        A dictionary of Gaussian Mixture Models, with the number of components as the key.
    output_folder : str
        The output folder to save the plot to.

    Returns
    -------
    None
    """
    bic_scores = []
    components = []
    for n_components, gmm in gmm_dict.items():
        bic_scores.append(gmm["bic"])
        components.append(n_components)
    plt.plot(components, bic_scores)
    plt.xlabel("n_components")
    plt.ylabel("BIC score")
    if output_folder:
        plt.savefig(output_folder + "bic_scores.png")
    else:
        plt.show()
    plt.close()


def draw_ellipse(position, covariance, ax=None, **kwargs):
    """Draw an ellipse with a given position and covariance

    Parameters
    ----------
    position : array-like
        The mean of the distribution

    covariance : array-like
        The covariance matrix of the distribution

    Other Parameters
    ----------------
    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    Returns
    -------
    matplotlib.patches.Ellipse
        The ellipse drawn.
    """
    ax = ax or plt.gca()

    # Convert covariance to principal axes
    if covariance.shape == (2, 2):
        U, s, Vt = np.linalg.svd(covariance)
        angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
        width, height = 2 * np.sqrt(s)
    else:
        angle = 0
        width, height = 2 * np.sqrt(covariance)

    # Draw the Ellipse
    nsig = 1
    ax.add_patch(Ellipse(position, nsig * width, nsig * height, angle=angle, **kwargs))


def plot_gmm(interaction_matrix, gmm, combination, filename, ax=None):
    """
    Plot the Gaussian Mixture Model.

    Parameters
    ----------
    interaction_matrix : array-like
        The interaction matrix to plot.

    gmm : GaussianMixture
        The Gaussian Mixture Model to plot.

    combination : str
        The combination of the two segments being plotted.

    filename : str
        The filename to save the plot to.

    ax : matplotlib.axes.Axes
        The axes object to plot the Gaussian Mixture Model into.

    Returns
    -------
    None
    """
    plt.imshow(np.log10(interaction_matrix + 1), cmap="PiYG")
    plt.colorbar(label="log10(counts + 1)")
    plt.xlabel(f"{combination[1]}")
    plt.ylabel(f"{combination[0]}")
    for pos, covar in zip(gmm.means_, gmm.covariances_):
        draw_ellipse(pos, covar, ax=ax, alpha=0.5, color="black")
    plt.savefig(filename)
    plt.clf()


def fit_optimal_gmm(
    density_array,
    min_components,
    max_components,
    max_iter=500,
    expected_delta=0.000001,
):
    """
    Using BIC score, fit a Gaussian Mixture Model to the array, and decide the optimal number of components.

    Parameters
    ----------
    """
    optimal_number_of_components = False
    components = min_components
    gmm_dict = {}
    while not optimal_number_of_components and components <= max_components:
        print(f"Fitting GMM with {components} components")
        gmm_dict[components] = {}
        gmm_dict[components] = mix.GaussianMixture(
            n_components=components,
            max_iter=max_iter,
            covariance_type="full",
            init_params="k-means++",
        ).fit(density_array)
        if len(gmm_dict) > 1:
            bic_delta = np.absolute(
                gmm_dict[components].bic(density_array)
                - gmm_dict[components - 1].bic(density_array)
            )
            if bic_delta < expected_delta:
                optimal_number_of_components = True
                break
            elif components == max_components:
                break
            else:
                components += 1
        else:
            components += 1
    return gmm_dict[components]


def fit_gmms(array_dict, min_components, max_components, max_value=2000000):
    """ """
    gmms_dict = {}
    for combination, array in array_dict.items():
        print(f"Fitting GMM for {combination}")
        normalized_array = normalize_array(array, max_value)
        density_array = convert_to_density_array(normalized_array)
        gmm = fit_optimal_gmm(density_array, min_components, max_components)
        gmms_dict[combination] = gmm
        print(f"Optimal number of components for {combination} is {gmm.n_components}")
    return gmms_dict


def parse_rectangular_regions(gmm, sigma, combination, output_file=None):
    """
    Parse the regions that are covered by the Gaussian Mixture Model's
    components and extract the coordinates of a rectangle that covers the region.

    Parameters
    ----------
    gmm : GaussianMixture
        The Gaussian Mixture Model to parse.

    sigma : float
        The number of standard deviations to use to calculate the rectangle.

    combination : tuple
        The combination of the two segments being plotted.

    output_file : str
        The output file to save the regions to.

    Returns
    -------
    pandas.DataFrame
        A table containing the annotations for the rectangular regions.
    """
    regions = []
    for pos, covar in zip(gmm.means_, gmm.covariances_):
        x, y = pos
        x_std, y_std = np.sqrt(covar[0][0]), np.sqrt(covar[1][1])
        x_min, x_max = x - sigma * x_std, x + sigma * x_std
        y_min, y_max = y - sigma * y_std, y + sigma * y_std
        regions.append(
            {
                "segment01": combination[1],
                "start01": int(x_min),
                "end01": math.ceil(x_max),
                "segment02": combination[0],
                "start02": int(y_min),
                "end02": math.ceil(y_max),
            }
        )
    if output_file:
        pd.DataFrame(regions).to_csv(output_file, mode="a", index=False)
        return pd.DataFrame(regions)
    return pd.DataFrame(regions)


def plot_regions(
    regions, interaction_matrix, gmm=None, plot_gmms=False, output_file=None
):
    """
    Plot the rectangular regions on top of a 2d interaction matrix.

    Parameters
    ----------
    regions : pandas.DataFrame
        The table containing the annotations for the rectangular regions.

    interaction_matrix : array-like
        The interaction matrix to plot the rectangles on top of.

    gmm : GaussianMixture
        The Gaussian Mixture Model to plot.

    plot_gmms : bool
        Whether to plot the Gaussian Mixture Model.

    output_file : str
        The output file to save the plot to.

    Returns
    -------
    None
    """
    ax = plt.gca()
    plt.imshow(interaction_matrix)
    plt.colorbar()
    for region in regions.itertuples():
        ax.add_patch(
            Rectangle(
                (region.start01, region.start02),
                region.end01 - region.start01,
                region.end02 - region.start02,
                linewidth=1,
                edgecolor="r",
                facecolor="none",
            )
        )
        if plot_gmms and gmm:
            for pos, covar in zip(gmm.means_, gmm.covariances_):
                draw_ellipse(pos, covar, ax=plt.gca(), alpha=0.5, color="black")
    if output_file:
        plt.savefig(output_file)
        plt.clf()
    else:
        plt.show()


def fill_count_table(interaction, count_table, annotation_table_dict):
    """
    Fills the count table with the number of interactions between the given
    regions.

    Parameters
    ----------
    interaction : list
        A list containing the two regions of the interaction and start and end positions.

    count_table : dict
        A dictionary containing the count table.

    annotation_table_dict : dict
        A dictionary containing the annotation table.

    Returns
    -------
    None
    """
    for index, row in annotation_table_dict.items():
        if row["segment01"] == interaction[0] and row["segment02"] == interaction[3]:
            if (
                row["start01"] <= interaction[1] <= row["end01"]
                or row["start01"] <= interaction[2] <= row["end01"]
                and row["start02"] <= interaction[4] <= row["end02"]
                or row["start02"] <= interaction[5] <= row["end02"]
            ):
                count_table[index] += 1
            else:
                continue
        elif row["segment02"] == interaction[0] and row["segment01"] == interaction[3]:
            if (
                row["start02"] <= interaction[1] <= row["end02"]
                or row["start02"] <= interaction[2] <= row["end02"]
                and row["start01"] <= interaction[4] <= row["end01"]
                or row["start01"] <= interaction[5] <= row["end01"]
            ):
                count_table[index] += 1
            else:
                continue


def make_count_table(
    annotation_table, trns_file, interaction_arrays, intra_combinations=False
):
    """
    Creates a count table from a given annotation table and segemehl trns file.
    """
    annotation_table_dict = annotation_table.to_dict(orient="index")
    count_table = {index: 0 for index in range(len(annotation_table))}
    with open(trns_file) as input_stream:
        for line in input_stream:
            line = line.strip().split()
            firstRead = line[0].split(",")
            secondRead = line[1].split(",")
            currentRow = hc.__extract_start_stop_segemehl(
                firstRead
            ) + hc.__extract_start_stop_segemehl(secondRead)
            interaction = hc.__check_interaction(currentRow, interaction_arrays)
            if intra_combinations:
                fill_count_table(interaction, count_table, annotation_table_dict)
            else:
                if interaction[0] != interaction[3]:
                    fill_count_table(interaction, count_table, annotation_table_dict)
                else:
                    continue
    return count_table


def main():
    # Parse the command line arguments
    arguments = docopt(__doc__)
    genome_file_path = arguments["--genome"]
    trns_file = arguments["--input"]
    output_folder = arguments["--output"]
    min_components = int(arguments["--min_components"])
    max_components = int(arguments["--max_components"])

    # Process input files
    genome_dict = hp.parse_fasta(genome_file_path)
    combination_array = hp.make_combination_array(genome_dict)

    # Fill the array with the number of reads that map to each combination
    hc.segemehlTrans2heatmap(trns_file, combination_array)

    # Fit GMMs to the arrays
    gmms_dict = fit_gmms(combination_array, min_components, max_components)

    # Parse the regions covered by the GMMs
    for combination, gmm_dict in gmms_dict.items():
        regions = parse_rectangular_regions(
            gmm_dict[combination],
            sigma=1,
            combination=combination,
            output_file=f"{output_folder}/{trns_file.split('/')[-1]}_annotation_table.csv",
        )
        plot_regions(
            regions,
            combination_array[combination],
            gmm=gmm_dict[combination],
            plot_gmms=True,
        )

    # Create a count table
    annotation_table = pd.read_csv(
        f"{output_folder}/{trns_file.split('/')[-1]}_annotation_table.csv"
    )
    interaction_arrays = hp.make_interaction_array(genome_dict)
    count_table = pd.DataFrame(
        make_count_table(
            annotation_table, trns_file, interaction_arrays, intra_combinations=True
        ),
        index=[0],
    )
    count_table.to_csv(
        f"{output_folder}/{trns_file.split('/')[-1]}_count_table.csv", index=False
    )


if __name__ == "__main__":
    main()
