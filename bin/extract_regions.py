#!/usr/bin/env python3

"""interaction_finder.py

Takes an arbitrary number of trns files, finds and merges interactions, outputing an annotation table.

Usage:
  interaction_finder.py -g <genome> -i <input_file> -o <output_folder>
  interaction_finder.py -g <genome> -i <input_file> -o <output_folder> -m <min_components> -M <max_components>
  interaction_finder.py -g <genome> -i <input_file> -o <output_folder> --ignore_intra
  interaction_finder.py -g <genome> -i <input_file> -o <output_folder> -m <min_components> -M <max_components> --ignore_intra

Options:
  -h --help                             Show this screen.
  -g --genome=<genome>                  The genome filepath.
  -i --input=<input_file>               The input trns.txt file.
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
import scipy.stats as stats
import pandas as pd
import sklearn.mixture as mix
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle


def normalize_array(array, max_value=200000, mode="number_of_data_points", round=True):
    """
    Normalize an array to a range of 0 to max_value (default 200000).

    Parameters
    ----------
    array : array-like
        The array to normalize.

    round : bool, optional
        Whether to round the values to integers. The default is True.

    max_value : int
        The maximum value to normalize to.

    mode : str
        The mode to normalize the array in. Options are "peak_height" and "number_of_data_points".

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
        if round:
            if np.sum(array) < max_value:
                raise Exception("Number of datapoints is too small")
            else:
                return np.rint((array / np.sum(array)) * max_value)
        else:
            if np.sum(array) < max_value:
                raise Exception("Number of datapoints is too small")
            else:
                return (array / np.sum(array)) * max_value
    else:
        raise ValueError("Invalid mode")


def combine_arrays(arrays, max_value=200000, mode="number_of_data_points"):
    """
    Combine an arbitrary number of arrays into a single array, by normalizing each array to the same number of data points, summing them, and then rounding to the nearest integer.

    Parameters
    ----------
    arrays : list
        A list of arrays to combine.

    normalize : bool
        Whether to normalize the arrays to the same number of data points.

    max_value : int
        The maximum value to normalize to.

    mode : str
        The mode to normalize the array in. Options are "peak_height" and "number_of_data_points".

    Returns
    -------
    array-like
        The combined array.
    """
    normalized_arrays = []
    for array in arrays:
        normalized_arrays.append(
            normalize_array(array, max_value=max_value, mode=mode, round=False)
        )
    return np.rint(np.sum(normalized_arrays, axis=0))


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


def draw_ellipse(position, covariance, ax=None):
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
        # This is the case when the covariance_type is 'full'
        U, s, Vt = np.linalg.svd(covariance)
        angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
        width, height = 2 * np.sqrt(s)
    elif covariance.shape == (1, 1):
        # This is the case when the covariance_type is 'tied'
        angle = 0
        width, height = 2 * np.sqrt(covariance)
    elif covariance.shape == (2,):
        # This is the case when the covariance_type is 'diag'
        angle = 0
        width, height = 2 * np.sqrt(covariance)
    elif covariance.shape == ():
        # This is the case when the covariance_type is 'spherical'
        angle = 0
        width = height = 2 * np.sqrt(covariance)

    # Draw the Ellipse
    nsig = 1
    ax.add_patch(
        Ellipse(
            position,
            nsig * width,
            nsig * height,
            angle=angle,
            edgecolor="black",
            facecolor="none",
            antialiased=True
        )
    )


def plot_gmm(interaction_matrix, gmm, combination, output_file=None, ax=None, label_components=False):
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
    component = 0
    for pos, covar in zip(gmm.means_, gmm.covariances_):
        if label_components:
            plt.text(
                pos[0],
                pos[1],
                component,
                horizontalalignment="center",
                verticalalignment="center",
                color="black",
                fontsize=10,
            )
        draw_ellipse(pos, covar, ax=ax)
        component += 1
    if output_file:
        plt.savefig(output_file)
        plt.clf()
    else:
        plt.show()


def fit_optimal_gmm(
    density_array,
    min_components,
    max_components,
    max_iter=500,
    expected_delta=0.000001,
    get_all_gmms=False,
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
    if get_all_gmms:
        return gmm_dict
    else:
        return gmm_dict[components]


def calculate_individual_pdfs(gmm, density_array, weights=None):
    """
    Calculate the probability density function for each component of the Gaussian Mixture Model.

    Parameters
    ----------
    gmm : GaussianMixture
        The Gaussian Mixture Model to calculate the probability density function for.

    density_array : array-like
        The array to calculate the probability density function for.
    
    weights : array-like
        The weights of each component of the Gaussian Mixture Model.

    Returns
    -------
    pdfs : array-like
        The probability density function for each component of the Gaussian Mixture Model.
    """
    pdfs = {}
    means, covariances = gmm.means_.copy(), gmm.covariances_.copy()
    if weights is None:
        for component in range(gmm.n_components):
            pdfs[component] = stats.multivariate_normal.pdf(
                density_array, mean=means[component], cov=covariances[component]
            )
    else:
        for component in range(gmm.n_components):
            pdfs[component] = weights[component] * stats.multivariate_normal.pdf(
                density_array, means[component], covariances[component]
            )
    return pdfs


def calculate_residual_log_likelihoods(gmm, density_array, rebalance_weights=False):
    """
    Calculate the residual log likelihood of the Gaussian Mixture Model, when
    the probability density function of each component is calculated using
    the weights of the other components.
    
    Parameters
    ----------
    gmm : GaussianMixture
        The Gaussian Mixture Model to calculate the residual log likelihood for.

    density_array : array-like
        The array to calculate the residual log likelihood for.

    rebalance_weights : bool
        Whether to rebalance the weights of the Gaussian Mixture Model.

    Returns
    -------
    residual_log_likelihoods : array-like
        The residual log likelihood of the Gaussian Mixture Model.
    """
    residual_log_likelihoods = {}
    if rebalance_weights:
        for component in range(gmm.n_components):
            weights = gmm.weights_.copy()
            weights[component] = 0
            weights = weights / np.sum(weights)
            pdfs = calculate_individual_pdfs(gmm, density_array, weights=weights)
            residual_log_likelihoods[component] = np.sum(np.log(np.sum([pdfs[other_component] for other_component in pdfs.keys() if other_component != component], axis=0)))
    else:
        weights = gmm.weights_.copy()
        pdfs = calculate_individual_pdfs(gmm, density_array, weights=weights)
        for component in range(gmm.n_components):
            # Calculate the residual log likelihood for each component, except the current one
            residual_log_likelihoods[component] = np.sum(np.log(np.sum([pdfs[other_component] for other_component in pdfs.keys() if other_component != component], axis=0)))
    # Calculate the residual log likelihood for the whole model
    residual_log_likelihoods["total"] = np.sum(np.log(np.sum(list(pdfs.values()), axis=0)))
    residual_log_likelihoods["sklearn"] = np.sum(gmm.score_samples(density_array))
    return residual_log_likelihoods


def calculate_individual_log_likelihoods(gmm, density_array, rebalance_weights=False, refit_gmm=False):
    """
    Calculate the log likelihood of each component of the Gaussian Mixture Model.

    Parameters
    ----------
    gmm : GaussianMixture
        The Gaussian Mixture Model to calculate the log likelihood for.

    density_array : array-like
        The array to calculate the log likelihood for.

    rebalance_weights : bool
        Whether to rebalance the weights of the Gaussian Mixture Model.

    Returns
    -------
    log_likelihoods : array-like
        The log likelihood of each component of the Gaussian Mixture Model.
    """
    log_likelihoods = {}
    if refit_gmm:
        refitted_gmms = {}
        gmm_loglikelihood = np.sum(gmm.score_samples(density_array))
        print(f"{gmm.n_components}")
        for component in range(gmm.n_components):
            print(f"Refitting GMM with component {component} removed")
            # Refit the Gaussian Mixture Model, excluding the current component
            weights = np.delete(gmm.weights_.copy(), component) / np.sum(np.delete(gmm.weights_.copy(), component))
            means = np.delete(gmm.means_.copy(), component, axis=0)
            precisions = np.delete(gmm.precisions_.copy(), component, axis=0)
            gmm_refitted = mix.GaussianMixture(
                n_components=gmm.n_components - 1,
                max_iter=gmm.max_iter,
                covariance_type="full",
                init_params="k-means++",
                means_init=means,
                precisions_init=precisions,
                weights_init=weights,
                warm_start=True,
            ).fit(density_array)
            # Calculate the log likelihood of the current component
            log_likelihoods[component] = gmm_loglikelihood - np.sum(gmm_refitted.score_samples(density_array))
            refitted_gmms[component] = gmm_refitted
        return log_likelihoods, refitted_gmms
    else:
        residual_log_likelihoods = calculate_residual_log_likelihoods(
            gmm, density_array, rebalance_weights=rebalance_weights
        )
        for component in range(gmm.n_components):
            log_likelihoods[component] = residual_log_likelihoods["total"] - residual_log_likelihoods[component]
        return log_likelihoods


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


def parse_overlaping_elipses(gmm):
    """ """
    means = gmm.means_
    covariances = gmm.covariances_
    overlaps = []
    for i, mean_i in enumerate(means):
        for j, mean_j in enumerate(means):
            if i != j:
                if np.linalg.norm(mean_i - mean_j) < np.sqrt(
                    covariances[i][0][0] + covariances[j][0][0]
                ):
                    overlaps.append((i, j))
    return overlaps


def parse_rectangular_regions(gmm, sigma, combination, output_file=None):
    """
    Parse the regions that are covered by the Gaussian Mixture Model's
    components and extract the coordinates of a rectangle that covers the region.

    Parameters
    ----------
    gmm: GaussianMixture
        The Gaussian Mixture Model to parse.

    sigma: float
        The number of standard deviations to use to calculate the rectangle.

    combination: tuple
        The combination of the two segments being plotted.

    output_file: str
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
                "end01": int(math.ceil(x_max)),
                "segment02": combination[0],
                "start02": int(y_min),
                "end02": int(math.ceil(y_max)),
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

    Parameters
    ----------
    annotation_table : pandas.DataFrame
        The annotation table containing the regions to count the interactions for.

    trns_file : str
        The path to the segemehl trns file.

    intra_combinations : bool
        Whether to include the intra-segment combinations.
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
