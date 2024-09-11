#!/usr/bin/env python3

"""annotate_interactions.py

Takes an arbitrary number of trns files, finds and merges interactions, outputing an annotation table and the GMMs for each combination.

Usage:
    annotate_interactions.py -t <trns_file>... -g <genome> -o <output_file> [-m <min_components> -M <max_components> --step_size <step_size> --sigma <sigma>]
    annotate_interactions.py -d <array_dir> -g <genome> -o <output_file> [-m <min_components> -M <max_components> --step_size <step_size> --sigma <sigma>]
    annotate_interactions.py --filter --only_partner01 <partner_segment01> --only_partner02 <partner_segment02> -a <annotation_table> -o <output_file>

Options:
    -h --help                             Show this screen.
    -t --trns_file=<trns_file>...         The trns files (space-separated).
    -d --array_dir=<array_dir>...         The array directories (space-separated).
    -g --genome=<genome>                  The genome filepath.
    -o --output=<output_file>             The output folder.
    -m --min_components=<min_components>  The minimum number of components to use
                                          for the Gaussian Mixture Model [default: 70].
    -M --max_components=<max_components>  The maximum number of components to use.
                                          for the Gaussian Mixture Model [default: 80].
    --step_size=<step_size>               The step size to use for each iteration of the
                                          Gaussian Mixture Model optimization [default: 5].
    --sigma=<sigma>                       The number of standard deviations to use to [default: 1].

"""
import helper as hp
import trns_handler as th
import array_handler as ah
from docopt import docopt
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pandas as pd
import sklearn.mixture as mix
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
import os
import pickle


def plot_bic_scores(gmm_dict, density_array, output_folder=None):
    """
    Plot the BIC scores for each number of components.

    Parameters
    ----------
    gmm_dict : dict
        A dictionary of Gaussian Mixture Models, with the number of components as the key.
    density_array : array-like
        The array to plot the BIC scores for.
    output_folder : str
        The output folder to save the plot to.

    Returns
    -------
    None
    """
    bic_scores = []
    components = []
    for n_components, gmm in gmm_dict.items():
        bic_scores.append(gmm.bic(density_array))
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
            antialiased=True,
        )
    )


def plot_gmm(
    interaction_matrix,
    gmm,
    combination,
    output_file=None,
    ax=None,
    label_components=False,
):
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
    step_size=1,
):
    """
    Using BIC score, fit a Gaussian Mixture Model to the array, and decide the optimal number of components.
    """
    if min_components > max_components or min_components < 1:
        raise ValueError("min_components must be less than or equal to max_components and both greater than 0")

    optimal_number_of_components = False
    gmm_dict = {}

    for components in range(min_components, max_components + 1, step_size):
        print(f"Fitting GMM with {components} components")
        gmm = mix.GaussianMixture(
            n_components=components,
            max_iter=max_iter,
            covariance_type="full",
            init_params="k-means++",
        ).fit(density_array)
        gmm_dict[components] = gmm

        if components > min_components:  # Ensure there's a previous model to compare
            bic_delta = np.absolute(
                gmm.bic(density_array) - gmm_dict[components - step_size].bic(density_array)
            )
            if bic_delta < expected_delta:
                optimal_number_of_components = True
                break

    if get_all_gmms:
        return gmm_dict
    elif optimal_number_of_components:
        return gmm_dict[components]
    else:
        return gmm_dict[max_components]  # Return the model with max_components if no optimal found


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
            residual_log_likelihoods[component] = np.sum(
                np.log(
                    np.sum(
                        [
                            pdfs[other_component]
                            for other_component in pdfs.keys()
                            if other_component != component
                        ],
                        axis=0,
                    )
                )
            )
    else:
        weights = gmm.weights_.copy()
        pdfs = calculate_individual_pdfs(gmm, density_array, weights=weights)
        for component in range(gmm.n_components):
            # Calculate the residual log likelihood for each component, except the current one
            residual_log_likelihoods[component] = np.sum(
                np.log(
                    np.sum(
                        [
                            pdfs[other_component]
                            for other_component in pdfs.keys()
                            if other_component != component
                        ],
                        axis=0,
                    )
                )
            )
    # Calculate the residual log likelihood for the whole model
    residual_log_likelihoods["total"] = np.sum(
        np.log(np.sum(list(pdfs.values()), axis=0))
    )
    residual_log_likelihoods["sklearn"] = np.sum(gmm.score_samples(density_array))
    return residual_log_likelihoods


def calculate_individual_log_likelihoods(
    gmm, density_array, rebalance_weights=False, refit_gmm=False
):
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
            weights = np.delete(gmm.weights_.copy(), component) / np.sum(
                np.delete(gmm.weights_.copy(), component)
            )
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
            log_likelihoods[component] = gmm_loglikelihood - np.sum(
                gmm_refitted.score_samples(density_array)
            )
            refitted_gmms[component] = gmm_refitted
        return log_likelihoods, refitted_gmms
    else:
        residual_log_likelihoods = calculate_residual_log_likelihoods(
            gmm, density_array, rebalance_weights=rebalance_weights
        )
        for component in range(gmm.n_components):
            log_likelihoods[component] = (
                residual_log_likelihoods["total"] - residual_log_likelihoods[component]
            )
        return log_likelihoods


def fit_gmms(array_dict, min_components, max_components, max_value=2000000):
    """ """
    gmms_dict = {}
    for combination, array in array_dict.items():
        print(f"Fitting GMM for {combination}")
        normalized_array = ah.normalize_array(array, max_value)
        density_array = ah.convert_to_density_array(normalized_array)
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


def parse_rectangular_regions(gmm, combination, sigma=1, output_file=None, with_header=False):
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
        if with_header:
            pd.DataFrame(regions).to_csv(output_file, mode="a", index=False)
            return pd.DataFrame(regions)
        else:
            pd.DataFrame(regions).to_csv(output_file, mode="a", header=False, index=False)
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


def main():
    # Parse the command line arguments
    # interaction_finder.py -g <genome> -i <input_file> -o <output_folder> [-m <min_components> -M <max_components> --make_plots --ignore_intra]
    arguments = docopt(__doc__)
    genome_file_path = arguments["--genome"]
    array_folder = arguments["--array_dir"]
    output_folder = arguments["--output"]
    min_components = int(arguments["--min_components"])
    max_components = int(arguments["--max_components"])
    step_size = int(arguments["--step_size"])
    sigma =  float(arguments["--sigma"])

    # Process input files
    genome_dict = hp.parse_fasta(genome_file_path)
    combination_arrays = {}
    
    # Get the name of the current array folder
    array_folder_name = os.path.basename(array_folder)
    array_folder_name = array_folder_name.split(".")[0]

    # Import  arrays
    combination_arrays = hp.make_combination_array(genome_dict)
    ah.import_combination_arrays(combination_arrays, array_folder)

    density_arrays = {
        combination: ah.convert_to_density_array(combination_array)
        for combination, combination_array in combination_arrays.items()
    }

    # Use BIC score to fit optimal GMMs to the density arrays
    gmms_dict = {}
    for combination, density_array in density_arrays.items():
        print(f"Fitting GMMs for {combination}")
        gmms_dict[combination] = fit_optimal_gmm(
            density_array,
            min_components,
            max_components,
            max_iter=500,
            expected_delta=0.000001,
            get_all_gmms=False,
            step_size=step_size,
        )


    # Save the gmms dict to a pickle file
    gmms_pickle = f"{output_folder}/{output_folder}_gmms.pickle"
    with open(gmms_pickle, "wb") as handle:
        pickle.dump(gmms_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


    # Export regions as a table
    for combination, gmm in gmms_dict.items():
        # Parse the regions
        parse_rectangular_regions(
            gmms_dict[combination], combination, sigma, f"{output_folder}/{output_folder}.csv"
        )


    refit_gmms=False
    if refit_gmms:
        # Calculate the individual log likelihood for each component (refitting the model)
        individual_log_likelihoods_refitted = {}
        refitted_gmms_dict = {}

        for combination, gmms in gmms_dict.items():
            individual_log_likelihoods_refitted[combination] = {}
            refitted_gmms_dict[combination] = {}
            for number_of_components, gmm in gmms.items():
                (
                    individual_log_likelihoods_refitted[combination][number_of_components],
                    refitted_gmms_dict[combination][number_of_components],
                ) = calculate_individual_log_likelihoods(gmm, density_array, refit_gmm=True)

        # Create bar plots of the individual_log_likelihoods_refitted
        for (
            combination,
            individual_log_likelihoods,
        ) in individual_log_likelihoods_refitted.items():
            for number_of_components, log_likelihoods in individual_log_likelihoods.items():
                plt.bar(range(number_of_components), log_likelihoods)
                plt.title(
                    f"Individual log likelihoods for {combination} with {number_of_components} components"
                )
                plt.xlabel("Component")
                plt.ylabel("Log likelihood")
                # Save the plot on the combination folder
                combination_folder = f"{output_folder}/{combination}"
                plt.savefig(
                    f"{combination_folder}/{combination}_{number_of_components}_individual_log_likelihoods.pdf"
                )
                plt.close()

        # Plot the refitted GMMs for each combination, for each number of components, for each component
        for combination, refitted_gmms in refitted_gmms_dict.items():
            for number_of_components, refitted_gmm in refitted_gmms.items():
                # Save plots from each combination on a new folder
                combination_folder = f"{output_folder}/{combination}"
                if not os.path.exists(combination_folder):
                    os.makedirs(combination_folder)
                # Plot the GMM
                plot_gmm(
                    combination_arrays[combination],
                    refitted_gmm,
                    combination,
                    combination_folder,
                )


if __name__ == "__main__":
    main()
