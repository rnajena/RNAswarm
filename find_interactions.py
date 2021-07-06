import numpy as np
from operator import itemgetter


def std_deviation_filter(arr):
    """
    Given a nd.array, returns a binary nd.array in which the value for a cell is true if
    it is greater than the standard deviation.
    """
    filtered_arr = arr > np.std(arr.flatten())
    return filtered_arr


def mean_filter(arr):
    """
    Given a nd.array, returns a binary nd.array in which the value for a cell is true if
    it is greater than the mean.
    """
    filtered_arr = arr > np.mean(arr.flatten())
    return filtered_arr


def arbitrary_filter(arr, threshold):
    """
    Given a nd.array, returns a binary nd.array in which the value for a cell is true if
    it is smaller than a given int.
    """
    filtered_arr = arr < threshold
    return filtered_arr


# ref: https://stackoverflow.com/questions/20528328/numpy-logical-or-for-more-than-two-arguments
def combine_filters(binarry_arrays_list):
    combined_array = np.logical_or.reduce((np.array(binarry_arrays_list)))
    return combined_array


def extract_coordinates(binarry_array):
    """ Return a list of lists containing coordinates for True values on a np.array, 
    clustering those if they are neighboring values.

    Parameters
    ----------
    binarry_array: np.aray with binary values (True or False)

    Returns
    -------
    res: coordinate_list
        coordinates in the form of lists of tuples grouped if values are neighboring
    """
    coord_to_cluster = {}
    idx = 0
    for (i, j), value in np.ndenumerate(binarry_array):
        if value:
            previous_neighbors = [
                neighbor
                for coord in [(i - 1, j), (i, j - 1), (i - 1, j - 1), (i + 1, j - 1)]
                if coord[0] > 0 and coord[1] > 0 and coord[0] < binarry_array.shape[0]
            ]
            if any([coord_to_cluster[neighbor] for neighbor in previous_neighbors]):
                cluster_idx = set(
                    [coord_to_cluster[neighbor] for neighbor in previous_neighbors]
                )
                if len(cluster_idx) == 1:
                    coord_to_cluster[(i, j)] = cluster_idx[0]
                else:
                    lowest = min(cluster_idx)
                    for neighbor in previous_neighbors:
                        if binarry_array[neighbor]:
                            coord_to_cluster[neighbor] = lowest
                    coord_to_cluster[(i, j)] = lowest
            else:
                idx += 1
                for coord in previous_neighbors:
                    if coord_to_cluster[coord]:
                        coord_to_cluster[coord]
                coord_to_cluster[(i, j)] = idx

    return coord_to_cluster


def extract_regions(coord_to_cluster):
    """ Return a list of dictionaries containing the coordinates for a rectangle that
    contains the interactions of a given cluster.

    Parameters
    ----------
    coordinate_list: coordinates in the form of a list of tuples grouped if values are
    neighboring

    Returns
    -------
    res: regions_dict
        dict of dicts containing the coordinate, height and width of each
        interaction region.
    """
    regions_dict = {}
    clusters = set([idx for coordinate, idx in coord_to_cluster.items()])
    for cluster in clusters:
        cluster_coords = [
            coordinate for coordinate, idx in coord_to_cluster.items() if idx == cluster
        ]

        i_sorted_coord_l = sorted(cluster_coords, key=itemgetter(0))
        j_sorted_coord_l = sorted(cluster_coords, key=itemgetter(1))

        coordinate = (i_sorted_coord_l[0][0], j_sorted_coord_l[0][1])
        height = i_sorted_coord_l[-1][0] - coordinate[0] + 1
        width = j_sorted_coord_l[-1][1] - coordinate[1] + 1

        regions_dict[cluster] = {
            "coordinate": coordinate,
            "height": height,
            "width": width,
        }


def readcounts_to_means(regions_dict, readcount_aray):
    """ Returns a dictionary with interaction id as key and mean of the values inside
    the interaction region extracted with extract_regions function.

    Parameters
    ----------
    regions_dict: dict of dicts containing the coordinate, height and width of each
    interaction region.

    Returns
    -------
    res: readcount_dict
    """
    readcount_dict = {}
    for id, region in regions_dict.items():
        if region["width"] > 1 and region["height"] > 1:
            readcount_dict[id] = readcount_aray[
                region["coordinate"][0] : region["width"],
                region["coordinate"][1] : region["height"],
            ].mean()
        elif region["width"] > 1:
            readcount_dict[id] = readcount_aray[
                region["coordinate"][0],
                region["coordinate"][1] : region["coordinate"][1] + region["width"] - 1,
            ].mean()
        elif region["height"] > 1:
            readcount_dict[id] = readcount_aray[
                region["coordinate"][0] : region["coordinate"][0]
                + region["height"]
                - 1,
                region["coordinate"][1],
            ].mean()
        else:
            readcount_dict[id] = readcount_aray[
                region["coordinate"][0], region["coordinate"][1],
            ]
    return readcount_dict


def format_to_table(readcounts, sep=",", output_path=None):
    """ Returns ...

    Parameters
    ----------
    readcounts: ...
    path: ...

    Returns
    -------
    res: DEseq2_input
    """
    table = ""
    # we should check if all dicts inside readcounts have the same size
    for id in range(1, len(readcounts[0] + 1)):
        # We should generalize this to experimental settings with n replicates
        table = f"{table}\n{id}{sep}{readcounts[0][id]}{sep}{readcounts[1][id]}{sep}{readcounts[2][id]}"
    if output_path:
        with open(output_path, "w") as file:
            file.write(table)
    return table
