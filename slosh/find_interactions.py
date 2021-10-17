import numpy as np
from operator import is_, itemgetter


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
    combined_array = np.logical_or.reduce(np.array(binarry_arrays_list))
    return combined_array


def is_valid_coordinate(coord, array):
    """Check if a coordinate is present in a given array.
    Parameters
    ----------
    binarry_array: np.aray with binary values (True or False)

    Returns
    -------
    res: coordinate_list
        coordinates in the form of lists of tuples grouped if values are neighboring
    """
    if coord[0] > 0 and coord[1] > 0 and coord[0] < array.shape[0]:
        return True


def extract_coordinates(binarry_array):
    """Return a list of lists containing coordinates for True values on a np.array,
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
    # Dictionary that has coordinate as key and cluster index as value
    idx = 0  # Iterator that holds cluster index value
    for (i, j), value in np.ndenumerate(binarry_array):
        # Iterates through each element of the array
        if value:  # If the element is True...
            neighbors = [
                (i - 1, j),
                (i, j - 1),
                (i - 1, j - 1),
                (i - 1, j + 1),
            ]  # Put all neighbors on a list
            neighbors_filtered = list(
                filter(
                    lambda coord: coord[0] > 0
                    and coord[1] > 0
                    and coord[1] < binarry_array.shape[1],
                    neighbors,
                )
            )  # Filter for just valid neighbors
            values_of_neighbors_filtered = [
                binarry_array[neighbor] for neighbor in neighbors_filtered
            ]
            cluster_idx = set([])
            if any(values_of_neighbors_filtered):
                # Case where there is a neighbor holding the value True
                for neighbor in neighbors_filtered:
                    if neighbor in coord_to_cluster.keys():
                        cluster_idx.add(coord_to_cluster[neighbor])
                cluster_idx = list(cluster_idx)
                cluster_idx.sort()
                # Creates a set that holds the indexes of all valid neighbors
                if len(cluster_idx) == 1:
                    # Case where all neighbors have the same cluster index
                    coord_to_cluster[(i, j)] = cluster_idx[0]
                else:
                    # Case where there are neighbors with different cluster indexes
                    lowest = min(cluster_idx)
                    for neighbor in neighbors_filtered:
                        if binarry_array[neighbor]:
                            coord_to_cluster[neighbor] = lowest
                    coord_to_cluster[(i, j)] = lowest
            else:
                # Case where there is no neighbor holding the value True
                idx += 1
                for neighbor in neighbors_filtered:
                    # Why are we checking previous neighbors here if they don't hold the value True?
                    if (
                        neighbor in coord_to_cluster.keys()
                        and coord_to_cluster[neighbor]
                    ):
                        coord_to_cluster[neighbor]
                coord_to_cluster[(i, j)] = idx

    return coord_to_cluster


def extract_regions(coord_to_cluster):
    """Return a list of dictionaries containing the coordinates for a rectangle that
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
    return regions_dict


def readcounts_to_means(regions_dict, readcount_aray):
    """Returns a dictionary with interaction id as key and mean of the values inside
    the interaction region as value.

    Parameters
    ----------
    regions_dict: dict of dicts containing the coordinate, height and width of each
    interaction region.

    Returns
    -------
    res: readcount_dict
    """
    readcount_dict = {}
    for idx, region in regions_dict.items():
        if region["width"] > 1 and region["height"] > 1:
            readcount_dict[idx] = readcount_aray[
                region["coordinate"][0] : region["coordinate"][0] + region["width"] - 1,
                region["coordinate"][1] : region["coordinate"][1]
                + region["height"]
                - 1,
            ].mean()
        elif region["width"] > 1:
            readcount_dict[idx] = readcount_aray[
                region["coordinate"][0],
                region["coordinate"][1] : region["coordinate"][1] + region["width"] - 1,
            ].mean()
        elif region["height"] > 1:
            readcount_dict[idx] = readcount_aray[
                region["coordinate"][0] : region["coordinate"][0]
                + region["height"]
                - 1,
                region["coordinate"][1],
            ].mean()
        else:
            readcount_dict[idx] = readcount_aray[
                region["coordinate"][0],
                region["coordinate"][1],
            ]
    return readcount_dict
