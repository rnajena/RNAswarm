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


def extract_coordinates(binarry_array):
    """ Return a list of lists containing coordinates for True values on a np.array, 
    clustering those if they are neighbouring values.

    Parameters
    ----------
    binarry_array: np.aray with binary values (True or False)

    Returns
    -------
    res: coordinate_list
        coordinates in the form of a list of tuples grouped if values are neighbouring
    """
    # Iterating through array and clustering together neighboring True values on the
    # same line
    coordinate_list = []
    for (i, j), value in np.ndenumerate(binarry_array):
        if value == True:
            if coordinate_list == []:
                coordinate_list.append([(i, j)])
            elif coordinate_list[-1][-1] == (i, j - 1):
                coordinate_list[-1].append((i, j))
            else:
                coordinate_list.append([(i, j)])

    # Iterating through array and clustering together neighboring True values on the
    # same column
    for (i, j), value in np.ndenumerate(binarry_array):
        if value == True:
            for cluster in coordinate_list:
                for coordinate in cluster:
                    if (i - 1, j) == coordinate:
                        for cluster in coordinate_list:
                            if (i, j) in cluster:
                                cluster.remove((i, j))
                        cluster.append((i, j))
            if [] in coordinate_list:
                coordinate_list.remove([])

    return coordinate_list


def extract_regions(coordinate_list):
    """ Return a list of dictionaries containing the coordinates for a rectangle that
    contains the interactions of a given cluster.

    Parameters
    ----------
    coordinate_list: coordinates in the form of a list of tuples grouped if values are
    neighbouring

    Returns
    -------
    res: regions_list
        list of dictionaries containing the coordinates for the four extremities of the
        triangle.
    """
    regions_dict = {}
    for region_id in range(len(coordinate_list)):
        i_sorted_coord_l = sorted(coordinate_list[region_id], key=itemgetter(0))
        j_sorted_coord_l = sorted(coordinate_list[region_id], key=itemgetter(1))
        height = i_sorted_coord_l[-1][0] - i_sorted_coord_l[0][0]
        width = j_sorted_coord_l[-1][1] - j_sorted_coord_l[0][1]
        if height == 0:
            height = 1
        if width == 0:
            width = 1
        regions_dict[region_id] = {
            "coordinate": (i_sorted_coord_l[0][0], j_sorted_coord_l[0][1]),
            "height": height,
            "width": width,
        }

    return regions_dict
