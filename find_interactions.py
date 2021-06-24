import numpy as np
from operator import itemgetter

# import networkx
# from networkx.algorithms.components.connected import connected_components


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
    unclustered_list = []
    for (i, j), value in np.ndenumerate(binarry_array):
        if value:
            if not unclustered_list:
                unclustered_list.append([(i, j)])
            else:
                is_neighbour = False
                for cluster in unclustered_list:
                    tmp = []
                    for coordinate in cluster:
                        if coordinate == (i - 1, j) or coordinate == (i, j - 1):
                            tmp.append((i, j))
                            is_neighbour = True
                            break
                    if is_neighbour:
                        cluster.extend(tmp)
                if not is_neighbour:
                    unclustered_list.append([(i, j)])

    # ref of the following block of code:
    # https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements
    coordinate_list = []
    while unclustered_list:
        first, *rest = unclustered_list
        first = set(first)

        lf = -1
        while len(first) > lf:
            lf = len(first)

            rest2 = []
            for r in rest:
                if first.intersection(set(r)):
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2

        coordinate_list.append(first)
        unclustered_list = rest

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
        list of dictionaries containing the coordinate, height and width of each
        interaction region.
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
