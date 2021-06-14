import numpy as np


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
    """
    
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
                        coordinate_list.remove([(i, j)])
                        cluster.append((i, j))

    return coordinate_list


def extract_regions(coordinate_list):
    """
    """
    regions_dict = {}
    for region_id in range(len(coordinate_list)):
        regions_dict[region_id] = {
            "start": coordinate_list[region_id][0],
            "end": coordinate_list[region_id][-1],
        }
