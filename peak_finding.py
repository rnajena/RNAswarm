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


def extract_regions(binarry_array):
    """
    
    """
    # First iterate through all lines of the array
    
    _list = []
    for (i, j), value in np.ndenumerate(binarry_array):
        if value == True:
            if str(coordinate_list).endswith(f"({i}, {j - 1})]]"):
                coordinate_list[-1].append((i, j))
            else:
                coordinate_list.append([(i, j)])
    # Then iterate through the 
    for element in list:
        
    return coordinate_list

