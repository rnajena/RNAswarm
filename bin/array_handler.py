import numpy as np
import os


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


def combine_arrays(combination_arrays, normalise_array=True, max_value=2000000, mode="number_of_data_points"):
    """
    Merges arrays for each combination of segments in the combination_arrays dictionary.

    Parameters
    ----------
    combination_arrays : dict
        A dictionary of arrays, with the keys being the combination of segments.

    normalise_array : bool, optional
        Whether to normalise the array. The default is True.

    max_value : int, optional
        The maximum value to normalise to. The default is 2000000.

    mode : str, optional
        The mode to normalise the array in. Options are "peak_height" and "number_of_data_points". The default is "number_of_data_points".

    Returns
    -------
    dict
        A dictionary of arrays, with the keys being the combination of segments.
    """
    merged_combination_arrays = {}
    for combinations in combination_arrays.values():
        for combination, array in combinations.items():
            if combination in merged_combination_arrays:
                merged_combination_arrays[combination] += array
            else:
                merged_combination_arrays[combination] = array
    if normalise_array:
        for combination, array in merged_combination_arrays.items():
            merged_combination_arrays[combination] = normalize_array(
                array, max_value=max_value, mode=mode
            )
    return merged_combination_arrays


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


def save_combination_arrays(combination_arrays, output_folder):
    """
    Save the combination arrays as a numpy array.

    Parameters
    ----------
    combination_arrays : dict
        A dictionary of arrays, with the keys being the combination of segments.

    output_folder : str
        The output folder to save the arrays to.
    """
    for combination, array in combination_arrays.items():
        output_file = os.path.join(output_folder, f"{combination[0]}-{combination[1]}.npy")
        np.save(output_file, array)

def import_combination_arrays(combination_arrays, input_folder):
    """
    Import the combination arrays as a numpy array.

    Parameters
    ----------
    combination_arrays : dict
        A dictionary of arrays, with the keys being the combination of segments.

    input_folder : str
        The input folder to import the arrays from.

    Returns
    -------
    dict
        A dictionary of arrays, with the keys being the combination of segments.
    """
    for combination, array in combination_arrays.items():
        combination_arrays[combination] = np.load(os.path.join(input_folder, combination + ".npy"))
    return combination_arrays