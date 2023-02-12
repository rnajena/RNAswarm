import numpy as np


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


def get_peak_cell_from_annotation_table(combination_arrays, annotation, genome_dict):
    """
    Get the peak cell of the interaction matrix from the annotation table.

    Parameters
    ----------
    combination_arrays : dict
        A dictionary of arrays, with the keys being the combination of segments.

    annotation : pandas.DataFrame
        The annotation dataframe.

    genome_dict : dict
        A dictionary with keys being genome segments and values being the genome sequence.

    Returns
    -------
    dict
        A dictionary of peak cells, with the keys being the combination of segments.
    """
    combination = (annotation['segment01'], annotation['segment02'])
    combination_reverse = (annotation['segment02'], annotation['segment01'])
    start01 = int(annotation['start01'])
    end01 = int(annotation['end01'])
    start02 = int(annotation['start02'])
    end02 = int(annotation['end02'])
    combination_array = combination_arrays[combination]
    if combination in combination_arrays:
        # slice the array to the region of interest
        region_of_interest = combination_array[
            start01:end01,
            start02:end02,
        ]
        # get the peak cell
        segment01_peak, segment02_peak = np.unravel_index(region_of_interest.argmax(), region_of_interest.shape)
        valuePeak = region_of_interest[segment01_peak, segment02_peak]
    elif combination_reverse in combination_arrays:
        # slice the array to the region of interest
        region_of_interest = combination_array[
            start02:end02,
            start01:end01,
        ]
        # get the peak cell
        segment01_peak, segment02_peak = np.unravel_index(region_of_interest.argmax(), region_of_interest.shape)
        valuePeak = region_of_interest[segment01_peak, segment02_peak]
    else:
        raise ValueError("Combination not found")
    peak_dict = {
        'segment01_peak': segment01_peak + start01,
        'segment02_peak': segment02_peak + start02,
        'valuePeak': valuePeak
    }
    return peak_dict


def get_peak_cell_from_daniels_table(combination_arrays, annotation, genome_dict):
    """
    Get the peak cell of the interaction matrix from daniel's annotation table.

    Parameters
    ----------
    combination_arrays : dict
        A dictionary of arrays, with the keys being the combination of segments.

    annotation : pandas.DataFrame
        The annotation dataframe.

    Returns
    -------
    dict
        A dictionary of peak cells, with the keys being the combination of segments.
    """
    combination = (annotation['aSeq'], annotation['bSeq'])
    combination_reverse = (annotation['bSeq'], annotation['aSeq'])
    transposed_annotation = negative_to_positive_strand(genome_dict, annotation['aSeq'], annotation['cai'], annotation['caj'], annotation['bSeq'], annotation['cbi'], annotation['cbj'])
    if combination in combination_arrays:
        # slice the array to the region of interest
        region_of_interest = combination_arrays[combination][
            transposed_annotation[1]:transposed_annotation[2],
            transposed_annotation[4]:transposed_annotation[5],
        ]
        # get the peak cell
        aPeak, bPeak = np.unravel_index(region_of_interest.argmax(), region_of_interest.shape)
        valuePeak = region_of_interest[aPeak, bPeak]
        # convert the peak cell to the original coordinates
        aPeak += transposed_annotation[1]
        bPeak += transposed_annotation[4]
    elif combination_reverse in combination_arrays:
        # slice the array to the region of interest
        region_of_interest = combination_arrays[combination_reverse][
            transposed_annotation[4]:transposed_annotation[5],
            transposed_annotation[1]:transposed_annotation[2],
        ]
        # get the peak cell
        bPeak, aPeak = np.unravel_index(region_of_interest.argmax(), region_of_interest.shape)
        valuePeak = region_of_interest[aPeak, bPeak]
        # convert the peak cell to the original coordinates
        aPeak += transposed_annotation[4]
        bPeak += transposed_annotation[1]
    else:
        raise ValueError("Combination not found")
    peak_dict = {
        'aPeak': positive_to_negative_strand_point(genome_dict, annotation['aSeq'], aPeak),
        'bPeak': positive_to_negative_strand_point(genome_dict, annotation['bSeq'], bPeak),
        'valuePeak': valuePeak
    }
    return peak_dict

def positive_to_negative_strand_point(genome_dict, segment, position):
    """Transpose a point from the positive to the negative strand.

    Args:
        genome_dict (dict): A dictionary containing the genome.
        aSeq (str): The name of the segment.
        a (int): The position of the point on the segment.

    Returns:
        int:
            The position of the point on the negative strand.
    """
    # Get the length of the segment
    aLen = len(genome_dict[segment])
    # Get the position on the negative strand
    return aLen - position + 1

def negative_to_positive_strand(genome_dict, aSeq, cai, caj, bSeq, cbi, cbj):
    """Transpose the region of interest from the negative to the positive strand.

    Args:
        genome_dict (dict): A dictionary containing the genome.
        aSeq (str): The name of the first segment.
        bSeq (str): The name of the second segment.
        cai (int): The start position of the first segment.
        caj (int): The end position of the first segment.
        cbi (int): The start position of the second segment.
        cbj (int): The end position of the second segment.

    Returns:
        tuple: The transposed region of interest.
    """
    # Get the length of the first segment
    aLen = len(genome_dict[aSeq])

    # Get the length of the second segment
    bLen = len(genome_dict[bSeq])

    # Transpose the first segment
    caj_pos = aLen - cai + 1
    cai_pos = aLen - caj

    # Transpose the second segment
    cbj_pos = bLen - cbi + 1
    cbi_pos = bLen - cbj

    return aSeq, cai_pos, caj_pos, bSeq, cbi_pos, cbj_pos

def positive_to_negative_strand(genome_dict, aSeq, cai, caj, bSeq, cbi, cbj):
    """Transpose the region of interest from the positive to the negative strand.

    Args:
        genome_dict (dict): A dictionary containing the genome.
        aSeq (str): The name of the first segment.
        bSeq (str): The name of the second segment.
        cai (int): The start position of the first segment.
        caj (int): The end position of the first segment.
        cbi (int): The start position of the second segment.
        cbj (int): The end position of the second segment.

    Returns:
        tuple: The transposed region of interest.
    """
    # Get the length of the first segment
    aLen = len(genome_dict[aSeq])

    # Get the length of the second segment
    bLen = len(genome_dict[bSeq])

    # Transpose the first segment
    caj_neg = aLen - cai
    cai_neg = aLen + 1 - caj

    # Transpose the second segment
    cbj_neg = bLen - cbi
    cbi_neg = bLen + 1 - cbj

    return aSeq, cai_neg, caj_neg, bSeq, cbi_neg, cbj_neg
