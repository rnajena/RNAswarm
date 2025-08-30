import csv
import trns_handler as th


def __extract_complementary_region(line, interaction_arrays=None):
    """Extracts complementary region information from a ChimericFragments line.

    Parameters
    ----------
    line : list
        A list of fields from a ChimericFragments file line.
        name1,type1,ref1,strand1,pos1,name2,type2,ref2,strand2,pos2,ligationpoint_count,complementarity_pvalue,complementarity_score,complementarity1_left,complementarity1_right,complementarity1_length,complementarity2_left,complementarity2_right,complementarity2_length
    interaction_arrays : dict, optional
        A dictionary to store interaction data.

    Returns
    -------
    list
    """
    complementary_region = [line[2], int(line[13]), int(line[14]), line[7], int(line[16]), int(line[17])]
    if interaction_arrays:
        # Checks if the cobination is not in the interaction_arrays
        if (complementary_region[0], complementary_region[3]) not in interaction_arrays:
            # If not, it reverses the interaction
            complementary_region = complementary_region[3:] + complementary_region[0:3]
    return complementary_region


def chimericFragments2heatmap(cf_file, interaction_arrays, intra_only=False):
    """Parses the ChimericFragments file and fills the interaction_arrays

    Parameters
    ----------
    cf_file : str
        Path to the ChimericFragments file.
    interaction_arrays : dict
        Dictionary to fill with interaction data.
    intra_only : bool, optional
        If True, only plot intra-segment interactions.

    Returns
    -------
    None
    """
    with open(cf_file) as inputStream:
        for line in inputStream:
            with open(cf_file) as inputStream:
                reader = csv.reader(inputStream)
                header = next(reader)  # Skip header line if present
                for line in reader:
                    complementary_region = __extract_complementary_region(line)
                    ligationpoint_count = int(line[10])
                    if intra_only:
                        if complementary_region[0] == complementary_region[3]:
                            th.fill_heatmap(complementary_region, interaction_arrays, copies=ligationpoint_count, intra=True)
                    else:
                        if complementary_region[0] != complementary_region[3]:
                            th.fill_heatmap(complementary_region, interaction_arrays, copies=ligationpoint_count)
