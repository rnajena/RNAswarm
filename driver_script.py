# %%
# Import functions for analysis of SPLASH np.arrays
import variance_by_element


# Some scripts to analyse the variances across arrays:
# %%
# Define file paths and viral features
DIRECTORY = "/data/dessertlocal/projects/gl_iav-splash_freiburg"
INPUT = f"{DIRECTORY}/data/arrays"
RESULT = f"{DIRECTORY}/results/202104/20210512"
iav_segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]


# %%
# Read the arrays from file and put them into dictionaries
wt_d_repDir2Combinations, wt_d_combination2Array = variance_by_element.read_arrays(
    INPUT, iav_segments
)


# %%
# Calculate the variances between the arrays
wt_d_comb2variances = variance_by_element.calculate_variances(wt_d_combination2Array)


# %%
#
wt_d_comb2characteristics = variance_by_element.save_characteristics(
    wt_d_comb2variances, RESULT
)
# print(wt_d_comb2characteristics)


# %%
wt_d_comb2variances["NP_NS"].flatten()


# %%
