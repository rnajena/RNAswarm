# %%
# Import functions for analysis of SPLASH np.arrays
import variance_by_element


# Some scripts to analyse the variances across arrays:
# %%
# Define file paths and viral features
DIRECTORY = "/data/dessertlocal/projects/gl_iav-splash_freiburg"
INPUT = f"{DIRECTORY}/data/arrays"
RESULT = f"{DIRECTORY}/results/202104/20210427"
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


# Some scripts for direclty ploting the read-count arrays
# %%
# Define file paths and viral features
DIRECTORY = "/data/dessertlocal/projects/gl_iav-splash_freiburg"
INPUT = f"{DIRECTORY}/data/arrays"
RESULT = f"{DIRECTORY}/results/202104/20210429"
iav_segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]


# %%
# Read the arrays from file and put them into dictionaries
wt_d_repDir2Combinations, wt_d_combination2Array = variance_by_element.read_arrays(
    INPUT, iav_segments
)

# %%
# Take a look at the dictionary
print(wt_d_combination2Array)


# %%
# Extract just intercations between segments PB1 and PB2
l_PB1_PB2 = wt_d_combination2Array["PB1_PB2"]
print(l_PB1_PB2)


# %%
# Plot histograms for the three arrays of PB1 and PB2
i = 1
for array in l_PB1_PB2:
    variance_by_element.plot_heatmap(array, f"{RESULT}/PB1_PB2_0{i}")
    i += 1
