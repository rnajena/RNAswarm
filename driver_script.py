# %%
# Import functions for analysis of SPLASH np.arrays
import numpy as np
import variance_by_element as vbe
import peak_finding as pf

# %%
# Define file paths and viral features:
DIRECTORY = "/data/dessertlocal/projects/gl_iav-splash_freiburg"
INPUT = f"{DIRECTORY}/data/arrays"
RESULT = f"{DIRECTORY}/results/202106/20210611"
iav_segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]

# %%
# Read the arrays from file and put them into dictionaries
wt_d_repDir2Combinations, wt_d_combination2Array = vbe.read_arrays(INPUT, iav_segments)

# %%
# Extract one of the arrays for testing
sample01_NA_NP = wt_d_combination2Array["NA_NP"][0]

# %%
# Extract interaction regions from the arrays...
sample01_NA_NP_binary = pf.std_deviation_filter(sample01_NA_NP)

# %%
# ...and plot it.
vbe.plot_heatmap(sample01_NA_NP_binary, RESULT, "sample01-NA_NP_01-threshold_std")

# %%
# Or we can use a different filter to extract interaction regions from the arrays...
sample01_NA_NP_mean_filtered = pf.mean_filter(sample01_NA_NP)

# %%
# ...and plot it.
vbe.plot_heatmap(
    sample01_NA_NP_mean_filtered, RESULT, "sample01-NA_NP_01-threshold_mean"
)

# %%
# We can also plot only cells with value smaller than 10, for instance...
sample01_NA_NP_10_filtered = pf.arbitrary_filter(sample01_NA_NP, 10)

# %%
# ...and plot it.
vbe.plot_heatmap(sample01_NA_NP_10_filtered, RESULT, "sample01-NA_NP_01-threshold_10")

# %%
# As a test for the extract_regions functions, I'll create a binarry array:
test_bin_array = np.array(
    (
        (True, True, False, False, False),
        (False, False, False, True, False),
        (True, True, True, False, True),
        (True, False, True, False, True),
        (True, False, False, False, False),
    )
)

# %%
regions_test = pf.extract_coordinates(test_bin_array)
print(regions_test)
# %%

