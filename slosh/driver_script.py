# %%
# Import functions for analysis of SPLASH np.arrays
import variance_by_element as vbe
import find_interactions as fi

# %%
# Define file paths and viral features:
DIRECTORY = "/data/dessertlocal/projects/gl_iav-splash_freiburg"
INPUT = f"{DIRECTORY}/data/arrays"
RESULT = f"{DIRECTORY}/results/202108/20210824"
iav_segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]

# %%
# Read the arrays from file and put them into dictionaries
wt_d_repDir2Combinations, wt_d_combination2Array = vbe.read_arrays(INPUT, iav_segments)

# %%
# Extract arrays for testing, in this case interactions betwen NA and NP segments
NA_NP_arrays = wt_d_combination2Array["NA_NP"]
NA_NP_arrays

# %%
# Now we filter regions with readcounts greater than the mean of all values...
NA_NP_arrays_meanfiltered = [fi.mean_filter(array) for array in NA_NP_arrays]
NA_NP_arrays_meanfiltered

# %%
# ...and plot it.
iterator = 0
for array in NA_NP_arrays_meanfiltered:
    iterator += 1
    vbe.plot_heatmap(array, RESULT, f"NA_NP_{iterator}")

# %%
# Combine all regions in one array...
NA_NP_all_samples = fi.combine_filters(NA_NP_arrays_meanfiltered)
NA_NP_all_samples

# %%
# ...and plot it.
vbe.plot_heatmap(NA_NP_all_samples, RESULT, f"NA_NP_combined")

# %%
NA_NP_coordinates = fi.extract_coordinates(NA_NP_all_samples)

# %%
NA_NP_regions = fi.extract_regions(NA_NP_coordinates)
NA_NP_regions

# %%
NA_NP_means = [fi.readcounts_to_means(NA_NP_regions, array) for array in NA_NP_arrays]
fi.format_to_table(NA_NP_means, output_path=f"{RESULT}/NA_NP_interactions.csv")

# %%
