# %%
# Import functions for analysis of SPLASH np.arrays
import variance_by_element as vbe
import find_interactions as fi

# %%
# Define file paths and viral features:
DIRECTORY = "/home/gabriellovate/Projects/gl_iav-splash_freiburg"
INPUT = f"{DIRECTORY}/data/arrays"
RESULT = f"{DIRECTORY}/results/202110/20211002"
iav_segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]

# %%
# Read the arrays from file and put them into dictionaries
wt_d_repDir2Combinations, wt_d_combination2Array = vbe.read_arrays(INPUT, iav_segments)

# %%
# Unpack the arrays and filter the regions with readcounts greater than the mean of all values...
wt_d_combination_2_array_filtered = {}
for segment, arrays in wt_d_combination2Array.items():
    for i in range(len(arrays)):
        if i == 0:
            wt_d_combination_2_array_filtered[segment] = [fi.mean_filter(arrays[i])]
        else:
            wt_d_combination_2_array_filtered[segment].append(fi.mean_filter(arrays[i]))

# %%
wt_d_combination_2_array_filtered

# %%
# Unpack the filtered arrays and plot them
for segment, arrays in wt_d_combination_2_array_filtered.items():
    for i in range(len(arrays)):
        vbe.plot_heatmap(arrays[i], RESULT, f"{segment}_{i}")

# %%
# Combine the arrays of each segment combination in one array
wt_d_combination_2_array_combined = {}
for segment, arrays in wt_d_combination_2_array_filtered.items():
    wt_d_combination_2_array_combined[segment] = fi.combine_filters(arrays)

# %%
# ...and plot it.
for segment, array in wt_d_combination_2_array_combined.items():
    vbe.plot_heatmap(array, RESULT, f"{segment}_combined")

# %%
wt_d_coordinates = {}
for segment, array in wt_d_combination_2_array_combined.items():
    wt_d_coordinates[segment] = fi.extract_coordinates(array)

# %%
wt_d_regions = {}
for segment, coordinates in wt_d_coordinates.items():
    wt_d_coordinates[segment] = fi.extract_coordinates(array)

# %%

wt_d_means = {}
for segment, regions in wt_d_regions.items():
    wt_d_means[segment] = [fi.readcounts_to_means(regions, array) for array in wt_d_combination_2_array_filtered[segment]]
    fi.format_to_table(wt_d_means[segment], output_path=f"{RESULT}/{segment}_interactions.csv")