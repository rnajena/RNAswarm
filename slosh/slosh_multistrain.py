# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
# Import functions for analysis of SPLASH np.arrays
import io_ops as iops
import find_interactions as fi
import pandas as pd


# %%
# Define file paths and viral features:
DIRECTORY = "/home/ru27wav/Projects/gl_iav-splash_freiburg"
INPUT = f"{DIRECTORY}/data/arrays_test"
RESULT = f"{DIRECTORY}/results/test"
iav_segments = ["HA","M"]
strains = ["wt", "mut4xh4", "muts23"]


# %%
# Read the arrays from file and put them into dictionaries
wt_d_repDir2Combinations, wt_d_combinations2arrays = iops.read_arrays(
    f'{INPUT}/wt', iav_segments
)

muts23_d_repDir2Combinations, muts23_d_combinations2arrays = iops.read_arrays(
    f'{INPUT}/muts23', iav_segments
)

mut4xh4_d_repDir2Combinations, mut4xh4_d_combinations2arrays = iops.read_arrays(
    f'{INPUT}/mut4xh4', iav_segments
)


# %%
# Unpack the arrays and filter the regions with readcounts greater than the mean of all values...
wt_d_combinations2arrays_filtered = {}
for combination, arrays in wt_d_combinations2arrays.items():
    for i in range(len(arrays)):
        if i == 0:
            wt_d_combinations2arrays_filtered[combination] = [fi.mean_filter(arrays[i])]
        else:
            wt_d_combinations2arrays_filtered[combination].append(
                fi.mean_filter(arrays[i])
            )

muts23_d_combinations2arrays_filtered = {}
for combination, arrays in muts23_d_combinations2arrays.items():
    for i in range(len(arrays)):
        if i == 0:
            muts23_d_combinations2arrays_filtered[combination] = [fi.mean_filter(arrays[i])]
        else:
            muts23_d_combinations2arrays_filtered[combination].append(
                fi.mean_filter(arrays[i])
            )

mut4xh4_d_combinations2arrays_filtered = {}
for combination, arrays in mut4xh4_d_combinations2arrays.items():
    for i in range(len(arrays)):
        if i == 0:
            mut4xh4_d_combinations2arrays_filtered[combination] = [fi.mean_filter(arrays[i])]
        else:
            mut4xh4_d_combinations2arrays_filtered[combination].append(
                fi.mean_filter(arrays[i])
            )


# %%
# Merge the binary arrays
# This is kind of a hack, it doesn't generalize very well

d_combinations2arrays_combined = {}
for combination, arrays in wt_d_combinations2arrays_filtered.items():
    arrays.append(muts23_d_combinations2arrays_filtered[combination][0])
    arrays.append(mut4xh4_d_combinations2arrays_filtered[combination][0])
    d_combinations2arrays_combined[combination] = fi.combine_filters(arrays)


# %%
# Create a dictionary of coordinates
d_combinations2coordinates = {}
for combination, array in d_combinations2arrays_combined.items():
    d_combinations2coordinates[combination] = fi.extract_coordinates(array)


# %%
# Create a dictionary of regions
d_combinations2regions = {}
for combination, coordinates in d_combinations2coordinates.items():
    d_combinations2regions[combination] = fi.extract_regions(coordinates)


# %%
# Create a dictionary of the mean countvalue for each region
wt_d_combinations2means = {}
for combination, regions in d_combinations2regions.items():
    wt_d_combinations2means[combination] = [
        fi.readcounts_to_means(regions, array)
        for array in wt_d_combinations2arrays[combination]
    ]

muts23_d_combinations2means = {}
for combination, regions in d_combinations2regions.items():
    muts23_d_combinations2means[combination] = [
        fi.readcounts_to_means(regions, array)
        for array in muts23_d_combinations2arrays[combination]
    ]

mut4xh4_d_combinations2means = {}
for combination, regions in d_combinations2regions.items():
    mut4xh4_d_combinations2means[combination] = [
        fi.readcounts_to_means(regions, array)
        for array in mut4xh4_d_combinations2arrays[combination]
    ]


# %%
# Format each of the mean dictionaries to a .csv file and save it
for combination, means in wt_d_combinations2means.items():
    iops.format_means_to_table(
        wt_d_combinations2means[combination],
        output_path=f"{RESULT}/wt_{combination}_interactions.csv",
    )

for combination, means in muts23_d_combinations2means.items():
    iops.format_means_to_table(
        muts23_d_combinations2means[combination],
        output_path=f"{RESULT}/muts23_{combination}_interactions.csv",
    )

for combination, means in mut4xh4_d_combinations2means.items():
    iops.format_means_to_table(
        mut4xh4_d_combinations2means[combination],
        output_path=f"{RESULT}/mut4xh4_{combination}_interactions.csv",
    )


# %%
slosh_dataset = {}
for combination in d_combinations2arrays_combined.keys():
        # This loop goest through the same combinations lots of times, this has to be fixed
            slosh_dataset[f"{combination}"] = {
                strain : pd.read_csv(f"{RESULT}/{strain}_{combination}_interactions.csv") for strain in strains
                }


