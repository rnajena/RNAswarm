# Import functions for analysis of SPLASH np.arrays
import io_ops as iops
import find_interactions as fi

# Define file paths and viral features:
DIRECTORY = "/home/ru27wav/Projects/gl_iav-splash_freiburg/"
INPUT = f"{DIRECTORY}/data/arrays_test"
RESULT = f"{DIRECTORY}/results/test"
iav_segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]

# Read the arrays from file and put them into dictionaries
wt_d_repDir2Combinations, wt_d_combinations2arrays = iops.read_arrays(
    INPUT, iav_segments
)

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

# Merge the binary arrays
wt_d_combinations2arrays_combined = {}
for combination, arrays in wt_d_combinations2arrays_filtered.items():
    wt_d_combinations2arrays_combined[combination] = fi.combine_filters(arrays)

# Create a dictionary of coordinates
wt_d_combinations2coordinates = {}
for combination, array in wt_d_combinations2arrays_combined.items():
    wt_d_combinations2coordinates[combination] = fi.extract_coordinates(array)

# Create a dictionary of regions
wt_d_combinations2regions = {}
for combination, coordinates in wt_d_combinations2coordinates.items():
    wt_d_combinations2regions[combination] = fi.extract_regions(coordinates)

# Create a dictionary of the mean countvalue for each region
wt_d_combinations2means = {}
for combination, regions in wt_d_combinations2regions.items():
    wt_d_combinations2means[combination] = [
        fi.readcounts_to_means(regions, array)
        for array in wt_d_combinations2arrays[combination]
    ]

# Format each of the mean dictionaries to a .csv file and save it
for combination, means in wt_d_combinations2means.items():
    fi.format_means_to_table(
        wt_d_combinations2means[combination],
        output_path=f"{RESULT}/{combination}_interactions.csv",
    )
