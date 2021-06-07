# Import functions for analysis of SPLASH np.arrays
import variance_by_element as vbe

# Import functions from the SciPy stack
import numpy as np
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
import seaborn as sns

# Define file paths and viral features:
DIRECTORY = "/data/dessertlocal/projects/gl_iav-splash_freiburg"
INPUT = f"{DIRECTORY}/data/arrays"
RESULT = f"{DIRECTORY}/results/202104/20210607"
iav_segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]


# Read the arrays from file and put them into dictionaries:
wt_d_repDir2Combinations, wt_d_combination2Array = vbe.read_arrays(INPUT, iav_segments)

# Visualising arrays:
# We can visualize one of our arrays, in this case the array representing the # interaction between NA and NP segments for sample 1
sample01_NA_NP = wt_d_combination2Array["NA_NP"][1]

# We can also generate a heatmap of this array:
vbe.plot_heatmap(sample01_NA_NP, RESULT)

# Finding local maxima:
# In order to find a local maxima we need to define a treshold:
treshold = 100

# Plotting variance versus read-count:
