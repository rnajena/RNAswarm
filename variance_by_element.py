import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def degenerate_array(array):
    """ (ndarray) -> ndarray

    """
    return array + np.random.randint(10, size=array.shape)

# Import the array for HA/HA interactions in the WT sample R1_wt
wt_r1_ha_ha = np.load('data/reads/R1_wt_npy/splash_q20_noN_min6_R1_1120_wtI_cs1_HA_HA_inf.npy')

# "Degenerate" arrays for further variance calculation
wt_r1_ha_ha_dg01 = degenerate_array(wt_r1_ha_ha)
wt_r1_ha_ha_dg02 = degenerate_array(wt_r1_ha_ha)
wt_r1_ha_ha_dg03 = degenerate_array(wt_r1_ha_ha)
wt_r1_ha_ha_dg04 = degenerate_array(wt_r1_ha_ha)
wt_r1_ha_ha_dg05 = degenerate_array(wt_r1_ha_ha)

# Calculate variance by element and create variance array
wt_r1_ha_ha_variance = (np.var([wt_r1_ha_ha, wt_r1_ha_ha_dg01,
                        wt_r1_ha_ha_dg02, wt_r1_ha_ha_dg03, 
                        wt_r1_ha_ha_dg04, wt_r1_ha_ha_dg05], axis=0))

# Calculate the histogram of the variance array
wt_r1_ha_ha_histogram = np.histogram(wt_r1_ha_ha_variance, bins=100)
print(wt_r1_ha_ha_histogram)

# Tests on histogram plotting
histogram_plot = sns.displot(wt_r1_ha_ha_histogram[0])
histogram_plot.savefig('results/figures/histogram_plot_test.png')

# Tests on heatmap plotting
heatmap_plot = sns.heatmap(wt_r1_ha_ha_variance)
heatmap_plot.figure.savefig('results/figures/heatmap_plot_test.png')

