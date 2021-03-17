import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Import the array for HA/HA interactions in the WT sample R1_wt
wt_r1_ha_ha = np.load('data/reads/R1_wt_npy/splash_q20_noN_min6_R1_1120_wtI_cs1_HA_HA_inf.npy')

def degenerate_array(array):
    """
    (ndarray) -> ndarray

    """
    return array + np.random.randint(10, size=array.shape)

ax = sns.heatmap(wt_r1_ha_ha)
ax_figure = ax.get_figure()
ax_figure.savefig('results/figures/output.png')



