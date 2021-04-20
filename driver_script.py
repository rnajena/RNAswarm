# %%
# Import functions for analysis of SPLASH np.arrays
import variance_by_element


# %%
# Define file paths and viral features
DIRECTORY = '/data/dessertlocal/projects/gl_iav-splash_freiburg'  
INPUT = f'{DIRECTORY}/data/arrays'
RESULT = f'{DIRECTORY}/results/202104/20210420'
iav_segments = ['PB2','PB1','PA','HA','NP','NA','M','NS']


# %%
# 
wt_d_repDir2Combinations, wt_d_combination2Array = variance_by_element.read_arrays(INPUT, iav_segments)


# %%
# 
wt_d_comb2variances = variance_by_element.calculate_variances(wt_d_combination2Array)


# %%
# 
variance_by_element.save_heatmaps(wt_d_comb2variances, RESULT)


# %%



