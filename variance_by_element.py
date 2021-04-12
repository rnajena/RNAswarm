import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy import stats


def plot_heatmap(array, output_dir):
    heatmap = sns.heatmap(array)
    plt.figure()
    heatmap.figure.savefig(output_dir)


def plot_histogram(array, output_dir):
    histogram = plt.figure()
    ax = histogram.add_subplot(111)
    ax.hist(array.flatten(), bins='auto')
    histogram.tight_layout()
    histogram.savefig(output_dir)
    plt.close(histogram)


def read_arrays(data_directory):
    """
    Takes a directory containing subdirectories related to diferent replicates and returns
    two dictionaries, one describing the paths to the arrays and another containing the 
    arrays for each segment combination.

    Input:
        
    Return:
        
    """
    d_repDir2Combinations = {}
    d_combination2Array = {}

    replicateDirs = os.listdir(data_directory)
    for repDir in replicateDirs:
        allCombinations = os.listdir(repDir)
        # dict = {'replicate01_path' : [replicate01_arrays]}
        d_repDir2Combinations[repDir] = allCombinations
        

        for combination in allCombinations:
            # /data/dessertlocal/projects/gl_iav-splash_freiburg/data/arrays/wt0120/ --> 
            # wt0120 (basename)
            uniqueID = f'{os.basename(repDir)}/{combination[39:-8]}'
            array = np.load(f'{repDir}/{combination}')
            # dict = {'Seg1-Seg2' : [np.array(array01), np.array(array02), 
            #                        np.array(array03)]}
            if d_combination2Array[uniqueID]:
                d_combination2Array[uniqueID].extend(array)
            else:
                d_combination2Array[uniqueID] = [array]
            

    return (d_repDir2Combinations, d_combination2Array)


def calculate_variances(d_arrays):
    """
    Takes a dictionary with numpy arrays and calculates the variance for each unique entry.

    Input:
        d_arrays -- {'uniqueName' : [np.arrays]}

    Return:
        d_comb2variance -- {'uniqueName' : np.array}
    
    """
    d_comb2variance = {}
    for combination, l_countTable in d_arrays.items():
        d_comb2variance[combination] = np.var(l_countTable, axis=0)
    
    return d_comb2variance


def check_heatmaps(variances, output_dir):
    """
    Takes a dictionary of np.arrays and plots heatmaps for the np.arrays and saves them
    to the output_dir

    Input:

    Return:
    """
    # Retrieves the key (comb) and the value (variance_array)
    for comb, variance_array in variances.items():
        plot_heatmap(variance_array, f'{output_dir}/{comb}_hist')


def check_normality(variances):
    normality_summary = {}

    for comb, variance_array in variances.items():
        k2, p = stats.normaltest(variance_array.flatten())
        if p < 0.05:
            normality_summary[comb] = 'normal'       
        else:
            normality_summary[comb] = 'non-normal'
    return normality_summary  


def check_histograms(variances, output_dir):
    for comb, variance_array in variances.items():
        plot_histogram(variance_array, f'{output_dir}/{comb}_hist')


DIRECTORY = '/data/dessertlocal/projects/gl_iav-splash_freiburg/'  
INPUT = f"{DIRECTORY}/data/arrays/"
RESULT = f"{DIRECTORY}/results/202104/20200412"

# 
wt_d_repDir2Combinations, wt_d_combination2Array = read_arrays(INPUT)
# 
wt_d_comb2variances = calculate_variances(wt_d_combination2Array)
# 
check_heatmaps(wt_d_comb2variances, RESULT)
