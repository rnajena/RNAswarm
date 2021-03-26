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

def check_heatmaps(replicate01_path, replicate02_path, replicate03_path, output_dir):
    replicate01_arrays = os.listdir(replicate01_path)
    replicate02_arrays = os.listdir(replicate02_path)
    replicate03_arrays = os.listdir(replicate03_path)
    
    replicate01_arrays.sort()
    replicate02_arrays.sort()
    replicate03_arrays.sort()

    for i in range(36):
        array01 = np.load(replicate01_path + replicate01_arrays[i])
        array02 = np.load(replicate02_path + replicate02_arrays[i])
        array03 = np.load(replicate03_path + replicate03_arrays[i])
        variance_array = np.var([array01, array02, array03], axis=0)
        plot_heatmap(variance_array, (output_dir + '/' + replicate01_arrays[i][39:-8]))

def check_normality(replicate01_path, replicate02_path, replicate03_path):
    replicate01_arrays = os.listdir(replicate01_path)
    replicate02_arrays = os.listdir(replicate02_path)
    replicate03_arrays = os.listdir(replicate03_path)
    
    replicate01_arrays.sort()
    replicate02_arrays.sort()
    replicate03_arrays.sort()
    
    normality_summary = {}

    for i in range(36):
        array01 = np.load(replicate01_path + replicate01_arrays[i])
        array02 = np.load(replicate02_path + replicate02_arrays[i])
        array03 = np.load(replicate03_path + replicate03_arrays[i])
        variance_array = np.var([array01, array02, array03], axis=0)
        k2, p = stats.normaltest(variance_array.flatten())
        if p < 0.05:
            normality_summary[replicate01_arrays[i][39:-8]] = 'normal'       
        else:
            normality_summary[replicate01_arrays[i][39:-8]] = 'non-normal'
    return normality_summary  

def check_histograms(replicate01_path, replicate02_path, replicate03_path, output_dir):
    replicate01_arrays = os.listdir(replicate01_path)
    replicate02_arrays = os.listdir(replicate02_path)
    replicate03_arrays = os.listdir(replicate03_path)
    
    replicate01_arrays.sort()
    replicate02_arrays.sort()
    replicate03_arrays.sort()

    for i in range(36):
        array01 = np.load(replicate01_path + replicate01_arrays[i])
        array02 = np.load(replicate02_path + replicate02_arrays[i])
        array03 = np.load(replicate03_path + replicate03_arrays[i])
        variance_array = np.var([array01, array02, array03], axis=0)
        plot_histogram(variance_array, (output_dir + '/' + replicate01_arrays[i][39:-8] + '_hist'))

wt_histograms = check_histograms('/data/dessertlocal/projects/gl_iav-splash_freiburg/data/arrays/wt0120/',
                                 '/data/dessertlocal/projects/gl_iav-splash_freiburg/data/arrays/wt1120_I/',
                                 '/data/dessertlocal/projects/gl_iav-splash_freiburg/data/arrays/wt1120_II/',
                                 '/data/dessertlocal/projects/gl_iav-splash_freiburg/results/202003/20200325')

wt_heatmaps = check_heatmaps('/data/dessertlocal/projects/gl_iav-splash_freiburg/data/arrays/wt0120/',
                             '/data/dessertlocal/projects/gl_iav-splash_freiburg/data/arrays/wt1120_I/',
                             '/data/dessertlocal/projects/gl_iav-splash_freiburg/data/arrays/wt1120_II/',
                             '/data/dessertlocal/projects/gl_iav-splash_freiburg/results/202003/20200325')

wt_normality = check_normality('/data/dessertlocal/projects/gl_iav-splash_freiburg/data/arrays/wt0120/',
                               '/data/dessertlocal/projects/gl_iav-splash_freiburg/data/arrays/wt1120_I/',
                               '/data/dessertlocal/projects/gl_iav-splash_freiburg/data/arrays/wt1120_II/')

print(wt_normality)
