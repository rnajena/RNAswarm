import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import configparser

def degenerate_array(array):
    """ (ndarray) -> ndarray
    Returns a 'degenerate' copy of an ndarray by adding random int from 0 to 9
    to each of the elements in the array.
    """
    return array + np.random.randint(10, size=array.shape)

def get_config(configuration_file):
    """
    """
    if configuration_file is str:
        config = configparser.ConfigParser()
        config.read(configuration_file)


def get_filenames():
    os.listdir()

def calculate_variance_array(array01, array02, array03):
    """
    """
    return np.var(array01, array02, array03)

def plot_heatmap(array):
    heatmap_plot = sns.heatmap(wt_r1_ha_ha_variance)
    heatmap_plot.figure.savefig('results/figures/heatmap_plot_test.png')

def plot_histogram(array):
    np.histogram(wt_r1_ha_ha_variance, bins=100)
    histogram_plot = sns.displot(wt_r1_ha_ha_histogram[0])
    histogram_plot.savefig('results/figures/histogram_plot_test.png')

def analyse_experiment():
    print("Hey")

