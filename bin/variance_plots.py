# %%
from statistics import mean
import numpy as np
import matplotlib.pyplot as plt
import os
import itertools
import sklearn.linear_model as linear_model
import sklearn.metrics as metrics

# %%
def import_arrays(dir_name):
    """
    Import all arrays in a directory

    input:
        dir_name: directory name

    output:
        array_list: dictionarry of arrays
    """
    files = os.listdir(dir_name)
    arrays = {}
    for file in files:
        arrays[file] = np.load(dir_name + '/' + file).flatten()
    return arrays

def plot_heatmaps(arrays, dir_name):
    """
    Plot heatmaps for all arrays and output to directory

    input:
        arrays: list of arrays
        dir_name: directory name

    output:
        None
    """
    for array in arrays:
        plt.hist(array, bins=100)
        plt.savefig(dir_name + '/' + str(array.shape[0]) + '_' + str(array.shape[1]) + '.png')

def calculate_variance_array(arrays):
    """
    Calculate variance array for all arrays

    input:
        arrays: dict of arrays

    output:
        variance_array: variance array
    """
    variance_array = np.var([arrays[0], arrays[1], arrays[2]], axis=0)
    return variance_array


def calculate_mean_array(arrays):
    """
    Calculate mean array for all arrays

    input:
        arrays: dict of arrays

    output:
        mean_array: mean array
    """
    mean_array = np.mean([arrays[0], arrays[1], arrays[2]], axis=0)
    return mean_array

def plot_variance_vs_mean(mean_array, variance_array, model, filename):
    """
    Plot variance vs mean for all arrays, plot a  and output to directory

    input:
        mean_array: mean array
        variance_array: variance array
        linear_model: linear model
        filename: filename

    output:
        None
    """
    plt.scatter(mean_array, variance_array, alpha=0.1)
    plt.plot(mean_array, model.predict(mean_array), color='k')
    plt.xlabel('Mean')
    plt.ylabel('Variance')
    plt.savefig(f'{filename}.png', dpi=300)
    plt.close()


# %%
dir_4xHA = '/home/ru27wav/Projects/gl_iav-splash_freiburg/data/arrays/4xHA'
dir_S23 = '/home/ru27wav/Projects/gl_iav-splash_freiburg/data/arrays/S23'
dir_wt0120 = '/home/ru27wav/Projects/gl_iav-splash_freiburg/data/arrays/wt0120'
dir_wt1120_I = '/home/ru27wav/Projects/gl_iav-splash_freiburg/data/arrays/wt1120_I'
dir_wt1120_II = '/home/ru27wav/Projects/gl_iav-splash_freiburg/data/arrays/wt1120_II'

arrays_4xHA = import_arrays(dir_4xHA)
arrays_S23 = import_arrays(dir_S23)
arrays_wt0120 = import_arrays(dir_wt0120)
arrays_wt1120_I = import_arrays(dir_wt1120_I)
arrays_wt1120_II = import_arrays(dir_wt1120_II)

segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
segment_combinations = [segment_combination for segment_combination in itertools.combinations_with_replacement(segments, 2)]

# %%
dict_replicates = {}
dict_replicates['merged'] = {}
for segment_combination in segment_combinations:
    dict_replicates[segment_combination] = {}
    replicate_iterator = 0
    for replicate in [arrays_wt0120, arrays_wt1120_I, arrays_wt1120_II]:
        for filename in replicate.keys():
            if f'{segment_combination[0]}_{segment_combination[1]}' in filename:
                dict_replicates[segment_combination][replicate_iterator] = replicate[filename]
                if replicate_iterator in dict_replicates['merged']:
                    dict_replicates['merged'][replicate_iterator] = np.concatenate((dict_replicates['merged'][replicate_iterator], replicate[filename]))
                else:
                    dict_replicates['merged'][replicate_iterator] = replicate[filename]
            elif f'{segment_combination[1]}_{segment_combination[0]}' in filename:
                dict_replicates[segment_combination][replicate_iterator] = replicate[filename]
                if replicate_iterator in dict_replicates['merged']:
                    dict_replicates['merged'][replicate_iterator] = np.concatenate((dict_replicates['merged'][replicate_iterator], replicate[filename]))
                else:
                    dict_replicates['merged'][replicate_iterator] = replicate[filename]
        replicate_iterator += 1

# %%
for key in dict_replicates.keys():
    dict_replicates[key]['mean'] = calculate_mean_array(dict_replicates[key])
    dict_replicates[key]['variance'] = calculate_variance_array(dict_replicates[key])

# %%
for key in dict_replicates.keys():
    model = linear_model.LinearRegression()
    model.fit(dict_replicates[key]['mean'].reshape(-1, 1), dict_replicates[key]['variance'].reshape(-1, 1))
    dict_replicates[key]['linear_regression'] = model
    dict_replicates[key]['r_squared'] = metrics.r2_score(dict_replicates[key]['variance'].reshape(-1, 1), model.predict(dict_replicates[key]['mean'].reshape(-1, 1)))
    dict_replicates[key]['mean_squared_error'] = metrics.mean_squared_error(dict_replicates[key]['variance'].reshape(-1, 1), model.predict(dict_replicates[key]['mean'].reshape(-1, 1)))

# %%
print('combination,coeficients,intercept,r_squared,mean_squared_error')
for key in dict_replicates.keys():
    print(f'{key},{dict_replicates[key]["linear_regression"].coef_},{dict_replicates[key]["linear_regression"].intercept_},{dict_replicates[key]["r_squared"]},{dict_replicates[key]["mean_squared_error"]}')

# %%
for key in dict_replicates.keys():
    plot_variance_vs_mean(
        dict_replicates[key]['mean'].reshape(-1, 1),
        dict_replicates[key]['variance'].reshape(-1, 1),
        dict_replicates[key]['linear_regression'],
        f'/home/ru27wav/Projects/gl_iav-splash_freiburg/results/variance_plots/variance_plots/{str(key[0])}_{str(key[1])}'
        )
