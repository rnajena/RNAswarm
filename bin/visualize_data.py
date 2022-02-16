import seaborn as sns
import matplotlib.pyplot as plt


def plot_heatmap(array, output_dir, filename):
    heatmap = sns.heatmap(array)
    plt.figure()
    heatmap.figure.savefig(f"{output_dir}/{filename}")
    plt.close("all")


def plot_histogram(array, output_dir):
    plt.hist(array.flatten())
    plt.close()
