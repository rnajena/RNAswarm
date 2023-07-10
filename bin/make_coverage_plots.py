#!/usr/bin/env python3

"""make_coverage_plots.py

Takes a fasta file and a trns.txt file from segemehl and makes coverage plots for each segment

Usage:
    make_coverage_plots.py <fasta> <trns> -o <output_folder> [--DMS_data <DMS_data>]
    make_coverage_plots.py <fasta> <trns> -o <output_folder> [--DMS_data <DMS_data>] [--interactive]
    make_coverage_plots.py -h | --help

Options:
    -h --help                       Show this screen.
    <fasta>                         Fasta file
    <trns>                          Trns file
    -o --output=<output_folder>     Output folder
    --DMS_data=<DMS_data>           DMS data sheet
    --interactive                   Create interactive plots
                                    using plotly
"""

from docopt import docopt
import helper as hp
import trns_handler as th
import numpy as np
import os
import matplotlib.pyplot as plt
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def make_coverage_dict(trnsFile, genome_dict):
    """
    Makes a dictionary of coverage arrays
    
    Parameters
    ----------
    trns : str
        Path to trns.txt file
    genome_dict : dict
        Dictionary of genome sequences

    Returns
    -------
    cov_dict : dict
        Dictionary of coverage arrays
    """
    # Generate empty coverage dictionary
    cov_dict = {}
    for segment, sequence in genome_dict.items():
        cov_dict[segment] = np.zeros(len(sequence))
    # Fill coverage dictionary
    with open(trnsFile) as inputStream:
        for line in inputStream:
            line = line.strip().split()
            firstRead = line[0].split(",")
            secondRead = line[1].split(",")
            currentRow = th.__extract_start_stop_segemehl(
                firstRead
            ) + th.__extract_start_stop_segemehl(secondRead)
            interaction = th.__check_interaction(currentRow)
            cov_dict[interaction[0]][interaction[1]:interaction[2]] += 1
            cov_dict[interaction[3]][interaction[4]:interaction[5]] += 1
    return cov_dict


def parse_DMS_data(DMS_data):
    """
    Parses DMS data sheet
    
    Parameters
    ----------
    DMS_data : str
        Path to DMS data sheet

    Returns
    -------
    DMS_dict : dict
        Dictionary of DMS data
    """
    DMS_dict = {}
    with open(DMS_data) as inputStream:
        for line in inputStream:
            line = line.strip().split(",")
            # line[0] = segment name
            # line[1] = DMS reactivity file
            # Now we parse the DMS reactivity file as a numpy array
            # each line is the reactivity at a given position (there is only one column)
            DMS_dict[line[0]] = np.loadtxt(line[1])
            # If the value is equal to -1 we set it to NaN
            DMS_dict[line[0]][DMS_dict[line[0]] == -1] = np.nan
    return DMS_dict


def plot_coverage(cov_dict, output_folder, DMS_dict=None):
    """
    Plots coverage for each segment
    
    Parameters
    ----------
    cov_dict : dict
        Dictionary of coverage arrays
    output_folder : str
        Path to output folder
    DMS_dict : dict
        Dictionary of DMS data

    Returns
    -------
    None
    """
    for segment, coverage in cov_dict.items():
        # Plot coverage
        plt.plot(coverage)
        plt.xlabel("Position")
        plt.ylabel("Coverage")
        plt.title(segment)
        plt.savefig(os.path.join(output_folder, f"{segment}_coverage.png"))
        plt.close()
        # Plot coverage and DMS data
        if DMS_dict:
            # Now plot the two together (with the same y axis)
            fig, ax1 = plt.subplots()
            ax1.plot(coverage, color="blue")
            ax1.set_xlabel("Position")
            ax1.set_ylabel("Coverage", color="blue")
            ax1.tick_params(axis="y", labelcolor="blue")
            ax2 = ax1.twinx()
            ax2.plot(DMS_dict[segment], color="red")
            ax2.set_ylabel("DMS reactivity", color="red")
            ax2.tick_params(axis="y", labelcolor="red")
            plt.title(segment)
            plt.savefig(os.path.join(output_folder, f"{segment}_coverage_DMS.png"))
            plt.close()
            # Plot the same thing again, coverage in the front this time
            fig, ax1 = plt.subplots()
            ax1.plot(DMS_dict[segment], color="red")
            ax1.set_xlabel("Position")
            ax1.set_ylabel("DMS reactivity", color="red")
            ax1.tick_params(axis="y", labelcolor="red")
            ax2 = ax1.twinx()
            ax2.plot(coverage, color="blue")
            ax2.set_ylabel("Coverage", color="blue")
            ax2.tick_params(axis="y", labelcolor="blue")
            plt.title(segment)
            plt.savefig(os.path.join(output_folder, f"{segment}_DMS_coverage.png"))
            plt.close()


def plotly_coverage(cov_dict, genome_dict, output_folder, DMS_dict=None):
    # iterate over each segment on cov_dict and genome_dict at the same time
    for (segment, coverage), (_, sequence) in zip(cov_dict.items(), genome_dict.items()):
        # Create object to hold the figure
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        # Add coverage trace
        fig.add_trace(
            go.Scatter(x=np.arange(len(coverage)), y=coverage, name="Coverage"),
            secondary_y=False,
        )
        # Add DMS trace if DMS_dict is not None
        if DMS_dict:
            fig.add_trace(
                go.Scatter(
                    x=np.arange(len(DMS_dict[segment])),
                    y=DMS_dict[segment],
                    name="DMS reactivity",
                ),
                secondary_y=True,
            )
        # Add labels to the x axes
        fig.update_xaxes(title_text="position")
        # Plot the sequence parallel to the x axis
        fig.add_trace(
            go.Scatter(
                x=np.arange(len(sequence)),
                y=[-1] * len(sequence),
                mode="text",
                text=list(sequence),
                textposition="bottom center",
                # make sure there is no legend for this trace
                showlegend=False,
            ),
            secondary_y=False,
        )
        # Add labels to the y axes
        fig.update_yaxes(title_text="chimeric read coverage", secondary_y=False)
        fig.update_yaxes(title_text="DMS reactivity", secondary_y=True)
        # Limit zooming to the x axis
        fig.update_layout(xaxis=dict(rangeslider=dict(visible=True)))
        # Add title
        fig.update_layout(title_text=segment)
        # Save figure
        fig.write_html(os.path.join(output_folder, f"{segment}_coverage.html"))

def main():
    args = docopt(__doc__)
    fasta = args["<fasta>"]
    trns = args["<trns>"]
    output_folder = args["--output"]
    DMS_data = args["--DMS_data"]

    # Process input files
    genome_dict = hp.parse_fasta(fasta)
    cov_dict = make_coverage_dict(trns, genome_dict)
    if DMS_data:
        DMS_dict = parse_DMS_data(DMS_data)
    else:
        DMS_dict = None

    # Plot coverage
    plot_coverage(cov_dict, output_folder, DMS_dict)

    # Plot coverage with plotly
    plotly_coverage(cov_dict, genome_dict, output_folder, DMS_dict)


if __name__ == "__main__":
    main()