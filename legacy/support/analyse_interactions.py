#-------------------------------------------------------------------------------
# Name:        analyse Interactions
# Purpose:     plot alignments 5' and 3' reads to genome
#
# Author:      RS
#
# Created:     02/10/2019
# Updated:     19/02/2020 updated to work with simple interactions file
#-------------------------------------------------------------------------------

import sys, os, re, argparse
from itertools import islice
from time import time
from itertools import zip_longest
from nxviz.plots import CircosPlot
from subprocess import call

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import math
import scipy.stats as stats
import plotly.express as px

def ParseArg():
    p=argparse.ArgumentParser(description = 'Align to genome', epilog = 'Library dependency: Bio, itertools')
    p.add_argument('interactions',type=str,help='interaction file')
    p.add_argument('-cs',type=str, nargs='+',help='chunk size')
    p.add_argument('-t',type=str, nargs='+',help='threshold for minimum reads')
    p.add_argument('-segments',type=str, nargs='+',help='segment name:segment size')
    p.add_argument('-plots',type=str, nargs='+',help='plots: im, ic, ng, circ')

    if len(sys.argv)==1:
        print (p.print_help())
        exit(0)

    return p.parse_args()

def enumerate_double_loop(start1, end1, start2, end2, chunk_size):

    "Produce pairs of indexes in range(n)"
    for i in range(start1, end1, chunk_size):
        for j in range(start2, end2, chunk_size):
            yield i, j

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)

def loadInteractionList(interaction_file, segments):
    "Loads an interaction file into a list structure"
    segment_interactions = []

    for interaction in interaction_file:
        segment1_name, segment1_start, segment1_end, segment2_name, segment2_start, segment2_end = interaction.split("\t")
        try:
            segment1_name = [i for i in segments if "_"+i in segment1_name][0]
            segment2_name = [i for i in segments if "_"+i in segment2_name][0]
        except:
            pass #no match
        else:
            segment1_start = int(segment1_start)
            segment1_end = int(segment1_end)
            segment1_mid = int(segment1_start+(segment1_end-segment1_start)/2)
            segment1 = {"name":segment1_name, "start":segment1_start, "end":segment1_end, "mid":segment1_mid}
            segment2_start = int(segment2_start)
            segment2_end = int(segment2_end)
            segment2_mid = int(segment2_start+(segment2_end-segment2_start)/2)
            segment2 = {"name":segment2_name, "start":segment2_start, "end":segment2_end, "mid":segment2_mid}

            segment_interactions.append(sorted([segment1, segment2], key=lambda segment: segment["name"]))

    return segment_interactions

def convertInteractionListToInteractionArray(interaction_list, segments, segment_lengths, chunk_size=1):
    """Converts a list of interactions into a list of numpy arrays

    chunkSize parameter determines the resolution of analysis :
    that is, a chunk size of 100, means that all fragments that fall into 100 nt windows
    will be combined

    segment names and lengths must be passed as a list to enable chunking
    """

    # Initialise interaction array
    list_of_interaction_arrays = []
    chunks=[math.ceil(x/chunk_size) for x in segment_lengths]

    for i, segment1 in enumerate(segments):
        row = []
        for j, segment2 in enumerate(segments):
            row.append(np.zeros((chunks[i], chunks[j]))) # a numpy array sized for the segment - segment interaction
        list_of_interaction_arrays.append(row)

    # Fill interaction array
    for interaction in interaction_list: # go through each segment interaction in the list
        #extract name, start and end position for each fragment
        seg1 = segments.index(interaction[0]['name'])
        start1 = interaction[0]['start']
        end1 = interaction[0]['end']
        seg2 = segments.index(interaction[1]['name'])
        start2 = interaction[1]['start']
        end2 = interaction[1]['end']

        interaction_array = list_of_interaction_arrays[seg1][seg2] # try to pull out the interaction array corresponding to the two segments

        for i in range(start1, end1, chunk_size): # assign the read to position in the array (increment position), taking into account chunk-size.
            for j in range(start2, end2, chunk_size):
                x = math.ceil(i/chunk_size)-1
                y = math.ceil(j/chunk_size)-1
                interaction_array[x,y]+=1

        list_of_interaction_arrays[seg1][seg2]=interaction_array # put incremented array back into the list

    return list_of_interaction_arrays

def applyThresholdToInteractionArray(interaction_array, threshold):
    interaction_array[interaction_array<threshold]=0

    return interaction_array

def applyThresholdToInteractionList(interaction_list, list_of_interaction_arrays, segments, chunk_size, threshold):
    """Make a new list of interactions, containing only interactions at positions that meet a defined threshold

    interaction_list is a list of the interactions
    list_of_interaction_arrays is a list of arrays containing interaction counts for a given chunk size
    threshold is the minimum number counts at that position for a read to be taken"""

    interaction_list_with_threshold = []

    for interaction in interaction_list: # go through each segment interaction in the list
        #extract name, start and end position for each fragment
        seg1 = segments.index(interaction[0]['name'])
        start1 = interaction[0]['start']
        end1 = interaction[0]['end']
        seg2 = segments.index(interaction[1]['name'])
        start2 = interaction[1]['start']
        end2 = interaction[1]['end']

        interaction_array = list_of_interaction_arrays[seg1][seg2]

        for i, j in enumerate_double_loop(start1, end1, start2, end2, chunk_size):
            x = math.ceil(i/chunk_size)-1
            y = math.ceil(j/chunk_size)-1

            if interaction_array[x,y]>=threshold:
                interaction_list_with_threshold.append(interaction)
                break

    return interaction_list_with_threshold

def countInterSegmentInteractions(list_of_interactions, segments):
    array_size = len(segments)
    interaction_counts = np.zeros((array_size, array_size))

    for interaction in list_of_interactions:
        seg1 = segments.index(interaction[0]['name'])
        seg2 = segments.index(interaction[1]['name'])
        interaction_counts[seg1, seg2] += 1

    # the code below places all intersegment interaction counts on the diagonal
    interaction_counts_transposed = interaction_counts.transpose()

    for i in range(0, array_size):
        for j in range(i+1, array_size):
            interaction_counts[i, j]=interaction_counts[i, j]+interaction_counts_transposed[i,j]

    for i in range(array_size-1, -1, -1):
        for j in range(i-1, -1, -1):
            interaction_counts[i, j]=0

    return interaction_counts

def saveGenomicLinkFile(interactions, filename):
    df_interactions_file = open(filename, "w")

    for interaction in interactions:

        segment1_name = interaction[0]['name']
        segment1_start = interaction[0]['start']
        segment1_end = interaction[0]['end']

        segment2_name = interaction[1]['name']
        segment2_start = interaction[1]['start']
        segment2_end = interaction[1]['end']

        print(segment1_name,segment1_start, segment1_end, segment2_name, segment2_start, segment2_end, file=df_interactions_file)

    df_interactions_file.close()

def saveHotSpotLinkFile(hotspots, filename, chunk_size=100):
    segments = ["HA", "M", "NA", "NP", "NS", "PA", "PB1", "PB2", "18S", "5S", "28S"]
    segments_length = [1775, 1027, 1413, 1565, 890, 2233, 2341, 2341, 1823, 157, 4412]
    chunks=[math.ceil(x/chunk_size) for x in segments_length]

    df_hotspot_file = open(filename, "w")

    for i, segment1 in enumerate(segments):
        for j, segment2 in enumerate(segments):
            segment_specific_hotspot = hotspots[i][j]

            segment1_name = segments[i]
            segment1_start = 1
            segment1_end = segments_length[i]

            segment2_name = segments[j]
            segment2_start = 1
            segment2_end = segments_length[j]

            for k in range(segment1_start, segment1_end, chunk_size):
                for l in range(segment2_start, segment2_end, chunk_size):
                    count = segment_specific_hotspot[math.ceil(k/chunk_size)-1,math.ceil(l/chunk_size)-1]

                    if k+chunk_size < segments_length[i]:
                        pos1_end = k+chunk_size-1
                    else:
                        pos1_end = segments_length[i]

                    if l+chunk_size < segments_length[j]:
                        pos2_end = l+chunk_size-1
                    else:
                        pos2_end = segments_length[j]

                    if count:
                        print(segment1_name, k, pos1_end, segment2_name, l, pos2_end, count, file=df_hotspot_file)

    df_hotspot_file.close()

def loadDistanceMap(distanceMapFile, chunk_size=100, avg=True):
    # Load a distance map file. This is a simple 2d matrix

    previous_pos1 = 0
    distanceMatrix = []
    row = []

    for line in distanceMapFile:
        pos1, pos2, distance = line.split(" ")

        if int(pos1) == int(previous_pos1):
            row.append(distance)

        else: #new row
            distanceMatrix.append(row)
            row = []

        previous_pos1 = pos1

    dimension_size = len(distanceMatrix)
    chunkedDistanceMatrix = []
    row = []

    for i in range(0, dimension_size, chunk_size):
        for j in range(0, dimension_size, chunk_size):
            count = 0
            total = 0
            min = 9999

            for k in range(j, j+chunk_size):
                    try:
                        if avg:
                            total+=float(distanceMatrix[i][k])
                            count+=1
                        else:
                            if float(distanceMatrix[i][k])<min:
                                min=float(distanceMatrix[i][k])
                    except:
                        pass # for when entry is blank

            try:
                avg = total/count
            except:
                avg=0

            if avg:
                row.append(avg)
            else:
                if min==9999:
                    row.append(0)
                else:
                    row.append(min)

        chunkedDistanceMatrix.append(row)
        row = []

    return chunkedDistanceMatrix

def PlotInteractionMatrix(interaction_array, condition, segment1, segment2, chunk_size, threshold):
    plt.figure(figsize = (200,200))
    im = plt.imshow(interaction_array)

#    try: # sometimes error in calculating x axis labels for small segments when chunk size is too big
    x_dimension = np.size(interaction_array, 1)
    y_dimension = np.size(interaction_array, 0)

    xticks = np.arange(0, x_dimension, 25)
    yticks = np.arange(0, y_dimension, 25)

    ax = plt.gca()
    ax.set_title(condition+"_"+segment1+"/"+segment2+": chunk size "+str(chunk_size)+" threshold: "+str(threshold))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, rotation='vertical')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)

    cbar = ax.figure.colorbar(im)
    cbar.ax.set_ylabel(condition+"_"+segment1+"/"+segment2, rotation=-90, va="bottom")
#    except:
#        print("error in plotting:", condition, segment1, segment2, chunk_size, threshold)

    plt.savefig(open("../results/figures/"+condition+"_"+segment1+"-"+segment2+"_"+str(chunk_size)+"_"+str(threshold)+".pdf", 'wb'), format='pdf')
    plt.close()

def PlotHtmlInteractionMatrix(interaction_array, condition, segment1, segment2, chunk_size, threshold):
    fig = px.density_heatmap(interaction_array)
    fig.write_html(open("../results/figures/"+condition+"_"+segment1+"-"+segment2+"_"+str(chunk_size)+"_"+str(threshold)+".html", 'wb'))

def PlotInterSegmentInteractionMatrix(interaction_counts, segments, condition, prefix):
    fig, ax = plt.subplots()
    im = ax.imshow(interaction_counts, aspect='equal')

    cbar = ax.figure.colorbar(im)
    cbar.ax.set_ylabel("interaction strength: "+condition+" "+prefix, rotation=-90, va="bottom")

    # ... and label them with the respective list entries
    ax.set_xticks(np.arange(interaction_counts.shape[1]))
    ax.set_yticks(np.arange(interaction_counts.shape[0]))
    ax.set_xticklabels(segments)
    ax.set_yticklabels(segments)
    #
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    ax.set_xticks(np.arange(interaction_counts.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(interaction_counts.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    ax.set_title("Intersegment Interaction Matrix")
    fig.tight_layout()
    segment_list = ""
    for segment in segments:
        segment_list+=(segment+" ")
    plt.savefig(open("../results/figures/"+condition+" intersegment interaction matrix:"+segment_list+prefix+".pdf", 'wb'), format='pdf')
    plt.close()

def PlotInteractionGraph(interactions, segments, condition, chunk_size, threshold, layout="sfdp"):
    "Makes an interaction graph visualising the predicted disposition of segments"

    G = nx.Graph()
    G.add_nodes_from(segments)

    # Initialise array
    for node in G.nodes:
        G.nodes[node]['weight']=0

    for interaction in interactions:
        seg1=interaction[0]['name']
        seg2=interaction[1]['name']

        if G.has_edge(seg1, seg2):
            #edge already exists, increase weight by one
            G[seg1][seg2]['weight']+=1
        else:
            #new edge, add with weight 1
            G.add_edge(seg1,seg2, weight=1)

        G.nodes[seg1]['weight']+=1
        G.nodes[seg2]['weight']+=1

    nodes = G.nodes()
    edges = G.edges()
    pos = nx.nx_agraph.graphviz_layout(G, prog=layout)

    weights = [G[u][v]['weight'] for u,v in edges]
    weights = [x / 10 for x in weights]

    nodesizes = [n[1]['weight'] for n in G.nodes.data()]

    nx.draw_networkx_nodes(G, pos, node_size=nodesizes)
    nx.draw_networkx_edges(G, pos, edges=edges, width=weights)
    nx.draw_networkx_labels(G, pos)

    segment_list = ""
    for segment in segments:
        segment_list+=(segment+" ")

    plt.savefig(open("../results/figures/"+condition+" graph plot:"+segment_list+chunk_size+" "+threshold+" "+layout+".pdf", 'wb'), format='pdf')

def PlotCircos():

    G = nx.barbell_graph(m1=10, m2=3)
    for n, d in G.nodes(data=True):
        G.node[n]["class"] = choice(["one", "two", "three"])
    c = CircosPlot(G, node_color="class", node_order="class", node_labels=True)
    c.draw()
    plt.show()

def saveGenomicInteractionsFile(interactions, filename):
    interactions_file = open(filename, "w")

    for interaction in interactions:

        segment1_name = interaction[0]['name']
        segment1_start = interaction[0]['start']
        segment1_end = interaction[0]['end']

        segment2_name = interaction[1]['name']
        segment2_start = interaction[1]['start']
        segment2_end = interaction[1]['end']

        print(segment1_name, segment1_start, segment1_end, segment2_name, segment2_start, segment2_end, file=interactions_file)

    interactions_file.close()


def main():
    args=ParseArg()
    t0=time()
    chunkSize = args.cs
    threshold = args.t
    segment_args = args.segments
    plots = args.plots

    segments = []
    segment_lengths = []

    # segments = ["HA", "M", "NA", "NP", "NS", "PA", "PB1", "PB2", "18S", "5S", "28S"]
    # segment_lengths = [1775, 1027, 1413, 1565, 890, 2233, 2341, 2341, 1823, 157, 4412]
    for segment in segment_args:
        segments.append(segment.split(":")[0])
        segment_lengths.append(int(segment.split(":")[1]))

    interaction_file = open(args.interactions, "r")
    condition = args.interactions.split("/")[-1].split(".")[0] # takes the condition name as the string before the first '.'
    list_of_interactions = loadInteractionList(interaction_file, segments)

    # Draw intra and inter segment interaction map
    if "im" in plots:
        for cs in chunkSize:
            for t in threshold:
                list_of_interaction_arrays = convertInteractionListToInteractionArray(list_of_interactions, segments, segment_lengths, int(cs))

                for segment1 in segments:
                    for segment2 in segments:
                        print("Plotting:", segment1, segment2, cs, t)
                        seg1 = segments.index(segment1)
                        seg2 = segments.index(segment2)

                        interaction_array = list_of_interaction_arrays[seg1][seg2]
                        print(segment1, segment2)

                        interactions = applyThresholdToInteractionArray(interaction_array, int(t))

                        if segment1==segment2: #intra-molecular interactions, can transpose data
                            transposed_interactions = interactions.transpose()
                            merged_interactions = interactions+transposed_interactions
                            #log_interactions = np.where(merged_interactions>0, log(merged_interactions), 0)
                            log_interactions = np.log(merged_interactions)
                            log_interactions[log_interactions==-np.inf]=0
                            PlotInteractionMatrix(interactions, condition+"_raw", segment1, segment2, chunk_size=int(cs), threshold=int(t))
                            PlotInteractionMatrix(merged_interactions, condition+"_merged", segment1, segment2, chunk_size=int(cs), threshold=int(t))
                            PlotInteractionMatrix(log_interactions, condition+"_log", segment1, segment2, chunk_size=int(cs), threshold=int(t))

                            PlotHtmlInteractionMatrix(interactions, condition+"_raw", segment1, segment2, chunk_size=int(cs), threshold=int(t))
                            PlotHtmlInteractionMatrix(merged_interactions, condition+"_merged", segment1, segment2, chunk_size=int(cs), threshold=int(t))
                            PlotHtmlInteractionMatrix(log_interactions, condition+"_log", segment1, segment2, chunk_size=int(cs), threshold=int(t))
                        else:
                            #log_interactions = np.where(log_interactions>0, log(log_interactions), 0)
                            log_interactions = np.log(interactions)
                            log_interactions[log_interactions==-np.inf]=0
                            PlotInteractionMatrix(interactions, condition+"_raw", segment1, segment2, chunk_size=int(cs), threshold=int(t))
                            PlotInteractionMatrix(log_interactions, condition+"_log", segment1, segment2, chunk_size=int(cs), threshold=int(t))

    # Draw intersegment count matrix

    # TODO #
    # should also plot these for varying thresholds and chunk size

    if "ic" in plots:
        intersegment_interaction_counts = countInterSegmentInteractions(list_of_interactions, segments)
        intersegment_interaction_counts_log = np.log(intersegment_interaction_counts)
        PlotInterSegmentInteractionMatrix(intersegment_interaction_counts, segments, condition, "raw")
        PlotInterSegmentInteractionMatrix(intersegment_interaction_counts_log, segments, condition, "log")

    # Draw network graph
    if "ng" in plots:
        for cs in chunkSize:
            for t in threshold:
                list_of_interaction_arrays = convertInteractionListToInteractionArray(list_of_interactions, segments, segment_lengths, int(cs))
                list_of_interactions_with_threshold = applyThresholdToInteractionList(list_of_interactions, list_of_interaction_arrays, segments, int(cs), int(t))
                PlotInteractionGraph(list_of_interactions_with_threshold, segments, condition, cs, t)

    # Draw circos plots
    if "circ" in plots:
        for cs in chunkSize:
            for t in threshold:
                list_of_interaction_arrays = convertInteractionListToInteractionArray(list_of_interactions, segments, segment_lengths, int(cs))
                list_of_interaction_arrays_with_threshold = [applyThresholdToInteractionArray(ia2, int(t)) for ia1 in list_of_interaction_arrays for ia2 in ia1]
                list_of_interactions_with_threshold = applyThresholdToInteractionList(list_of_interactions, list_of_interaction_arrays, segments, int(cs), int(t))

                # plot interaction intensities at each position (excluding interactions that don't meet threshold) with list_of_interaction_arrays_with_threshold
                saveGenomicInteractionsFile(list_of_interactions_with_threshold, "../results/figures/temp.df")
                RscriptArgs = ["Rscript", "5_circosPlot.R", "../results/figures/temp.df", "-segments"]+segment_args+["-plots", "circ", "-o", "../results/figures/"+condition+"_circos_"+cs+"_"+t+".pdf"]
                print(' '.join(map(str, RscriptArgs)))
                call(RscriptArgs)

                # plot interaction intensities at each position (excluding interactions that don't meet threshold) with list_of_interaction_arrays_with_threshold

                # plot actual interactions that meet threshold with list_of_interactions_with_threshold



    # ribosome_interaction_file = open(args.ribosome_interactions, "r")
    # ribosome_interactions = loadDistanceMap(ribosome_interaction_file,chunk_size=cs, avg=False)
    # max = np.amax(ribosome_interactions)
    # ribosome_interactions = np.abs(ribosome_interactions-max)
    #PlotIt(ribosome_interactions, "Ribosome distances (min)", chunk_size=cs, threshold="NA")

    #Stats
    # ribosome_interactions[ribosome_interactions==max]=0
    # s = stats.wilcoxon(np.ndarray.flatten(ribosome_interactions), np.ndarray.flatten(hotspot_log), zero_method='wilcox')
    # print(s)

    # Genome links
    # hotspots, hotspot_interactions = findHotspot(list_of_interactions, threshold=200)
    # saveGenomicLinkFile(hotspot_interactions, 'PR8_hotspot_interactions.df')
    # saveHotSpotLinkFile(hotspots, 'PR8_hotspots.df', chunk_size=100)

if __name__ == '__main__':
    main()
