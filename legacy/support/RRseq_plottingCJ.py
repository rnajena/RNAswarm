import sys
#sys.path.append('/home/fassakima/SPLASHanalysis/scripts')
from .helper_v1CJ import enumerate_double_loop, grouper, log

import sys, os, re, argparse
from itertools import islice
from time import time
from itertools import zip_longest
#from nxviz.plots import CircosPlot
from subprocess import call

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import scipy.stats as stats
import plotly.express as px

# loadChimeraFile
# chimeraListToInteractionArray

def loadChimeraFile(chimeraFileName, segments, maxToLoad=np.inf):
    """Loads a chimera file into a list structure
    this chimera file has a simple structure
    segment1_name, segment1_start, segment1_end, segment2_name, segment2_start, segment2_end"""
    
    chimera_list = []
    chimeraFile = open(chimeraFileName, "r")
    
    for i, interaction in enumerate(chimeraFile):
        segment1_name, segment1_start, segment1_end, segment2_name, segment2_start, segment2_end = interaction.rstrip().split("\t")
        
        try:
            segment1_name = [i for i in segments if "_"+i in segment1_name][0] #match _segment in array
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

            chimera_list.append(sorted([segment1, segment2], key=lambda segment: segment["name"]))

        if (i+1)%maxToLoad==0:
            break
    
    chimeraFile.close()
    
    return chimera_list

##########################
# Make interaction array #
##########################

def loadInteractionArray(fileName, wd=""):
    with open(os.path.join(wd, fileName), 'rb') as f:
        return np.load(f)
    
def saveInteractionArrays(interaction_files, segments, segment_lengths, cs, maxToLoad=np.Inf, outdirectory=""):
    for interaction_file in interaction_files:
        log("Loading "+interaction_file)
        list_of_interactions = loadChimeraFile(interaction_file, segments, maxToLoad=maxToLoad)
        list_of_interaction_arrays = chimeraListToInteractionArray(list_of_interactions, segments, segment_lengths, cs)
        
        for segment1 in segments:
            for segment2 in segments:
                # only plot for segments in alphbetical order (to prevent blank plots)
                #if segment1<segment2:
                if segment1<=segment2:
                    seg1 = segments.index(segment1)
                    seg2 = segments.index(segment2)
                    log("Analysing... %s and %s" % (segment1, segment2))

                    interaction_array = list_of_interaction_arrays[seg1][seg2]

                    fileName = os.path.join(outdirectory, os.path.split(interaction_file)[1].split(".chim")[0]+"_cs"+str(cs)+"_"+str(segment1)+"_"+str(segment2)+"_"+str(maxToLoad)+".npy")
                    
                    log("Saving... %s" % fileName)

                    with open(fileName, 'wb') as f:
                        np.save(f, interaction_array)
                    #np.save(fileName, interaction_array)

    log("Finished") 

def checkDebugChimeraFile(chimeraFileName, segments, segment_lengths, maxToLoad=np.inf):
    """Loads a chimera file into a list structure
    this chimera file has a simple structure
    segment1_name, segment1_start, segment1_end, segment2_name, segment2_start, segment2_end"""
    
    chimera_list = []
    chimeraFile = open(chimeraFileName, "r")
    #self.rname, self.pos, self.aend, self.seq, self.cigar, self.mapq
    for i, interaction in enumerate(chimeraFile):
        segment1_name, segment1_start, segment1_end, segment1_seq, segment1_cigar, segment1_mapq, segment2_name, segment2_start, segment2_end, segment2_seq, segment2_cigar, segment2_mapq = interaction.rstrip().split("\t")
        
        try:
            segment1_name = [i for i in segments if "_"+i in segment1_name][0] #match _segment in array
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
            
            if segment1_end > segment_lengths(segments.index(segment1_name)) or segment2_end > segment_lengths(segments.index(segment2_name)):
                print("debug: ", segment1_name, segment1_start, segment1_end, segment1_seq, segment1_cigar, segment1_mapq, segment2_name, segment2_start, segment2_end, segment2_seq, segment2_cigar, segment2_mapq)

        if (i+1)%maxToLoad==0:
            break
    
    chimeraFile.close()
    
    return chimera_list

def chimeraListToInteractionArray(chimera_list, segments, segment_lengths, chunk_size=1):
    """Converts a list of interactions into a list of numpy arrays

    chunkSize parameter determines the resolution of analysis :
    that is, a chunk size of 100, means that all fragments that fall into 100 nt windows
    will be combined

    segment names and lengths must be passed as a list to enable chunking
    """

    # Initialise interaction map array
    list_of_interaction_arrays = []
    chunks=[math.ceil(x/chunk_size) for x in segment_lengths]

    for i, segment1 in enumerate(segments):
        row = []
        for j, segment2 in enumerate(segments):
            row.append(np.zeros((chunks[i], chunks[j]))) # a numpy array sized for the segment - segment interaction
        list_of_interaction_arrays.append(row)

    # Fill interaction array
    for interaction in chimera_list: # go through each segment interaction in the list
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
                
                if x>=interaction_array.shape[0]:
                    x=interaction_array.shape[0]-1
                if y>=interaction_array.shape[1]:
                    y=interaction_array.shape[1]-1
                    
                interaction_array[x,y]+=1

        list_of_interaction_arrays[seg1][seg2]=interaction_array # put incremented array back into the list

    return list_of_interaction_arrays


def loadInteractionArray(fileName, wd=""):
    with open(os.path.join(wd, fileName), 'rb') as f:
        return np.load(f)

def ChimeraListToCleavageArray(chimera_list, segments, segment_lengths, chunk_size=1):
    """Converts a list of interactions into a list of numpy arrays

    chunkSize parameter determines the resolution of analysis :
    that is, a chunk size of 100, means that all fragments that fall into 100 nt windows
    will be combined

    segment names and lengths must be passed as a list to enable chunking
    """

    # Initialise
    cleavage_maps = {}
    chunks=[math.ceil(x/chunk_size) for x in segment_lengths]
    
    for i, segment in enumerate(segments):
        cleavage_maps[segment] = np.zeros(chunks[i]) # a numpy array sized for the segment 

    # Fill
    for interaction in chimera_list: # go through each segment interaction in the list
        #extract name, start and end position for each fragment
        seg1 = interaction[0]['name']
        start1 = interaction[0]['start']
        end1 = interaction[0]['end']
        seg2 = interaction[1]['name']
        start2 = interaction[1]['start']
        end2 = interaction[1]['end']

        # fill cleavage map
        try:
            x1 = math.ceil(start1/chunk_size)-1
            y1 = math.ceil(end1/chunk_size)-1
            cleavage_maps[seg1][x1]+=1
            cleavage_maps[seg1][y1]+=1
            
            x2 = math.ceil(start2/chunk_size)-1
            y2 = math.ceil(end2/chunk_size)-1
            cleavage_maps[seg2][x2]+=1
            cleavage_maps[seg2][y2]+=1
        except:
            pass
            #print("error: ", seg1, start1, end1, seg2, start2, end2)

    return cleavage_maps

def ChimeraListToMaps(chimera_list, segments, segment_lengths, chunk_size=1):
    """Converts a list of interactions into a list of numpy arrays

    chunkSize parameter determines the resolution of analysis :
    that is, a chunk size of 100, means that all fragments that fall into 100 nt windows
    will be combined

    segment names and lengths must be passed as a list to enable chunking
    """

    # Initialise
    interaction_maps = []
    cleavage_maps = []
    chunks=[math.ceil(x/chunk_size) for x in segment_lengths]
    
    for i, segment1 in enumerate(segments):
        row = []
        for j, segment2 in enumerate(segments):
            row.append(np.zeros((chunks[i], chunks[j]))) # a numpy array sized for the segment - segment interaction
        interaction_maps.append(row)
        cleavage_maps.append(row)

    # Fill
    for interaction in chimera_list: # go through each segment interaction in the list
        #extract name, start and end position for each fragment
        seg1 = segments.index(interaction[0]['name'])
        start1 = interaction[0]['start']
        end1 = interaction[0]['end']
        seg2 = segments.index(interaction[1]['name'])
        start2 = interaction[1]['start']
        end2 = interaction[1]['end']

        # pull out the interaction/cleavage map corresponding to the two segments
        interaction_map = interaction_maps[seg1][seg2] 
        cleavage_map = cleavage_maps[seg1][seg2]

        # fill interaction map
        for i in range(start1, end1, chunk_size): # assign the read to position in the array (increment position), taking into account chunk-size.
            for j in range(start2, end2, chunk_size):
                x = math.ceil(i/chunk_size)-1
                y = math.ceil(j/chunk_size)-1
                
                if x>=interaction_map.shape[0]:
                    x=interaction_map.shape[0]-1
                if y>=interaction_map.shape[1]:
                    y=interaction_map.shape[1]-1
                    
                interaction_map[x,y]+=1

        interaction_maps[seg1][seg2]=interaction_map # put incremented array back into the list
                
        # fill cleavage map
        x = math.ceil(start1/chunk_size)-1
        y = math.ceil(end1/chunk_size)-1
        cleavage_map[x,y]+=1
        x = math.ceil(start2/chunk_size)-1
        y = math.ceil(end2/chunk_size)-1
        cleavage_map[x,y]+=1

        cleavage_maps[seg1][seg2]=cleavage_map # put incremented array back into the list

    return interaction_maps, cleavage_maps

def applyThresholdToInteractionArray(interaction_array, threshold):
    interaction_array[interaction_array<threshold]=0

    return interaction_array

def applyPercentThresholdToInteractionArray(interaction_array, threshold):
    maxval = np.amax(interaction_array)
    log("max val: %i" % maxval)
    thresh = maxval*threshold
    log("threshold: %i" % thresh)
    
    interaction_array[interaction_array<thresh]=0

    # check that original array is not modified
    return interaction_array

def applyThresholdToInteractionList(chimera_list, list_of_interaction_arrays, segments, chunk_size, threshold):
    """Make a new list of interactions, containing only interactions at positions that meet a defined threshold

    chimera_list is a list of the interactions
    list_of_interaction_arrays is a list of arrays containing interaction counts for a given chunk size
    threshold is the minimum number counts at that position for a read to be taken"""

    chimera_list_with_threshold = []

    for interaction in chimera_list: # go through each segment interaction in the list
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
                chimera_list_with_threshold.append(interaction)
                break

    return chimera_list_with_threshold

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

##################
# Plot Functions #
##################

def makePlotsFromInteractionFiles(interaction_files, segments, segment_lengths, condition, cs, t, maxToLoad=np.inf, outdirectory=""):
    
    for interaction_file in interaction_files:
        log("Loading "+interaction_file)
        list_of_interactions = loadChimeraFile(interaction_file, segments, maxToLoad=maxToLoad)
        list_of_interaction_arrays = chimeraListToInteractionArray(list_of_interactions, segments, segment_lengths, cs)
        makePlotsFromInteractionArrays(list_of_interaction_arrays, segments, segment_lengths, condition, cs, t, outdirectory=outdirectory)

def makePlotsFromInteractionArrays(interaction_arrays, segments, segment_lengths, condition, cs, t, outdirectory=""):
    
    for segment1 in segments:
        for segment2 in segments:
            # only plot for segments in alphbetical order (to prevent blank plots)
            if segment1<segment2:
                log("Analysing...")
                seg1 = segments.index(segment1)
                seg2 = segments.index(segment2)
                interaction_array = interaction_arrays[seg1][seg2]
                interactions = applyThresholdToInteractionArray(interaction_array, int(t))
                log("Plotting...")
                #log_interactions = np.log(interactions)
                #log_interactions[log_interactions==-np.inf]=0
                plotTitle=plotname+"_"+segment1+"/"+segment2+": chunk size "+str(cs)+" threshold: "+str(t)
                plotFile=os.path.join(outdirectory, plotname+"_"+segment1+"_"+segment2+"_cs"+str(cs)+"_t"+str(t))
                PlotInteractionMatrix(interactions, title=plotTitle, fileName=plotFile)
                #PlotInteractionMatrix(log_interactions, condition+"_log", segment1, segment2, chunk_size=int(cs), threshold=int(t))


def PlotInteractionMatrix(interaction_array, title="", showPlot=True, fileName=None, figuresize=[20,20]):
    
    plt.figure(figsize=figuresize)
    im = plt.imshow(interaction_array)

    x_dimension = np.size(interaction_array, 1)
    y_dimension = np.size(interaction_array, 0)

    xticks = np.arange(0, x_dimension, 200)
    yticks = np.arange(0, y_dimension, 200)

    ax = plt.gca()
    ax.set_title(title)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, rotation='vertical')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)

    cbar = ax.figure.colorbar(im)
    cbar.ax.set_ylabel(title, rotation=-90, va="bottom")

    if fileName:
        plt.savefig(fileName, format='pdf')
    if showPlot:
        plt.show()
        
def PlotInteractionMatrixs(interaction_arrays, titles, showPlot=True, fileName=None, figuresize=[20,20]):
    
    fig, axs = plt.subplots(1, len(interaction_arrays), figsize=figuresize)
    
    for ax, title, data in zip(axs, titles, interaction_arrays):
        ax.grid(True)

        x_dimension = np.size(data, 1)
        y_dimension = np.size(data, 0)

        xticks = np.arange(0, x_dimension, 200)
        yticks = np.arange(0, y_dimension, 200)

        ax.set_title(title)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, rotation='vertical')
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)
        
        ax.imshow(data, cmap="hot")

    if fileName:
        plt.savefig(fileName, format='pdf')
    if showPlot:
        plt.show()
