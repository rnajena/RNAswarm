import os

import numpy as np
import matplotlib.pyplot as plt

from support import helper_v1CJ as helper
from support import RRseq_plottingCJ as plotter
from support import log

timeStamp = helper.timeStamp


class ChimerHandler:

    def __init__(self, chimfiles):
        self.chimfiles = chimfiles
        self.pdfdirectories = list()
        self.npydirectories = list()
        self.segments = ["HA", "M", "NA", "NP", "NS", "PA", "PB1", "PB2"]
        self.segment_lengths = [1740, 1027, 1461, 1565, 890, 2233, 2341, 2341]
        self.cs=1
        self.t=0
        #log.info("Initiating the chimer handler")
        #log.info("\tcs: %s" %self.cs)
        #log.info("\tt: %s" %self.t)
        #log.info("\tMaxToLoad: %s" %self.maxToLoad)
        #log.info()


    def save_arrays(self, plot=False):
        pass
        
        for _file in self.chimfiles:
            _npy_path = os.path.join(os.path.split(_file)[0], "npy")
            if not os.path.exists(_npy_path):
                os.makedirs(_npy_path)

            self.npydirectories.append(_npy_path)
 
            plotter.saveInteractionArrays(interaction_files=[_file],\
                                          segments=self.segments, segment_lengths=self.segment_lengths,\
                                          cs=self.cs, outdirectory=_npy_path)

            if plot:
                _pdf_path = os.path.join(os.path.split(_file)[0], "plots")
                if not os.path.exists(_pdf_path):
                    os.makedirs(_pdf_path)
                if _pdf_path not in self.pdfdirectories:
                    self.pdfdirectories.append(_pdf_path)
		
                _list_npys = os.listdir(_npy_path)

                for i, _npy in enumerate(_list_npys):
                    if os.path.splitext(_npy)[1] == ".npy":
                        log.info("[%s] Converting %s to pdf" %(timeStamp(), _npy))        
                        
                        interaction_array = np.load(os.path.join(_npy_path,_npy))
                        plotname = os.path.join(_pdf_path, os.path.split(_npy)[1].split(".npy")[0]+".pdf")
 
                        figuresize=[20,20]
                        plt.figure(figsize=figuresize)

                        im = plt.imshow(interaction_array)
                        
                        x_dimension = np.size(interaction_array, 1)
                        y_dimension = np.size(interaction_array, 0)
                        xticks = np.arange(0, x_dimension, 200)
                        yticks = np.arange(0, y_dimension, 200)
                        
                        ax = plt.gca()
                        #title = "segment (" + segment1+")/ segment ("+segment2+") : chunk size = "+str(cs)+" threshold = "+str(t)
                        # TODO
                        title = "TODO: not set"
                        ax.set_title(title)
                        ax.set_xticks(xticks)
                        ax.set_xticklabels(xticks, rotation='vertical')
                        ax.set_yticks(yticks)
                        ax.set_yticklabels(yticks)
                        cbar = ax.figure.colorbar(im)
                        cbar.ax.set_ylabel(title, rotation=-90, va="bottom")      
                        plt.savefig(plotname, format='pdf')
                        plt.close()



if __name__=="__main__":
    sample_chim_list = ["/home/fr/fr_fr/fr_cj59/lal_scratch/1_Data/data0120/20201118_030354/0_R1_SC35Mwt_cat/splash_q20_noN_min6_R1_SC35Mwt_cat.chim"]
    sample_handler = ChimerHandler(sample_chim_list)

    sample_handler.save_arrays(plot=True)
