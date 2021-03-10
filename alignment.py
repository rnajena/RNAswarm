import os

from support import helper_v1CJ as helper
from support import SPLASH_singCJ as SPLASHer
from support import RRseq_plottingCJ as plotter
from support import log

makedir = helper.makeDir
indexer = SPLASHer.generateIndexes
singer = SPLASHer.Splash_sing
trimmer = SPLASHer.trim
align_reader = SPLASHer.alignRead
chimfile_generator = SPLASHer.generateChimFile
timeStamp = helper.timeStamp

LOGGING=True


class Aligner:

    def __init__(self, directory, dataid, runname, indexfile, barcodes, quality=20, noN=True, minLength=6,\
                 old_timestamp=None, logging=LOGGING,):

        log.info()
        log.info("[%s] Initializing the Aligner object" %timeStamp())
    
        self.directory = directory
        self.indexdirectory = os.path.join(directory, "indexes")
        self.datadirectory = os.path.join(directory, dataid, "orig")
        self.outdirectory = os.path.join(directory, dataid, runname)
        
        if old_timestamp is not None:
            self.followupdirectory = os.path.join(directory, dataid, old_timestamp)
        else:
            self.followupdirectory = None
            

        if not os.path.isdir(self.outdirectory):
            os.makedirs(self.outdirectory)
        
        self.alignment_dictionary = dict()
        self.alignment_dictionary["general"] = dict()
        self.alignment_dictionary["general"]["directory"] = self.directory
        self.alignment_dictionary["general"]["indexdirectory"] = self.indexdirectory
        self.alignment_dictionary["data"] = dict()
        self.alignment_dictionary["data"]["datadirectory"] = self.datadirectory
        self.alignment_dictionary["out"] = dict()
        self.alignment_dictionary["out"]["outdirectory"] = self.outdirectory

        log.info("\tIndex directory: %s" %self.indexdirectory)
        log.info("\tData directory: %s" %self.datadirectory)
        log.info("\tOut directory: %s" %self.outdirectory)

        self.indexfile = indexfile
        self.barcodes = barcodes
        log.info()
        log.info("\tBarcodes: ")
        for i, barcode in enumerate(self.barcodes):
            log.info("\t\t%s. %s" %(i+1, barcode))
        log.info()
        
        self.runname = runname
        self.quality = quality
        self.noN = True
        self.minLength = minLength
        self.logging = logging

        self.alignment_dictionary["general"]["indexfile"] = self.indexfile
        self.alignment_dictionary["general"]["barcodes"] = self.barcodes
        self.alignment_dictionary["general"]["runname"] = self.runname
        self.alignment_dictionary["general"]["quality"] = self.quality
        self.alignment_dictionary["general"]["noN"] = self.noN
        self.alignment_dictionary["general"]["minLength"] = self.minLength
        self.alignment_dictionary["general"]["logging"] = self.logging

        log.info("\tFurther details:")
        log.info("\t\tName of run: %s" %self.runname)
        log.info("\t\tQuality: %s" %self.quality)
        log.info("\t\tnoN: %s" %str(self.noN))
        log.info("\t\tMin. length: %s" %self.minLength)
        log.info()
 
    
    def generate_indexes(self):
        log.info("[%s] Generating indexes." %timeStamp())
        try:
            indexer(indexfile=self.indexfile, indexdirectory=self.indexdirectory)
            log.info("[%s] Successfully generated the indexes")

        except:
            log.error("Failed to generate the indexes")


    def __custom_singer__(self, flag_trimmer=False, flag_align_reader=False, flag_chimfile_generator=False):

        if flag_trimmer or flag_align_reader or flag_chimfile_generator:
            for i, barcode in enumerate(self.barcodes):
                self.alignment_dictionary["out"][barcode] = dict()

                outdir = os.path.join(self.outdirectory, str(i) + "_" + barcode)
                self.alignment_dictionary["out"][barcode]["outdir"] = outdir

                makedir(outdir, datadirectory=self.outdirectory)
                
                if flag_trimmer:
                    try:
                        prefix = trimmer(barcode=barcode, quality=self.quality, noN=self.noN, minLength=self.minLength,\
                                         outdirectory=outdir, datadirectory=self.datadirectory, logging=self.logging)
                        
                        self.alignment_dictionary["out"][barcode]["prefix"] = prefix
                        self.alignment_dictionary["out"][barcode]["fastq"] = os.path.join(outdir, prefix + barcode + ".fastq")

                    except:
                        os.remove(os.path.join(outdir, prefix + barcode + ".fastq"))
                        log.error("[%s] Unsuccessful Trimming." %timeStamp())
                        
                else:
                    
                    if self.noN:
                        prefix="splash_q" + str(self.quality) + "_noN_min" + str(self.minLength) + "_"
                    else:
                        prefix="splash_q" + str(self.quality) + str(self.minLength) + "_"
                
                    self.alignment_dictionary["out"][barcode]["prefix"] = prefix
                    self.alignment_dictionary["out"][barcode]["fastq"] = os.path.join(outdir, prefix + barcode + ".fastq")


                if flag_align_reader:
                    try:
                        align_reader(barcode=barcode, indexfile=self.indexfile, indexdirectory=self.indexdirectory,\
                                     prefix=prefix, datadirectory=outdir, logging=self.logging)
                        self.alignment_dictionary["out"][barcode]["bam"] = os.path.join(outdir, prefix + barcode + ".bam")
                        self.alignment_dictionary["out"][barcode]["sam"] = os.path.join(outdir, prefix + barcode + ".sam")
                        log.info("[%s] The .bam and .sam files are ready." %timeStamp())
                        
                    except:
                        os.remove(os.path.join(outdir, prefix + barcode + ".bam"))
                        os.remove(os.path.join(outdir, prefix + barcode + ".sam"))
                        log.error("[%s] Unsuccessful alignement reading." %timeStamp())
                        
                else:
                    self.alignment_dictionary["out"][barcode]["bam"] = os.path.join(outdir, prefix + barcode + ".bam")
                    self.alignment_dictionary["out"][barcode]["sam"] = os.path.join(outdir, prefix + barcode + ".sam")


                if flag_chimfile_generator:
                    try:
                        chimfile_generator(barcode=barcode, prefix=prefix, fourReadChemistry=True, datadirectory=outdir)
                        self.alignment_dictionary["out"][barcode]["chim"] = os.path.join(outdir, prefix + barcode + ".chim")
                        log.info("[%s] The .chim file is ready." %timeStamp())

                    except:
                        os.remove(os.path.join(outdir, prefix + barcode + ".chim"))
                        log.error("[%s] Unsuccessful chimfile generation." %timeStamp())

    
    def trim(self):
        self.__custom_singer__(flag_trimmer=True)


    def read_alignment(self):
        self.__custom_singer__(flag_align_reader=True)


    def generate_chimfile(self):
        self.__custom_singer__(flag_chimfile_generator=True)


    def align(self):
        # singeer is actually :
        # 1. makeDir()
        # 2. trim()
        # 3. alignRead()
        # 4. generateChimFile()
        self.__custom_singer__(flag_trimmer=True, flag_align_reader=True, flag_chimfile_generator=True)


    def write(self):
        log.info("Starting to save the interaction arrays as npy files.")
        plotter.saveInteractionArrays(interaction_files, segments, segment_lengths, cs, maxToLoad,\
                                      outdirectory=self.npydirectory)
        log.info("Finished saving the npys.")

        log.info("Starting to save the interactions as plots.")
        plotter.makePlotsFromInteractionFiles()
        log.info("Finished saving the plots")


    def read(self):
        pass



if __name__=="__main__":
    aligner = Aligner(directory="/home/fr/fr_fr/fr_cj59/lal_scratch/1_Data", name="samplerun",\
                      indexfile="SC35MwtMulti", barcodes=["R1_SC35Mwt_cat"])

