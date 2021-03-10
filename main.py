import os
from datetime import datetime
import configparser
import argparse

import pytoml as toml

from support import log
from alignment import Aligner
from support.helper_v1CJ import timeStamp
from chim_handler import ChimerHandler


class Configuration:
    def __init__(self, config):
        self.project = config.get("General", "Project_name")
        self.user = config.get("General", "User")
        self.logdirectory = config.get("General", "LogDirectory")
        self.datadirectory = config.get("General", "DataDirectory")
        self.dataid = config.get("General", "DataID")
        self.indexfile = config.get("General", "IndexFile")
        self.barcodes = list(dict.fromkeys([barcode.strip()\
                             for barcode in config.get("General", "Barcodes").split(',')]))
        
        self.generate_index = config.getboolean("Run", "Generate_index")
        self.align = config.getboolean("Run", "Align")
        self.chim_handle = config.getboolean("Run", "Chim_handle")
        self.plot = config.getboolean("Run", "Plot")
        self.log = config.getboolean("Run", "Logging")
        self.log_level = config.get("Run", "Logging_level")

        self.old_timestamp = config.get("Re-run", "Old_timestamp")
    

def write_toml(file_fullname, dictionary):
    log.info("[%s] Writing toml %s" %(timeStamp(), file_fullname))
    try:
        with open(file_fullname, "w") as fi:
            toml.dump(dictionary, fi)
            log.info("[%s] Successfully written the toml." %timeStamp())
        return True
    except:
        log.info("[%s] Failure in writing the toml." %timeStamp())
        return False


def read_toml(file_fullname):
    log.info("[%s] Reading the toml %s" %(timeStamp(), file_fullname))
    try:
        with open(file_fullname, "rb") as fi:
            dictionary = toml.load(fi)
        log.info("[%s] Successfully read the toml." %timeStamp())

    except:
        raise Exception("[%s] Failure in reading the toml." %timeStamp())

    return dictionary



def setup_logger(logging, level, project, user, directory, old_timestamp=None):

    now = datetime.now()
    key_timestamp = now.strftime("%Y%m%d_%I%M%S")

    log.set_log_level(level)
    log.info()
    log.info()
    log.info("\tProject: %s" %project)
    log.info("\tUser: %s" %user)
    log.info()

    if logging:
        # then a log file is needed
        log_directory = os.path.join(directory, "log")
        
        if not os.path.exists(log_directory):
            # create the log file in the output directory
            os.makedirs(log_directory)

     
        logfile_name = os.path.join(log_directory, key_timestamp + ".txt")
        log.info("[%s] The logging is active" %timeStamp())
        log.info("[%s] Logging file %s" %(timeStamp(), logfile_name))
        log.add_log_file(logfile_name)

    else:
        # No logging
        log.info("[%s] The logging is not active" %timeStamp())


    return key_timestamp


def get_chim_files_list_from_dictionary(alignment_dictionary, barcodes):
    log.info()
    log.info("[%s] Collecting the Chimer files list from the alignement dictionary." %timeStamp())
    log.info()
    chimer_files = list()

    for i, barcode in enumerate(barcodes):
        log.info("\t%s. For barcode %s:" %(i, barcode))

        barcodes_in_toml = alignment_dictionary["out"].keys()
        if barcode in barcodes_in_toml:
            _file = alignment_dictionary["out"][barcode]["chim"]
            log.info("\t\tchimer file %s" %_file)
            chimer_files.append(_file)
        
        else:
            log.info("\t\tSKIPPING.. No chimer file in dictionary for %s" %barcode)
    
    return chimer_files


def init():
    parser = argparse.ArgumentParser(description="RNA Sequence analyzer")
    parser.add_argument('conf', type=str, help="run.conf is necessary")
    
    args = parser.parse_args()

    if os.path.exists(args.conf):
        config = configparser.ConfigParser()
        config.read(args.conf)
        
        configuration = Configuration(config)
        return configuration

    else:
        raise Exception("No configuration file.")


def main(configuration):
    # General configuration
    project = configuration.project
    user = configuration.user

    old_timestamp = configuration.old_timestamp
    logdirectory = configuration.logdirectory
    datadirectory = configuration.datadirectory
    dataid = configuration.dataid
    indexfile = configuration.indexfile
    barcodes = configuration.barcodes
    
    # Run configuration
    Index = configuration.generate_index
    Align = configuration.align
    Analyze_chim = configuration.chim_handle
    Plot = configuration.plot
    Logging = configuration.log
    Level = configuration.log_level
    
    # Re-run configuration
    old_timestamp = configuration.old_timestamp
 
    # Create the timestamp of this run
    # for folder creation and so on
    key_timestamp = setup_logger(logging=Logging, level=Level, project=project,\
                                 user=user, directory=logdirectory,\
                                 old_timestamp=old_timestamp) 

    # initiate the aligner object
    aligner = Aligner(directory=datadirectory, dataid=dataid, runname=key_timestamp,\
                      indexfile=indexfile, barcodes=barcodes, old_timestamp=old_timestamp)

    if Index:
        log.info("[%s] Generating the indexes." %timeStamp())
        aligner.generate_indexes()
    else:
        log.info("[%s] Skipping the index generation." %timeStamp())

    # starting to align the data
    if Align:
        log.info()
        log.info("[%s] Alignment starting." %timeStamp())
        
        #TODO
        aligner.align()

        log.info("[%s] Aligning ending." %timeStamp())
        log.info()
        alignment_dictionary = aligner.alignment_dictionary

        # fetch the chimer files from the dictionary
        chim_list = get_chim_files_list_from_dictionary(aligner.alignment_dictionary, barcodes=aligner.barcodes)
        
        # analyze the chimeras
        # and produce npy - this is not controllable
        if Analyze_chim:
            log.info("[%s] Initialting chimer files analysis." %timeStamp())

            analyzer = ChimerHandler(chim_list)

            # plotting can be controlled
            # currently only pdf format
            if Plot:
                log.info("\t With plotting.")
                #TODO
                analyzer.save_arrays(plot=True)

            else:
                log.info("\t Without plotting.")
                #TODO
                analyzer.save_arrays(plot=False)

    else:
        
        if old_timestamp != "None":
            log.info("\t\t\tOld timestamp provided is %s" %old_timestamp)

            old_run_directory = os.path.join(datadirectory, dataid, old_timestamp)

            if os.path.isdir(old_run_directory):
                log.info("[%s] Considering the old run from %s" %(timeStamp(), old_run_directory))
                log.info()
            else:
                old_run_directory = None
        else:
            old_run_directory = None

        if old_run_directory is not None:
            old_run_alignment_dictionary = os.path.join(old_run_directory, "alignment_dictionary.toml")

            if os.path.exists(old_run_alignment_dictionary):
                alignment_dictionary = read_toml(old_run_alignment_dictionary)
                # fetch the chimer files from the dictionary
                chim_list = get_chim_files_list_from_dictionary(alignment_dictionary, barcodes=barcodes)

                if Analyze_chim:
                    log.info("[%s] Initialting chimer files analysis." %timeStamp())

                    analyzer = ChimerHandler(chim_list)

                    # plotting can be controlled
                    # currently only pdf format
                    if Plot:
                        log.info("\t\tWith plotting.")
                        #TODO
                        analyzer.save_arrays(plot=True)

                    else:
                        log.info("\t\tWithout plotting.")
                        #TODO
                        analyzer.save_arrays(plot=False)
        else:
            raise Exception("Old timestamp for re-run is invalid")

    # Writing the toml
    # A must do task - not to be controlled through the run.conf
    toml_fullname = os.path.join(aligner.outdirectory, "alignment_dictionary.toml")
    write_toml(file_fullname=toml_fullname, dictionary=alignment_dictionary)


if __name__=="__main__":
    
    configuration = init()
    main(configuration)
 
