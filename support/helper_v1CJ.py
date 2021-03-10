####################
# Helper functions #
####################

from re import search
from datetime import datetime
from subprocess import Popen, call, PIPE
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

import sys
import os.path
import pysam

import numpy as np
from itertools import groupby
from heapq import merge

from . import log as _log

#####################
# Logging functions #
#####################

# Formatted time stamp
def timeStamp():
    return datetime.now().strftime("%c")

# Print command with time stamp
def log(command):
    _log.info("[%s] %s" % (timeStamp(), command))

##################
# Iter functions #
##################

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

    
###################
# Function calls  #
###################

# Non parallel function call
def callNonparallelFunction(f, args):
    for item in args:
        log("starting function call (%s) with args: %s" % (f.__name__, item))
        f(item)


# Call external scripts
def runCommand(command=None, stdin=PIPE, stderr=PIPE, stdout=PIPE, wd=None, pipe=False, backgroundProcess=False, encoding='utf8'):
    proc = Popen(command.split(" "), cwd=wd, stdin=stdin, stdout=stdout, stderr=stderr, encoding=encoding)
    
    if pipe:
        return proc.stdout
    elif backgroundProcess:
        return proc
    else:
        stdout, stderr = proc.communicate()
        return stdout, stderr

# Non parallel script call, with logging message
def callNonparallelScript(args, wd=None, logging=True, encoding='utf8'):
    stdouts = []
    stderrs = []
    
    for item in args:
        if logging: log("running command: %s" % item)
        stdout, stderr = runCommand(item, wd=wd, encoding=encoding) 
        
        stdouts.append(stdout)
        stderrs.append(stderrs)
        
        if logging:
            if stderr:
                log("stderr from command: %s\n\n%s" % (item, stderr))
            else:
                log("stdout from command: %s\n\n%s" % (item, stdout))
        else:
            return stdouts, stderrs
    
# Call scripts in parallel, with logging messages
def callParallelScript(commands, processname, maxProcesses=2, logging=True):
    log("starting parallel analysis: %s" % (processname))

    pool = Pool(maxProcesses)
    stdouts, stderrs = zip(*pool.map(runCommand, commands))
    pool.close()
    pool.join()
    
    if logging:            
        for command, stdout, stderr in zip(commands, stdouts, stderrs):
            if stderr:
                log("error from parallel command: %s\n\n%s" % (command, stderr))
            else:
                log("output from parallel command: %s\n\n%s" % (command, stdout))
    else:
        return stdouts, stderrs

def makeDir(name, datadirectory=""):
    command = " ".join(["mkdir", name])
    
    callNonparallelScript([command], wd=datadirectory)
    
#################
# Fastq and Sam #
#################

# Parse fasta file
def getFastaSeq(file, recordId, wd=""):
    """
        Read in the reference sequence used for the alignment
    """
    
    for record in SeqIO.parse(os.path.join(wd, file), "fasta"):
        if record.id == recordId:
            return str(record.seq).upper()

def bcl2fastq(wd, fastqForIndexReads=None):
    
    fastq = []
    
    if fastqForIndexReads:
        command = "/vol/biotools/bin/bcl2fastq --create-fastq-for-index-reads"
    else:
        command = "/vol/biotools/bin/bcl2fastq"

    stdout, stderr = runCommand(command, wd=wd)
    log("running command: %s" % command)
    log("output from  command: %s\n\n%s" % (command, stderr))
    
    for line in stderr.split("\n"):
        if search("Created FASTQ", line):
            fastq.append(line.split('"')[1])
    
    return fastq

def countReadsInFastQ(file, wd=""):
    extension = file.split(".")[-1]

    if extension=="gz":
        #log("decompressing and counting reads in fastq file: "+file)
        stdout = runCommand("gunzip -c "+file, pipe=True, encoding=None, wd=wd)
        stdout, stderr = runCommand("wc -l", stdin=stdout)
    elif extension=="fastq":
        #log("counting reads in fastq file: "+file)
        stdout, stderr = runCommand("wc -l "+file, wd=wd)
    
    numberOfReads = int(stdout.split(" ")[0])/4
    #log("number of reads in fastq file "+file+": "+str(numberOfReads))
    
    return numberOfReads

def getFastqReadLengths(file, wd=""):
    readLengths = []
    records = list(SeqIO.parse(os.path.join(wd, file), "fastq"))
    
    readLengths = [len(record) for record in records]

    return readLengths

def qualityTrim(rx, rb, ry, ra, maxN=0, minLength=10, prefix="", strict=True, wd=None):
    readLengths = []
    success, containsN, tooShort = 0, 0, 0
    
    rxrecords = SeqIO.parse(os.path.join(wd, rx), "fastq")
    ryrecords = SeqIO.parse(os.path.join(wd, ry), "fastq")
    rarecords = SeqIO.parse(os.path.join(wd, ra), "fastq")
    rbrecords = SeqIO.parse(os.path.join(wd, rb), "fastq")
    
    rxFile = open(os.path.join(wd, prefix+rx), 'w')
    ryFile = open(os.path.join(wd, prefix+ry), 'w')
    raFile = open(os.path.join(wd, prefix+ra), 'w')
    rbFile = open(os.path.join(wd, prefix+rb), 'w')
    
    for rxrecord, ryrecord, rarecord, rbrecord in zip(rxrecords, ryrecords, rarecords, rbrecords):
        
        if strict:
            if len(rxrecord)>=minLength and len(ryrecord)>=minLength and len(rarecord)>=minLength and len(rbrecord)>=minLength:
                if ("N" not in rxrecord.seq) and ("N" not in ryrecord.seq) and ("N" not in rarecord.seq) and ("N" not in rbrecord.seq):
                    SeqIO.write(rxrecord, rxFile, 'fastq')
                    SeqIO.write(ryrecord, ryFile, 'fastq')
                    SeqIO.write(rarecord, raFile, 'fastq')
                    SeqIO.write(rbrecord, rbFile, 'fastq')
                    success+=1
                else:
                    containsN+=1
            else:
                tooShort+=1
        else:
            pass
            # not implemented
                
    print("%i successful, %i containing N, %i too short" % (success, containsN, tooShort))

def writeFile(content, file, filetype="w", wd=""):
    with open(os.path.join(wd, file), filetype) as f:
        f.write(content)

def readFile(file, wd=""):
    with open(os.path.join(wd, file), 'r') as f:
        print(f.read())

def groupAlignments(samfiles): # group aligment identical headers into something we can iterate over
    sorted_samfiles = []
    
    for samfile in samfiles:
        sorted_samfiles.append(sorted(samfile, key=lambda samfile: samfile.qname))

    samfiles = merge(*sorted_samfiles, key=lambda samfile: samfile.qname)
    samfiles_grouped = groupby(samfiles, key=lambda samfile: samfile.qname)

    return samfiles_grouped

def annotateSamfiles(iter, read1=False, read2=False): #annnotate samfiles so they can be grouped by read
    for item in iter:
        item.is_read1 = read1
        item.is_read2 = read2
        yield item

def splitAlignments(iter): # split the groups according to the is_read1 and is_read2 label
    read1 = []
    read2 = []

    for item in iter:
       if item.is_read1:
          read1.append(item)
       elif item.is_read2:
          read2.append(item)
       else:
          raise ValueError('Alignment not annotated with is_read1 and is_read2')

    yield read1, read2
    
#############
# Debugging #
#############


# for debugging; print out fastQ reads for the given read id
def printFastqRead(files, readid, wd=""):
    for file in files:
        for record in SeqIO.parse(os.path.join(wd, file), "fastq"):
            if record.id == readid:
                print(record.seq)
                break                

# for debugging; write test fastQ files for the given read id
def writeFastqRead(filesIn, filesOut, readid, wd=""):
    for i, fileOut in enumerate(filesOut):
        filesOut[i] = open(os.path.join(wd, fileOut), 'w')
    
    for i, file in enumerate(filesIn):
        for record in SeqIO.parse(os.path.join(wd, file), "fastq"):
            if record.id == readid:
                SeqIO.write(record, filesOut[i], 'fastq')
                print(record.seq)
                break

# for debugging; print out sam file infor for the given read id
def printSamRead(samFile, readid, wd=""):
    samF = pysam.Samfile(os.path.join(wd, samFile), 'r')
    count = 0
    
    for alignedread in samF.fetch():
        if alignedread.query_name == readid:
            print("** sam")
            print("** refstart: %i" % alignedread.reference_start)
            print("** refend: %i" % alignedread.reference_end)
            print(alignedread)
            
            if count==1:
                break
            else:
                count+=1
##########
# Helper #
##########

def prepareDirectories():
    pass

def loadSampleSheet():
    pass

def noneIsInfinite(value):
    if value is None:
        return float("inf")
    else:
        return value

def noneIsZero(value):
    if value is None:
        return 0
    else:
        return value       
