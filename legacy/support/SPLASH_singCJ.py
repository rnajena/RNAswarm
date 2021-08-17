# import helper
import sys
#sys.path.append('/pfs/data5/home/fr/fr_fr/fr_cj59/SPLASHanalysis/scripts/')
from .helper_v1CJ import *
 
def makeDir(barcode, datadirectory=""):
    command = " ".join(["mkdir", barcode])
    
    callNonparallelScript([command], wd=datadirectory)
    
def generateIndexes(indexfile, indexdirectory=""):
    command = "/pfs/data5/home/fr/fr_fr/fr_cj59/miniconda3/bin/bwa index "+os.path.join(indexdirectory, indexfile)+".fasta"
    log(command)
    
    callNonparallelScript([command], wd=indexdirectory)  
    
def trim(barcode, quality, noN, minLength, datadirectory="", outdirectory="", logging=True):
    
    # quality trimming
    if logging: log("Quality trimming...q%i" % quality)

    if noN:
        options = "--nextseq-trim="+str(quality)+" --max-n 0 -u 3 -a AAAAAAAAAA"
        prefix="splash_q"+str(quality)+"_noN_min"+str(minLength)+"_"
    else:
        options = "--nextseq-trim="+str(quality)+" -a AAAAAAAAAA"
        prefix="splash_q"+str(quality)+str(minLength)+"_"

    outfile = os.path.join(outdirectory, prefix+barcode+".fastq")
           
    callNonparallelScript(["cutadapt "+options+" -o "+outfile+" "+barcode+".fastq"], wd=datadirectory, logging=logging)
    
    if logging: log("Finished quality trimming...")
        
    return prefix

def alignRead(barcode, indexfile,  indexdirectory="", prefix="", datadirectory="", logging=True):
    index = os.path.join(indexdirectory, indexfile)+".fasta"
    
    # align with bwa            
    command = " ".join(["/pfs/data5/home/fr/fr_fr/fr_cj59/miniconda3/bin/bwa",
               "mem",
               "-T 20",
               index,
               prefix+barcode+".fastq"])
    
    if logging: log("Aligning with command: "+command)
    
    stdout, stderr = callNonparallelScript([command], wd=datadirectory, logging=False)
    writeFile(stdout[0], prefix+barcode+".sam", wd=datadirectory)
    #writeFile(stderr[0], prefix+barcode+"_alignresults.txt", wd=datadirectory)
    
    
    # convert to bam
    if logging: log("convert to bam...")

    command = " ".join(["/pfs/data5/home/fr/fr_fr/fr_cj59/miniconda3/bin/samtools",
               "view",
               "-S",
               "-b",
               prefix+barcode+".sam"])
    
    stdout, stderr = callNonparallelScript([command], wd=datadirectory, logging=False, encoding=None)
    writeFile(stdout[0], prefix+barcode+".bam", wd=datadirectory, filetype="wb")

    #command = " ".join(["/vol/projects/rsmyth/tools/bwa",
    #           "mem",
    #           "-T 20",
    #           index,
    #           prefix+barcode+".fastq",
    #           ">",
    #           prefix+barcode+".sam"])

    
    #callNonparallelScript([command], wd=datadirectory, encoding=None)    
    #callNonparallelScript(["/vol/biotools/bin/samtools view -S -b "+prefix+barcode+".sam > "+prefix+barcode+".bam"], wd=datadirectory, encoding=None)

def generateChimFile(barcode, prefix="", fourReadChemistry=True, datadirectory=""):

    callNonparallelScript(["/home/fr/fr_fr/fr_cj59/miniconda3/envs/py2/bin/python /pfs/data5/home/fr/fr_fr/fr_cj59/lal_scratch/0_Working/newgenseq/support/find_chimeras_rs.py -i "+prefix+barcode+".bam -o "+prefix+barcode+".chim"], wd=datadirectory)

        
def Splash_sing(barcodes, indexfile, quality, noN, minLength, indexdirectory="", datadirectory="", outdirectory="", logging=True):
    
    
    #0
    generateIndexes(indexfile, indexdirectory)
    
    for i, barcode in enumerate(barcodes):
        outdir = os.path.join(outdirectory, str(i)+"_"+barcode)
        #1
        makeDir(outdir, datadirectory=datadirectory)  
        #2
        prefix = trim(barcode, quality, noN, minLength, outdirectory=outdir, datadirectory=datadirectory, logging=logging)
        #3
        alignRead(barcode, indexfile, indexdirectory=indexdirectory, prefix=prefix, datadirectory=outdir, logging=logging)
        #5
        generateChimFile(barcode, prefix=prefix, fourReadChemistry=True, datadirectory=outdir)
        #6
        #chimStats(genome_stats, match_stats)
