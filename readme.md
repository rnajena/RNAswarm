# WORK IN PROGRESS
While working with RNAswarm, you might fing bugs, some are known, while others remain undiscovered. We encourage the community to contribute by reporting any issues they encounter on GitHub. Feel free to reach out to me via email or open an issue directly. It's important to note that I cannot be held responsible for any results obtained using RNAswarm or any conclusions drawn from them.

***
# RNAswarm
RNAswarm is a tool for analyzing SPLASH data. It is a Nextflow pipeline that:
- Trims and aligns chimeric reads to a reference genome.
- Generate heatmaps of interactions between chromosomes and segments (in case of segmented viruses).
- Compare replicates and different datasets in order to identify differentially structured regions.
- Generate a summary of the results.
***

## Installation (only tested on Linux)
### Install Miniconda (or Conda)
- Download the latest version of Miniconda from [Miniconda](https://conda.io/miniconda.html).
```bash
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
- Run the installer.
```bash
$ bash Miniconda3-latest-Linux-x86_64.sh
```
- Add the channel `conda-forge`.
```bash
$ conda config --add channels conda-forge
```

### Create a Conda environment for Nextflow
- Create a new Conda environment for RNAswarm
```bash
$ conda create -n RNAswarm -c bioconda nextflow
```
- Activate the environment.
```bash
$ conda activate RNAswarm
```

### Clone this repository
- Clone the repository to a place where you can run it.
```bash
$ git clone https://github.com/gabriellovate/RNAswarm.git
```


## Usage

### Creating a sample sheet
- Create a sample sheet file with the following columns:

```
sample01,sample01.fastq,reference01.fasta,group01
sample02,sample02.fastq,reference01.fasta,group01
sample03,sample03.fastq,reference02.fasta,group02
sample04,sample04.fastq,reference02.fasta,group02
sample05,sample05.fastq,reference03.fasta,group03
sample06,sample06.fastq,reference03.fasta,group03
```

```bash
# Activate the environment
$ conda activate RNAswarm
# Run the pipeline
$ 
```

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Install [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines.

3. Download the pipeline and test it on a minimal dataset with a single command (under construction):

   ```bash
   nextflow run gabriellovate/RNAswarm -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `conda`, `local`, `slurm`.

4. Start running your own analysis (under construction)!

   ```bash
   nextflow run main.nf -profile slurm,conda\
                        --output <OUTDIR>\
                        --samples <SAMPLES_CSV_FILE>\
                        --annotation_table <ANNOTATION_TABLE_FILE>\
                        --slurm_queue <SLURM_QUEU_AVAILABLE>
   ```

## Documentation

RNAswarm comes with documentation about the pipeline [usage](docs/usage.md) and [output](docs/outpud.md) (under construction).

## Contributions and Support



## Citations

If you use RNAswarm for your analysis, please cite it using the following doi: 
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file (under construction).