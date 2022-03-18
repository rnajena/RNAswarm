# RNAswarm
RNAswarm is a tool for analyzing SPLASH data. It is a Nextflow pipeline that:
- Trims and aligns chimeric reads to a reference genome.
- Generate heatmaps of interactions between chromosomes and segments (in case of segmented viruses).
- Compare replicates and different datasets in order to identify differentially structured regions.
- Generate a summary of the results.

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
