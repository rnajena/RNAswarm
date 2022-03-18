# RNAswarm
RNAswarm is a tool for analyzing SPLASH data. It is a Nextflow pipeline that:
- Trims and aligns chimeric reads to a reference genome.
- Generate heatmaps of interactions between chromosomes and segments (in case of segmented viruses).
- Compare replicates and different datasets in order to identify differentially structured regions.
- Generate a summary of the results.

## Installation (only tested on Linux)
### Install Miniconda (or Conda)
- Download the latest version of Miniconda from [Miniconda](https://conda.io/miniconda.html).
- Run the installer.

### Setup channel order
- Add the channel `bioconda`.
- Add the channel `conda-forge`.

### Create a Conda environment for Nextflow
- Create a new Conda environment.
```bash
conda create -n rnaswarm
```
- Activate the environment.
```bash
source activate rnaswarm
```
- Install Nextflow.
```bash
conda install -c bioconda nextflow
```

### Clone this repository
- Clone the repository.
```bash
git clone https://github.com/gabriellovate/RNAswarm.git
```

## Usage
