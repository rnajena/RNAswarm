# WORK IN PROGRESS
While working with RNAswarm, you might fing bugs, some are known, while others remain undiscovered. We encourage the community to contribute by reporting any issues they encounter on GitHub. Feel free to reach out to me via email or open an issue directly. It's important to note that I cannot be held responsible for any results obtained using RNAswarm or any conclusions drawn from them.

***
# RNAswarm
RNAswarm is a tool for analyzing SPLASH data. It is a Nextflow pipeline that:
- Trims (with `fastp`) and maps (with `segemehl`) chimeric reads to a reference genome.
- Generate heatmaps of interactions between viral segments or/and host transcripts.
- Compare groups or mutants to identify differentially structured regions (with `DESeq2`)
- Generate circos plots of potential interactions and potential differentially structured regions.
***

## Dependencies and installation
The pipeline is written in Nextflow, which can be used on any POSIX compatible system (Linux, OS X, etc). Windows system is supported through WSL2. You need Nextflow installed, conda and apptainer.
1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
2. Install [miniconda](https://docs.conda.io/projects/miniconda/en/latest/)
3. Install [apptainer](https://apptainer.org/docs/admin/main/installation.html#install-unprivileged-from-pre-built-binaries) (from our experience, installing from the pre-built binaries is the easiest way to go)
4. With nextflow installed, you only have to run the pipeline, if you would like to cache the pipeline without running it you can use the following command:
```bash
nextflow pull gabriellovate/RNAswarm
```

In the future, only either conda or apptainer will be needed.

## Test your installation
This command retrieves the project from GitHub and executes the pipeline using the included test data. It utilizes the default test profile, which employs the provided [sample sheet](https://github.com/gabriellovate/RNAswarm/blob/main/data/samples.csv) and [comparisons sheet](https://github.com/gabriellovate/RNAswarm/blob/main/data/comparisons.csv), both located in the [`data/` directory](https://github.com/gabriellovate/RNAswarm/tree/main/data).

```bash
nextflow run gabriellovate/RNAswarm \
            -profile test,apptainer \
            --output <OUTDIR>
```

## Update RNAswarm to the latest version
After you `nextflow run`, the pipeline is cached by Nextflow. To update to the latest version, you can run the following command:
```bash
nextflow pull gabriellovate/RNAswarm
```

## Usage
RNAswarm takes as input a sample sheet and a comparisons sheet. The sample sheet contains the information about the samples to be analyzed, while the comparisons sheet contains the groups to be compared.

### Creating a sample sheet
Create a sample sheet file with the following columns:

```csv
sample01,sample01.fastq,reference01.fasta,group01
sample02,sample02.fastq,reference01.fasta,group01
sample03,sample03.fastq,reference02.fasta,group02
sample04,sample04.fastq,reference02.fasta,group02
sample05,sample05.fastq,reference03.fasta,group03
sample06,sample06.fastq,reference03.fasta,group03
```

### Creating a comparisons sheet
Create a comparisons sheet file with the groups for differential analysis:
```csv
group01,group02
group01,group03
group02,group03
```

### Running the pipeline
#### with local executor
```bash
nextflow run gabriellovate/RNAswarm \
            -profile local,apptainer \
            --samples <SAMPLES_CSV_FILE> \
            --comparisons <COMPARISONS_CSV_FILE> \
            --output <OUTDIR> \
```

#### with slurm
```bash
nextflow run gabriellovate/RNAswarm \
            -profile slurm,apptainer \
            --slurm_queue <SLURM_QUEU_AVAILABLE> \
            --samples <SAMPLES_CSV_FILE> \
            --comparisons <COMPARISONS_CSV_FILE> \
            --output <OUTDIR> \
```

## Cite us
If you use RNAswarm for your analysis, please cite our github repository.

```bibtex
@software{Lencioni_Lovate_RNAswarm,
author = {Lencioni Lovate, Gabriel and Lamkiewicz, Kevin},
license = {MIT},
title = {{RNAswarm}},
url = {https://github.com/gabriellovate/RNAswarm}
}
```
