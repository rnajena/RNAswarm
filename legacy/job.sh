#! /bin/bash

#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --mem=100000mb

module load devel/miniconda
conda activate ngs
srun python /home/fr/fr_fr/fr_cj59/lal_scratch/0_Working/newgenseq/main.py /home/fr/fr_fr/fr_cj59/lal_scratch/0_Working/newgenseq/run.conf
