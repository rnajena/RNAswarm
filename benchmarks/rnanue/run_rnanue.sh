#!/bin/sh

# Setup the environment with singularity
singularity pull docker://cobirna/rnanue

# Run RNAnue on pr8 reads
singularity run --no-home \
-B /data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/reads:/trtms,\
/data/dessertlocal/projects/gl_iav-splash_freiburg/results/benchmarks:/outdir,\
/data/dessertlocal/projects/gl_iav-splash_freiburg/src/scripts_and_notebooks/benchmarks/rnanue:/params,\
/data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/genomes:/genomes \
docker://cobirna/rnanue ./RNAnue/build/RNAnue complete --config /params/params_pr8.cfg

# Run RNAnue on udorn reads
singularity run --no-home \
-B /data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/reads:/trtms,\
/data/dessertlocal/projects/gl_iav-splash_freiburg/results/benchmarks:/outdir,\
/data/dessertlocal/projects/gl_iav-splash_freiburg/src/scripts_and_notebooks/benchmarks/rnanue:/params,\
/data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/genomes:/genomes \
docker://cobirna/rnanue ./RNAnue/build/RNAnue complete --config /params/params_udorn.cfg

# Run RNAnue on wsn reads
singularity run --no-home \
-B /data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/reads:/trtms,\
/data/dessertlocal/projects/gl_iav-splash_freiburg/results/benchmarks:/outdir,\
/data/dessertlocal/projects/gl_iav-splash_freiburg/src/scripts_and_notebooks/benchmarks/rnanue:/params,\
/data/fass2/reads/gl_iav-splash_freiburg/dadonaite_2019/genomes:/genomes \
docker://cobirna/rnanue ./RNAnue/build/RNAnue complete --config /params/params_wsn.cfg
