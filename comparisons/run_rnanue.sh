#!/bin/sh

# Run RNAnue
singularity pull docker://cobirna/rnanue
singularity run --no-home \
    -B /data/dessertlocal/projects/gl_iav-splash_freiburg/data:/data,/data/dessertlocal/projects/gl_iav-splash_freiburg/results:/results \
    docker://cobirna/rnanue

# Run SPLASH scripts


