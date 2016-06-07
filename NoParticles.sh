#!/usr/bin/env bash

set -e

if [ -f NoParticles.txt ]; then
    rm NoParticles.txt
fi

# Make sure Toycluster is compiled to exit after setup.
# Edit Toycluster's main.c to do so.
# Otherwise the sph regularisation ensures it'll take forever... :-)
# OMP_NUM_THREADS=4 ../Toycluster/Toycluster_SetupOnly ToyclusterTrial.par > NoParticles.txt
OMP_NUM_THREADS=4 ../Toycluster/Toycluster_SetupOnly_delta200 ToyclusterTrial.par > NoParticles.txt
