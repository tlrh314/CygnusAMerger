#!/usr/bin/env bash

set -e

if [ -f NoParticles.txt ]; then
    rm NoParticles.txt
fi

OMP_NUM_THREADS=4 ../Toycluster/Toycluster ToyclusterTrial.par > NoParticles.txt
