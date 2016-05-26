#!/usr/bin/env bash

# File: inventory.sh
# Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
# Date created: Thu May 19, 2016 10:31 am
# Last modified: Thu May 19, 2016 11:04 am
#
# Description: Make table of all run simulation parameters

set -e

LOGLEVEL="DEBUG"

setup_system() {
    # Better safe than sorry *_*... TODO: turn this off for production run?
    alias rm='rm -vi'
    alias cp='cp -vi'
    alias mv='mv -vi'
    SYSTYPE=`hostname`

    if [[ "${SYSTYPE}" == *".lisa.surfsara.nl" ]]; then
        BASEDIR="$HOME"
    elif [ "${SYSTYPE}" == "taurus" ]; then
        BASEDIR="/scratch/timo"
    elif [ "$(uname -s)" == "Darwin" ]; then
        SYSTYPE="MBP"
        BASEDIR="/Users/timohalbesma/Documents/Educatie/UvA/Master of Science Astronomy and Astrophysics/Jaar 3 (20152016)/Masterproject MScProj/Code"
    else
        echo "Unknown system. Exiting."
        exit 1
    fi

    GITHUBDIR="${BASEDIR}/CygnusAMerger"
    DATADIR="${BASEDIR}/runs"

    if [ "${LOGLEVEL}" == "DEBUG" ]; then
        echo "System settings"
        echo "Using machine   : ${SYSTYPE}"
        echo "Base directory  : ${BASEDIR}"
        echo "Github directory: ${GITHUBDIR}"
        echo "Sim output dir  : ${DATADIR}"
        echo
    fi

}

setup_system

for RUN in "${DATADIR}"/*
do
    if [ -d "${RUN}" ]; then
        DIRNAME="${RUN##*/}"
        TOYCLUSTERPARAMETERS="${RUN}/ICs/toycluster.par"

        echo "${DIRNAME}"

        if [ -d "${RUN}/ICs" ]; then
            echo "  ICs generated   : yes"
            # TODO: use ioparser.parse_toycluster_parms ?
        else
            echo "  ICs generated   : no"
        fi

        if [ -d "${RUN}/snaps" ]; then
            echo "  snaps generated : yes"
        else
            echo "  snaps generated : no"
        fi

        if [ -d "${RUN}/analysis" ]; then
            echo "  smac generated  : yes"
        else
            echo "  smac generated  : no"
        fi
        echo
        break
    fi
done
