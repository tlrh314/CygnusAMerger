#!/usr/bin/env bash

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

cd "${BASEDIR}/runs/20160526T1354/snaps"

j=25

# TODO: this magic works for leading 0's lol
# for i in 00{1..9} 0{10..99} {100..102}; do
for i in {0..77}; do
    if [[ ${#i} -lt 2 ]]; then
        snapnr="0${i}"
    else
        snapnr="${i}"
    fi

    if [[ ${#j} -lt 2 ]]; then
        snapnr_new="00${j}"
    elif [[ ${#j} -eq 3 ]]; then
        snapnr_new="${j}"
    else
        snapnr_new="0${j}"
    fi
    j=$((j+1))

    echo "snapshot_run2_0${snapnr}" "snapshot_${snapnr_new}"
    #mv "snapshot_run2_0${snapnr}" "snapshot_${snapnr_new}"
done
