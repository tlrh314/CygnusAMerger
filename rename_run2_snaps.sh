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

TIMESTAMP="20160617T1544"

cd "${BASEDIR}/runs/${TIMESTAMP}/snaps"

# Set the snapshot number to continue with 
# e.g. snapshot_011 is the last run1 snapshot, then j=12
j=12

# TODO: this magic works for leading 0's lol
# for i in 00{1..9} 0{10..99} {100..102}; do
# Set the number of snapshots produced by the restarted run
# e.g. is snapshot_run2_032 is the last snapshot, then use "for i in {0..32}"
for i in {0..32}; do
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

    # Use echo to check what is about to happen :-)..
    echo "snapshot_run2_0${snapnr}" "snapshot_${snapnr_new}"
    # TODO: hmm mv was not verbose in the last run...
    #mv "snapshot_run2_0${snapnr}" "snapshot_${snapnr_new}"
done
