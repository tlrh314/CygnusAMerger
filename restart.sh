#!/usr/bin/env bash
#PBS -lnodes=4:ppn=16:cores16
#PBS -lwalltime=03:00:00

setup_system() {
    # Better safe than sorry *_*... TODO: turn this off for production run?
    alias rm='rm -vi'
    alias cp='cp -vi'
    alias mv='mv -vi'
    SYSTYPE=`hostname`
    LOGLEVEL="DEBUG"

    if [[ "${SYSTYPE}" == *".lisa.surfsara.nl" ]]; then
        # TODO: check how multiple threas/nodes works on Lisa?
        # TODO: is the PBS situation the scheduler that also sets nodes/threads?
        # THREADS=$(grep -c ^processor /proc/cpuinfo)
        THREADS=64
        NICE=0  # default is 0
        BASEDIR="$HOME"  # TODO: look into the faster disk situation @Lisa?
        # TODO: I think home should not be used, instead use scratch??
        module load c/intel
        module load fftw2/sp/intel
        module load fftw2/dp/intel
        module load openmpi/intel
    elif [ "${SYSTYPE}" == "taurus" ]; then
        THREADS=8
        NICE=19
        BASEDIR="/scratch/timo"
    elif [ "${SYSTYPE}" == "Darwin" ]; then
        THREADS=$(sysctl -n hw.ncpu)  # OSX *_*
        NICE=0
        BASEDIR=""
        echo "Not implemented for OS X. Exiting."
        exit 1
    else
        echo "Unknown system. Exiting."
        exit 1
    fi

    # TODO match regex '/(?!-lnodes=)(?:\d*\.)?\d+/'
    # NODES=$(head $0 | grep "(?!-lnodes=)(?:\d*\.)?\d+")
    NODES=1
    GITHUBDIR="${BASEDIR}/CygnusAMerger"
    DATADIR="${BASEDIR}/runs"

    if [ "${LOGLEVEL}" == "DEBUG" ]; then
        echo "System settings"
        echo "Using machine   : ${SYSTYPE}"
        echo "Base directory  : ${BASEDIR}"
        echo "Github directory: ${GITHUBDIR}"
        echo "Niceness        : ${NICE}"
        echo "# MPI Nodes     : ${NODES}"
        echo -e "# OpenMP threads: ${THREADS}\n"
    fi

}


setup_system

cd /scratch/timo/runs/20160526T1354/snaps

# If restart files are present: only change parameterfile entries with a * in users-guide
# mpiexec Gadget2 gadget2.par 1

# If restart files are not present: restart from last snapshot and change Tbegin
nice -n 19 mpiexec.hydra -np $THREADS ./Gadget2 restart.par
