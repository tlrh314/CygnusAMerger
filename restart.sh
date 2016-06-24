#!/usr/bin/env bash
#PBS -lnodes=8:ppn=16:cores16
#PBS -l walltime=01:23:00:00

send_mail() {
    if [[ "${SYSTYPE}" == *".lisa.surfsara.nl" ]]; then
        msg="Dear Timo,\n\n\
I'm done with PBS JobID: ${PBS_JOBID}\n\n\
Succesfully ran job @ ${SYSTYPE}.\n\n\
Cheers,\n${SYSTYPE}"
        SUBJECT="Job @ ${SYSTYPE} is done executing :-)!"
        (echo -e $msg | mail $USER -s "${SUBJECT}") && echo "Mail Sent."
    elif [ "${SYSTYPE}" == "taurus" ]; then
        msg="To: timohalbesma@gmail.com\nFrom: tlrh@${SYSTYPE}\n\
Subject: ${0} @ ${SYSTYPE} is done executing :-)!\n\n\
Dear Timo,\n\n\
\"${0} $@\" is now done executing.\n\n\
Cheers,\n${SYSTYPE}"
        (echo -e $msg | sendmail -t timohalbesma@gmail.com) && echo "Mail Sent."
    fi
}

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
        THREADS=128  # set based on the nr of nodes requested
        NICE=0  # default is 0
        BASEDIR="$HOME"  # TODO: look into the faster disk situation @Lisa?
        MAIL=true
        # TODO: I think home should not be used, instead use scratch??
        module load c/intel
        module load fftw2/sp/intel
        module load fftw2/dp/intel
        module load openmpi/intel

        # Send mail once jobs starts.. it could queue god knows how long
        msg="Dear Timo,\n\n\
Gadget-2 run started for ${TIMESTAMP} with PBS JobID = ${PBS_JOBID}.\n\n\
The job runs at @ ${SYSTYPE}.\n\n\
Cheers,\n${SYSTYPE}"

        SUBJECT="Job @ ${SYSTYPE} has started :-)!"
    
        (echo -e $msg | mail $USER -s "${SUBJECT}") 
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


# TO ADJUST FOR RESTARTING
TIMESTAMP="20160617T1535"  # qsub no parse options..
setup_system

cp restart.par ${BASEDIR}/runs/${TIMESTAMP}/snaps
cd ${BASEDIR}/runs/${TIMESTAMP}/snaps

# If restart files are present: only change parameterfile entries with a * in users-guide
# mpiexec Gadget2 gadget2.par 1


# Make sure you copy restart.par to appropriate directory and adjust the
# required settings!

# If restart files are not present: restart from last snapshot and change Tbegin
if [[ "${SYSTYPE}" == *".lisa.surfsara.nl" ]]; then
    mpiexec ./Gadget2 restart.par
elif [ "${SYSTYPE}" == "taurus" ]; then
    nice -n ${NICE} mpiexec.hydra -np $THREADS ./Gadget2 restart.par
fi

if [ "$MAIL" = true ]; then
    send_mail
fi
