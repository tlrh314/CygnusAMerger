#!/usr/bin/env bash
#PBS -lnodes=1
#PBS -lwalltime=00:30:00

# File: run.sh
# Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
# Date created: Wed Apr 27, 2016 06:40 PM
# Last modified: Thu Apr 28, 2016 03:05 PM
#
# Description: run simulation pipeline

set -e

LOGLEVEL="DEBUG"

# TODO: set up for Lisa
send_mail() {
    msg="To: timohalbesma@gmail.com\nFrom: tlrh@${SYSTYPE}\n\
Subject: ${0} @ ${SYSTYPE} is done executing :-)!\n\n\
Dear Timo,\n\n\
\"${0} $@\" is now done executing.\n\n\
Cheers,\n${SYSTYPE}"
    (echo -e $msg | sendmail -t timohalbesma@gmail.com) && echo "Mail Sent."
}

setup_system() {
    SYSTYPE=`hostname`

    if [ "${SYSTYPE}" == "*.lisa.surfsara.nl" ]; then
        THREADS=$(grep -c ^processor /proc/cpuinfo)
        NICE=0  # default is 0
        BASEDIR="$HOME"
    elif [ "${SYSTYPE}" == "taurus" ]; then
        THREADS=4
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

    NODES=$(head $0 | grep "#PBS -lnodes" | cut -d'=' -f2)
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

setup_toycluster() {
    TOYCLUSTERDIR="${BASEDIR}/Toycluster"
    TOYCLUSTERLOGFILENAME="runToycluster.log"

    # If TIMESTAMP is not given we assume no initial conditions exist!
    if [ -z $TIMESTAMP ]; then
        TIMESTAMP=$(date +"%Y%m%dT%H%M")
        echo "No timestamp  --> generating new ICs"
        echo "Timestamp       : ${TIMESTAMP}"

        SIMULATIONDIR="${DATADIR}/${TIMESTAMP}"
        ICOUTDIR="${SIMULATIONDIR}/ICs"

        if [ ! -d "${SIMULATIONDIR}" ]; then
            mkdir "${SIMULATIONDIR}"
            echo "Created         : ${SIMULATIONDIR}"
        fi

        if [ ! -d "${ICOUTDIR}" ]; then
            mkdir "${ICOUTDIR}"
            echo "Created         : ${ICOUTDIR}"
        fi

        set_toycluster_compile_files
        # echo "Press enter to continue" && read enterKey
        compile_toycluster
        set_toycluster_runtime_files
        run_toycluster

        ICFILENAME=$(grep "Output_file" "${TOYCLUSTERPARAMETERS}" | cut -d' ' -f2)
        ICFILE="${ICOUTDIR}/${ICFILENAME:2}"

    else
        SIMULATIONDIR="${DATADIR}/${TIMESTAMP}"
        ICOUTDIR="${SIMULATIONDIR}/ICs"

        TOYCLUSTERMAKEFILE="${ICOUTDIR}/Makefile_Toycluster"
        TOYCLUSTEREXECNAME=$(grep "EXEC =" "${TOYCLUSTERMAKEFILE}" | cut -d' ' -f3)
        TOYCLUSTEREXEC="${ICOUTDIR}/${TOYCLUSTEREXECNAME}"
        TOYCLUSTERPARAMETERS="${ICOUTDIR}/toycluster.par"
        TOYCLUSTERLOGFILE="${ICOUTDIR}/${TOYCLUSTERLOGFILENAME}"
        ICFILENAME=$(grep "Output_file" "${TOYCLUSTERPARAMETERS}" | cut -d' ' -f2)
        ICFILE="${ICOUTDIR}/${ICFILENAME:2}"
    fi

    check_toycluster_run

    if [ "${LOGLEVEL}" == "DEBUG" ]; then
        echo -e "\nInitial Condition paths"
        echo "Toycluster dir  : ${TOYCLUSTERDIR}"
        echo "Simulation 'ID' : ${TIMESTAMP}"
        echo "Simulation dir  : ${SIMULATIONDIR}"
        echo "IC directory    : ${ICOUTDIR}"
        echo "Makefile        : ${TOYCLUSTERMAKEFILE}"
        echo "Toycluster exec : ${TOYCLUSTEREXEC}"
        echo "Parameterfile   : ${TOYCLUSTERPARAMETERS}"
        echo "Logfile         : ${TOYCLUSTERLOGFILE}"
        echo "ICs file        : ${ICFILE}"
        echo
    fi
}

set_toycluster_compile_files() {
    TOYCLUSTERMAKEFILE_GIT="${GITHUBDIR}/Makefile_Toycluster"
    if [ ! -f "${TOYCLUSTERMAKEFILE_GIT}" ]; then
        echo "Error: ${TOYCLUSTERMAKEFILE_GIT} does not exist!"
        exit 1
    fi
    cp -i "${TOYCLUSTERMAKEFILE_GIT}" "${TOYCLUSTERDIR}/Makefile"
}

compile_toycluster() {
    echo "Compiling Toycluster..."
    cd "${TOYCLUSTERDIR}"

    # Compile the code
    nice -n $NICE make clean
    nice -n $NICE make -j8

    echo "... done compiling Toycluster"
}

set_toycluster_runtime_files() {
    TOYCLUSTERMAKEFILE="${ICOUTDIR}/Makefile_Toycluster"
    mv -i "${TOYCLUSTERDIR}/Makefile" "${TOYCLUSTERMAKEFILE}"

    TOYCLUSTEREXECNAME=$(grep "EXEC =" "${TOYCLUSTERMAKEFILE}" | cut -d' ' -f3)
    mv -i "${TOYCLUSTERDIR}/${TOYCLUSTEREXECNAME}" "${ICOUTDIR}"
    TOYCLUSTEREXEC="${ICOUTDIR}/${TOYCLUSTEREXECNAME}"

    TOYCLUSTERPARAMETERS_GIT="${GITHUBDIR}/toycluster.par"
    if [ ! -f "${TOYCLUSTERPARAMETERS_GIT}" ]; then
        echo "Error: ${TOYCLUSTERPARAMETERS_GIT} does not exist!"
        exit 1
    fi
    cp -i "${TOYCLUSTERPARAMETERS_GIT}" "${ICOUTDIR}"
    TOYCLUSTERPARAMETERS="${ICOUTDIR}/toycluster.par"
    TOYCLUSTERLOGFILE="${ICOUTDIR}/${TOYCLUSTERLOGFILENAME}"
}

run_toycluster() {
    echo "Running Toycluster..."
    cd "${ICOUTDIR}"

    SECONDS=0
    OMP_NUM_THREADS=$THREADS nice --adjustment=$NICE "${TOYCLUSTEREXEC}" "${TOYCLUSTERPARAMETERS}" 2>&1 >> "${TOYCLUSTERLOGFILE}"
    RUNTIME=$SECONDS
    HOUR=$(($RUNTIME/3600))
    MINS=$(( ($RUNTIME%3600) / 60))
    SECS=$(( ($RUNTIME%60) ))
    printf "Runtime = %d s, which is %02d:%02d:%02d\n" "$RUNTIME" "$HOUR" "$MINS" "$SECS"

    echo "... done running Toycluster"

}

check_toycluster_run() {
    if [ ! -d "${SIMULATIONDIR}" ]; then
        echo "Error: ${SIMULATIONDIR} does not exist!"
        exit 1
    fi
    if [ ! -d "${ICOUTDIR}" ]; then
        echo "Error: ${ICOUTDIR} does not exist!"
        exit 1
    fi
    if [ ! -f "${TOYCLUSTERMAKEFILE}" ]; then
        echo "Error: ${TOYCLUSTERMAKEFILE} does not exist!"
        exit 1
    fi
    if [ ! -f "${TOYCLUSTEREXEC}" ]; then
        echo "Error: ${TOYCLUSTEREXEC} does not exist!"
        exit 1
    fi
    if [ ! -f "${TOYCLUSTERPARAMETERS}" ]; then
        echo "Error: ${TOYCLUSTERPARAMETERS} does not exist!"
        exit 1
    fi
    if [ ! -f "${TOYCLUSTERLOGFILE}" ]; then
        echo "Error: ${TOYCLUSTERLOGFILE} does not exist!"
        exit 1
    fi
    if [ ! -f "${ICFILE}" ]; then
        echo "Error: ${ICFILE} does not exist!"
        exit 1
    fi
}

setup_gadget() {
    GADGETDIR="${BASEDIR}/Gadget-2.0.7/Gadget2"
    #GADGETLOGFILENAME="runGadget.log"

    SIMOUTDIR="${SIMULATIONDIR}/snaps"

    if [ -z "${TIMESTAMP}" ]; then
        echo "Error: no timestamp given!"
        exit 1
    fi
    if [ ! -d "${SIMOUTDIR}" ]; then
        echo "No SIMOUTDIR  --> Running simulations!"
        mkdir "${SIMOUTDIR}"
        echo "Created         : ${SIMOUTDIR}"

        set_gadget_compile_files
        # echo "Press enter to continue" && read enterKey
        compile_gadget
        set_gadget_runtime_files
        run_gadget
    else
        GADGETMAKEFILE="${SIMOUTDIR}/Makefile_Gadget2"
        GADGETEXECNAME=$(grep "EXEC = " "${GADGETMAKEFILE}" | cut -d' ' -f3)
        GADGETEXEC="${SIMOUTDIR}/${GADGETEXECNAME}"
        GADGETPARAMETERS="${SIMOUTDIR}/gadget2.par"
        GADGETICFILE="${SIMOUTDIR}/${ICFILENAME}"
        #GADGETLOGFILE="${SIMOUTDIR}/${GADGETLOGFILENAME}"
    fi

    GADGETENERGY="${SIMOUTDIR}/energy.txt"
    GADGETINFO="${SIMOUTDIR}/info.txt"
    GADGETTIMINGS="${SIMOUTDIR}/timings.txt"
    GADGETCPU="${SIMOUTDIR}/cpu.txt"

    check_gadget_run

    # stat: -c --> format, where %Y --> last modification; %n --> filename
    # sort: -n --> numerical, -k1 --> key = 1 (==last modification)
    # so we do not let glob do its thing but we get a sorted list of snapshots
    GADGETSNAPSHOTS=($(stat -c '%Y %n' "${SIMOUTDIR}"/snapshot_*  \
        | sort -t ' ' -nk1 | cut -d ' ' -f2-))  # outer parenthesis --> array

    if [ "${LOGLEVEL}" == "DEBUG" ]; then
        echo -e "\nSimulation paths"
        echo "Gadget dir      : ${GADGETDIR}"
        echo "Simulation 'ID' : ${TIMESTAMP}"
        echo "Simulation dir  : ${SIMULATIONDIR}"
        echo "Snapshot dir    : ${SIMOUTDIR}"
        echo "Makefile        : ${GADGETMAKEFILE}"
        echo "Gadget exec     : ${GADGETEXEC}"
        echo "Parameterfile   : ${GADGETPARAMETERS}"
        # echo "Logfile         : ${GADGETLOGFILE}"
        echo "energy.txt      : ${GADGETENERGY}"
        echo "info.txt        : ${GADGETINFO}"
        echo "timings.txt     : ${GADGETTIMINGS}"
        echo "cpu.txt         : ${GADGETCPU}"

        for s in "${GADGETSNAPSHOTS[@]}"; do
            echo "snapshot        : ${s}"
        done
        echo
    fi
}

set_gadget_compile_files() {
    GADGETMAKEFILE_GIT="${GITHUBDIR}/Makefile_Gadget2"
    if [ ! -f "${GADGETMAKEFILE_GIT}" ]; then
        echo "Error: ${GADGETMAKEFILE_GIT} does not exist!"
        exit 1
    fi
    cp -i "${GADGETMAKEFILE_GIT}" "${GADGETDIR}/Makefile"
}

compile_gadget() {
    # If needed: unpack tar archive with source code
    # ( cd source ; tar xzv --strip-components=1 -f - ) < gadget-2.0.7.tar.gz

    echo "Compiling Gadget..."

    cd "${GADGETDIR}"
    nice -n $NICE make clean
    nice -n $NICE make -j8

    echo "... done compiling Gadget"
}

set_gadget_runtime_files() {
    GADGETMAKEFILE="${SIMOUTDIR}/Makefile_Gadget2"
    mv -i "${GADGETDIR}/Makefile" "${GADGETMAKEFILE}"
    GADGETEXECNAME=$(grep "EXEC = " "${GADGETMAKEFILE}" | cut -d' ' -f3)
    mv -i "${GADGETDIR}/${GADGETEXECNAME}" "${SIMOUTDIR}"
    GADGETEXEC="${SIMOUTDIR}/${GADGETEXECNAME}"

    GADGETPARAMETERS_GIT="${GITHUBDIR}/gadget2.par"
    if [ ! -f "${GADGETPARAMETERS_GIT}" ]; then
        echo "Error: ${GADGETPARAMETERS_GIT} does not exist!"
        exit 1
    fi
    cp -i "${GADGETPARAMETERS_GIT}" "${SIMOUTDIR}"
    GADGETPARAMETERS="${SIMOUTDIR}/gadget2.par"

    # Set correct BoxSize
    echo "Press enter to continue..." && read enterKey
    BOXSIZE=$(grep "Boxsize" "${TOYCLUSTERLOGFILE}" | cut -d'=' -f2 | cut -d' ' -f2)
    echo "Setting BoxSize in Gadget parameter file to: ${BOXSIZE}"
    perl -pi -e 's/BoxSize.*/BoxSize '${BOXSIZE}' % kpc/g' "${GADGETPARAMETERS}"
    grep -n --color=auto "BoxSize" "${GADGETPARAMETERS}"
    echo "Press enter to continue..." && read enterKey

    if [ ! -f "${ICFILE}" ]; then
        echo "Error: ${ICFILE} does not exist!"
        exit 1
    fi
    GADGETICFILE="${SIMOUTDIR}/${ICFILENAME}"
    cp -i "${ICFILE}" "${GADGETICFILE}"
    #GADGETLOGFILE="${SIMOUTDIR}/${GADGETLOGFILENAME}"
}

run_gadget() {
    echo "Running Gadget2..."
    cd "${SIMOUTDIR}"

    OMP_NUM_THREADS=$THREADS nice --adjustment=$NICE mpiexec.hydra -np $NODES "${GADGETEXEC}" "${GADGETPARAMETERS}" # 2&1> "${GADGETLOGFILE}"

    echo "... done running Gadget2"
}

check_gadget_run() {
    if [ ! -d "${SIMULATIONDIR}" ]; then
        echo "Error: ${SIMULATIONDIR} does not exist!"
        exit 1
    fi
    if [ ! -d "${SIMOUTDIR}" ]; then
        echo "Error: ${SIMOUTDIR} does not exist!"
        exit 1
    fi
    if [ ! -f "${GADGETMAKEFILE}" ]; then
        echo "Error: ${GADGETMAKEFILE} does not exist!"
        exit 1
    fi
    if [ ! -f "${GADGETEXEC}" ]; then
        echo "Error: ${GADGETEXEC} does not exist!"
        exit 1
    fi
    if [ ! -f "${GADGETPARAMETERS}" ]; then
        echo "Error: ${GADGETPARAMETERS} does not exist!"
        exit 1
    fi
    if [ ! -f "${GADGETICFILE}" ]; then
        echo "Error: ${GADGETICFILE} does not exist!"
        exit 1
    fi
    # if [ ! -f "${GADGETLOGFILE}" ]; then
    #     echo "Error: ${GADGETLOGFILE} does not exist!"
    #     exit 1
    # fi
    if [ ! -f "${GADGETENERGY}" ]; then
        echo "Error: ${GADGETENERGY} does not exist!"
        exit 1
    fi
    if [ ! -f "${GADGETINFO}" ]; then
        echo "Error: ${GADGETINFO} does not exist!"
        exit 1
    fi
    if [ ! -f "${GADGETTIMINGS}" ]; then
        echo "Error: ${GADGETTIMINGS} does not exist!"
        exit 1
    fi
    if [ ! -f "${GADGETCPU}" ]; then
        echo "Error: ${GADCPU} does not exist!"
        exit 1
    fi
}

setup_psmac2() {
    PSMAC2DIR="${BASEDIR}/P-Smac2"
    PSMAC2LOGFILENAME="runPSmac2.log"

    ANALYSISDIR="${SIMULATIONDIR}/analysis"

    if [ ! -d "${ANALYSISDIR}" ]; then  # compile :)
        echo "No ANALYSISDIR--> Compiling P-Smac2!"
        mkdir "${ANALYSISDIR}"
        echo "Created         : ${ANALYSISDIR}"

        PSMAC2MAKEFILE_GIT="${GITHUBDIR}/Makefile_PSmac2"
        if [ ! -f "${PSMAC2MAKEFILE_GIT}" ]; then
            echo "Error: ${PSMAC2MAKEFILE_GIT} does not exist!"
            exit 1
        fi
        cp -i "${PSMAC2MAKEFILE_GIT}" "${PSMAC2DIR}/Makefile"
        PSMAC2CONFIG_GIT="${GITHUBDIR}/Config_PSmac2"
        if [ ! -f "${PSMAC2CONFIG_GIT}" ]; then
            echo "Error: ${PSMAC2CONFIG_GIT} does not exist!"
            exit 1
        fi
        cp -i "${PSMAC2CONFIG_GIT}" "${PSMAC2DIR}/Config"

        compile_psmac2
        PSMAC2MAKEFILE="${ANALYSISDIR}/Makefile_PSmac2"
        mv -i "${PSMAC2DIR}/Makefile" "${PSMAC2MAKEFILE}"
        PSMAC2CONFIG="${ANALYSISDIR}/Config_PSmac2"
        mv -i "${PSMAC2DIR}/Config" "${PSMAC2CONFIG}"

        PSMAC2EXECNAME=$(grep "EXEC =" "${PSMAC2MAKEFILE}" | cut -d' ' -f3)
        mv -i "${PSMAC2DIR}/${PSMAC2EXECNAME}" "${ANALYSISDIR}"
        PSMAC2EXEC="${ANALYSISDIR}/${PSMAC2EXECNAME}"
    else
        PSMAC2MAKEFILE="${ANALYSISDIR}/Makefile_PSmac2"
        PSMAC2CONFIG="${ANALYSISDIR}/Config_PSmac2"
        PSMAC2EXECNAME=$(grep "EXEC =" "${PSMAC2MAKEFILE}" | cut -d' ' -f3)
        PSMAC2EXEC="${ANALYSISDIR}/${PSMAC2EXECNAME}"
    fi

    if [ ! -d "${SIMULATIONDIR}" ]; then
        echo "Error: ${SIMULATIONDIR} does not exist!"
        exit 1
    fi
    if [ ! -d "${ANALYSISDIR}" ]; then
        echo "Error: ${ANALYSISDIR} does not exist!"
        exit 1
    fi
    if [ ! -f "${PSMAC2MAKEFILE}" ]; then
        echo "Error: ${PSMAC2MAKEFILE} does not exist!"
        exit 1
    fi
    if [ ! -f "${PSMAC2CONFIG}" ]; then
        echo "Error: ${PSMAC2CONFIG} does not exist!"
        exit 1
    fi

    # Run
    PSMAC2PARAMETERS_GIT="${GITHUBDIR}/smac2.par"
    if [ ! -f "${PSMAC2PARAMETERS_GIT}" ]; then
        echo "Error: ${PSMAC2PARAMETERS_GIT} does not exist!"
        exit 1
    fi
    cp -i "${PSMAC2PARAMETERS_GIT}" "${ANALYSISDIR}"
    PSMAC2PARAMETERS="${ANALYSISDIR}/smac2.par"
    PSMAC2LOGFILE="${ANALYSISDIR}/${PSMAC2LOGFILENAME}"

    # TODO: set inputfiles to the number of generated snapshots
    # Use ${#GADGETSNAPSHOTS}, which is the nr of snapshots?
    SNAPMAX=${#GADGETSNAPSHOTS[@]}  # length of array
    SNAPMAX=$(printf "%03d" $(( 10#$SNAPMAX-1 )))  # 10# -> force decimal
    echo "Setting Input_File to: snapshot_000 snapshot_${SNAPMAX} 1"
    # Match line containing Input_File; set fits output name
    perl -pi -e 's/Input_File.*/Input_File snapshot_000 snapshot_'${SNAPMAX}' 1/g' "${PSMAC2PARAMETERS}" 
    grep -n --color=auto "Input_File" "${PSMAC2PARAMETERS}"

    # Match line containing Output_File; set fits output name
    OUTPUTFILE="dm.fits"
    echo "Setting Output_File to: ${OUTPUTFILE}"
    perl -pi -e 's/Output_File.*/Output_File '${OUTPUTFILE}'/g' "${PSMAC2PARAMETERS}" 
    grep -n --color=auto "Output_File" "${PSMAC2PARAMETERS}"

    # Effect_Module
    # Effect_Flag




    #if [ ! -f "${PSMAC2LOGFILE}" ]; then
    #    echo "Error: ${PSMAC2LOGFILE} does not exist!"
    #    exit 1
    #fi

    if [ "${LOGLEVEL}" == "DEBUG" ]; then
        echo -e "\nAnalysis paths"
        echo "P-Smac2 dir     : ${PSMAC2DIR}"
        echo "Simulation 'ID' : ${TIMESTAMP}"
        echo "Simulation dir  : ${SIMULATIONDIR}"
        echo "Analysis dir    : ${ANALYSISDIR}"
        echo "Makefile        : ${PSMAC2MAKEFILE}"
        echo "Config          : ${PSMAC2CONFIG}"
        echo "P-Smac2 exec    : ${PSMAC2EXEC}"
        echo "Parameterfile   : ${PSMAC2PARAMETERS}"
        echo -e "Logfile         : ${PSMAC2LOGFILE}\n"
    fi
}

compile_psmac2() {
    echo "Compiling P-Smac2"
    cd "${PSMAC2DIR}"

    # Compile the code
    nice -n $NICE make clean
    nice -n $NICE make -j8

    echo "... done compiling P-Smac2"
}

run_psmac2() {
    echo "Running P-Smac2..."
    cd "${ANALYSISDIR}"

    SECONDS=0
    OMP_NUM_THREADS=$THREADS nice --adjustment=$NICE mpiexec.hydra -np $NODES "${PSMAC2EXEC}" "${PSMAC2PARAMETERS}" 2>&1 >> "${PSMAC2LOGFILE}"
    RUNTIME=$SECONDS
    HOUR=$(($RUNTIME/3600))
    MINS=$(( ($RUNTIME%3600) / 60))
    SECS=$(( ($RUNTIME%60) ))
    printf "Runtime = %d s, which is %02d:%02d:%02d\n" "$RUNTIME" "$HOUR" "$MINS" "$SECS"

    echo "... done running P-Smac2"
}

# Main
echo -e "\nStart of program at $(date)\n"

# TODO: now only works if all functions are called
#TIMESTAMP="20160428T1149"
TIMESTAMP="20160428T1459"
setup_system
setup_toycluster
setup_gadget
#echo "Press enter to continue" && read enterKey
setup_psmac2

echo -e "\nEnd of program at $(date)\n"
