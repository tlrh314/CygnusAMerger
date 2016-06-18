#!/usr/bin/env bash
#PBS -lnodes=8:ppn=16:cores16
#PBS -l walltime=01:23:00:00

# File: run.sh
# Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
# Date created: Wed Apr 27, 2016 06:40 PM
#
# Description: run simulation pipeline

set -e

LOGLEVEL="DEBUG"

_usage() {
cat <<EOF
Usage: `basename $0` <options>
Compile&Run Toycluster, Compile&Run Gadget-2, Compile&Run P-Smac2

Make sure GITHUBDIR contains:
    - Makefile_Toycluster, toycluster.par
    - Makefile_Gadget2, gadget2.par
    - Makefile_PSmac2, Config_PSmac2, smac2.par
New simulations will use these files for the setup.

Simulations have a 'simulation ID' yyyymmddThhmm. This is used as output dir.
Each simulation run will create three directories.
    - 'ICs' contains the Toycluster Makefile/parameters, executable and ICs
       as F77 unformatted binary file together with a txt file with Toycluster
       runtime output.
    - 'snaps' contains the Gadget2 Makefile/parameters, executable and snapshots
       as F77 unformatted binary files together with all other Gadget2 runtime
       output (cpu/energy/timings/info text files and restart files).
    - 'analysis' contains the P-Smac2 Makefile/Config/parameters, executable
       and any observable generated from the snapshots in the 'snaps' dir.

Running additional P-Smac2 routines requires providing 'simulation ID' as
parameter to the runscript. Leaving it blank results in a new simulation.

Options:
  -e   --effect         select which observable to calculate using P-Smac2
                        valid options: "rho", "xray", "Tem", "Tspec", "Pth",
                                       "SZy", "SZdt", "DMrho", "DMan", "all"
  -h   --help           display this help and exit
  -l   --loglevel       set loglevel to 'ERROR', 'WARNING', 'INFO', 'DEBUG'
                        TODO: to implement loglevel (do I want this?)
  -m   --mail           send email when code execution is completed
  -r   --restart        resume Gadget-2 simulation (TODO: to implement!)
  -t   --timestamp      provide a timestamp/'simulation ID'. Do not run
                        new simulation but set paths for this specific simulation.

Examples:
  `basename $0` --timestamp="20160502T1553" --mail
  `basename $0` -mt "20160502T1553"
  `basename $0` -m -t "20160502T1553"
  `basename $0` -m -t "20160518T0348" -e "xray"
EOF
}


parse_options() {
    # https://stackoverflow.com/questions/192249
    # It is possible to use multiple arguments for a long option.
    # Specifiy here how many are expected.
    declare -A longoptspec
    longoptspec=( [loglevel]=1 [timestamp]=1 [effect]=1 )

    optspec=":e:hlmrt:-:"
    while getopts "$optspec" opt; do
    while true; do
        case "${opt}" in
            -) # OPTARG is long-option or long-option=value.
                # Single argument:   --key=value.
                if [[ "${OPTARG}" =~ .*=.* ]]
                then
                    opt=${OPTARG/=*/}
                    OPTARG=${OPTARG#*=}
                    ((OPTIND--))
                # Multiple arguments: --key value1 value2.
                else
                    opt="$OPTARG"
                    OPTARG=(${@:OPTIND:$((longoptspec[$opt]))})
                fi
                ((OPTIND+=longoptspec[$opt]))
                # opt/OPTARG set, thus, we can process them as if getopts would've given us long options
                continue
                ;;
            e|effect)
                EFFECT="${OPTARG}"
                echo "EFFECT          = ${EFFECT}"
                ;;
            h|help)
                _usage
                exit 2  # 2 means incorrect usage
                ;;
            l|loglevel)
                echo "This function is not implemented"
                loglevel="${OPTARG}"
                echo "The loglevel is $loglevel"
                exit 1
                ;;
            m|mail)
                MAIL=true
                echo "MAIL            = true"
                ;;
            r|restart)
                echo "This function is not implemented"
                echo "OPTARG=${OPTARG}$"
                exit 1
                ;;
            t|timestamp)
                TIMESTAMP="${OPTARG}"
                echo "Timestamp       = ${TIMESTAMP}"
                ;;
        esac
    break; done
    done

    # Not sure if this is needed...
    # shift $((OPTIND-1))
}

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

    if [[ "${SYSTYPE}" == *".lisa.surfsara.nl" ]]; then
        # TODO: check how multiple threas/nodes works on Lisa?
        # TODO: is the PBS situation the scheduler that also sets nodes/threads?
        # THREADS=$(grep -c ^processor /proc/cpuinfo)
        THREADS=128  # set based on the nr of nodes requested
        NICE=0  # default is 0
        BASEDIR="$HOME"  # TODO: look into the faster disk situation @Lisa?
        TIMESTAMP="20160617T1544"  # qsub no parse options..
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
    elif [ "$(uname -s)" == "Darwin" ]; then
        SYSTYPE="MBP"
        BASEDIR="/Users/timohalbesma/Documents/Educatie/UvA/Master of Science Astronomy and Astrophysics/Jaar 3 (20152016)/Masterproject MScProj/Code"
        THREADS=$(sysctl -n hw.ncpu)  # OSX *_*
        NICE=10
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
    cp "${TOYCLUSTERMAKEFILE_GIT}" "${TOYCLUSTERDIR}/Makefile"
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
    mv "${TOYCLUSTERDIR}/Makefile" "${TOYCLUSTERMAKEFILE}"

    TOYCLUSTEREXECNAME=$(grep "EXEC =" "${TOYCLUSTERMAKEFILE}" | cut -d' ' -f3)
    mv "${TOYCLUSTERDIR}/${TOYCLUSTEREXECNAME}" "${ICOUTDIR}"
    TOYCLUSTEREXEC="${ICOUTDIR}/${TOYCLUSTEREXECNAME}"

    TOYCLUSTERPARAMETERS_GIT="${GITHUBDIR}/toycluster.par"
    if [ ! -f "${TOYCLUSTERPARAMETERS_GIT}" ]; then
        echo "Error: ${TOYCLUSTERPARAMETERS_GIT} does not exist!"
        exit 1
    fi
    cp "${TOYCLUSTERPARAMETERS_GIT}" "${ICOUTDIR}"
    TOYCLUSTERPARAMETERS="${ICOUTDIR}/toycluster.par"
    TOYCLUSTERLOGFILE="${ICOUTDIR}/${TOYCLUSTERLOGFILENAME}"
}

run_toycluster() {
    echo "Running Toycluster..."
    cd "${ICOUTDIR}"

    SECONDS=0
    OMP_NUM_THREADS=$THREADS nice -n $NICE "${TOYCLUSTEREXEC}" "${TOYCLUSTERPARAMETERS}" 2>&1 >> "${TOYCLUSTERLOGFILE}"
    RUNTIME=$SECONDS
    HOUR=$(($RUNTIME/3600))
    MINS=$(( ($RUNTIME%3600) / 60))
    SECS=$(( ($RUNTIME%60) ))
    printf "Runtime = %d s, which is %02d:%02d:%02d\n" "$RUNTIME" "$HOUR" "$MINS" "$SECS"

    echo "... done running Toycluster"
}

check_toycluster_run() {
    # TODO: not all of these parameters exists. If they dont exist they cant be printed :-)...
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
    # And of course osx requires stat -f '%m %N' to do the exact same thing
    # sort: -n --> numerical, -k1 --> key = 1 (==last modification)
    # so we do not let glob do its thing but we get a sorted list of snapshots
    if [ "${SYSTYPE}" == "MBP" ]; then
        cd "${SIMOUTDIR}"
        GADGETSNAPSHOTS=($(stat -f '%m %N' snapshot_*\
            | sort -t ' ' -nk1 | cut -d ' ' -f2-))  # outer parenthesis --> array
        cd -
    else

        GADGETSNAPSHOTS=($(stat -c '%Y %n' "${SIMOUTDIR}"/snapshot_*  \
            | sort -t ' ' -nk1 | cut -d ' ' -f2-))  # outer parenthesis --> array
    fi

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
    cp "${GADGETMAKEFILE_GIT}" "${GADGETDIR}/Makefile"
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
    mv "${GADGETDIR}/Makefile" "${GADGETMAKEFILE}"
    GADGETEXECNAME=$(grep "EXEC = " "${GADGETMAKEFILE}" | cut -d' ' -f3)
    mv "${GADGETDIR}/${GADGETEXECNAME}" "${SIMOUTDIR}"
    GADGETEXEC="${SIMOUTDIR}/${GADGETEXECNAME}"

    GADGETPARAMETERS_GIT="${GITHUBDIR}/gadget2.par"
    if [ ! -f "${GADGETPARAMETERS_GIT}" ]; then
        echo "Error: ${GADGETPARAMETERS_GIT} does not exist!"
        exit 1
    fi
    cp "${GADGETPARAMETERS_GIT}" "${SIMOUTDIR}"
    GADGETPARAMETERS="${SIMOUTDIR}/gadget2.par"

    # Set correct BoxSize
    # echo "Press enter to continue..." && read enterKey
    BOXSIZE=$(grep "Boxsize" "${TOYCLUSTERLOGFILE}" | cut -d'=' -f2 | cut -d' ' -f2)
    echo "Setting BoxSize in Gadget parameter file to: ${BOXSIZE}"
    perl -pi -e 's/BoxSize.*/BoxSize '${BOXSIZE}' % kpc/g' "${GADGETPARAMETERS}"
    grep -n --color=auto "BoxSize" "${GADGETPARAMETERS}"
    # echo "Press enter to continue..." && read enterKey

    if [ ! -f "${ICFILE}" ]; then
        echo "Error: ${ICFILE} does not exist!"
        exit 1
    fi
    GADGETICFILE="${SIMOUTDIR}/${ICFILENAME}"
    cp "${ICFILE}" "${GADGETICFILE}"
    #GADGETLOGFILE="${SIMOUTDIR}/${GADGETLOGFILENAME}"
}

run_gadget() {
    echo "Running Gadget2..."
    cd "${SIMOUTDIR}"

    if [[ "${SYSTYPE}" == *".lisa.surfsara.nl" ]]; then
        # the pee-bee-es situation fixes nodes/threads?
        mpiexec "${GADGETEXEC}" "${GADGETPARAMETERS}" # 2&1> "${GADGETLOGFILE}"
    elif [ "${SYSTYPE}" == "taurus" ]; then
        nice -n $NICE mpiexec.hydra -np $THREADS "${GADGETEXEC}" "${GADGETPARAMETERS}" # 2&1> "${GADGETLOGFILE}"
    fi

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
        echo "Error: ${GADGETCPU} does not exist!"
        exit 1
    fi
}

setup_psmac2() {
    if [[ "${SYSTYPE}" == *".lisa.surfsara.nl" ]]; then
        return
    fi
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
        cp "${PSMAC2MAKEFILE_GIT}" "${PSMAC2DIR}/Makefile"
        PSMAC2CONFIG_GIT="${GITHUBDIR}/Config_PSmac2"
        if [ ! -f "${PSMAC2CONFIG_GIT}" ]; then
            echo "Error: ${PSMAC2CONFIG_GIT} does not exist!"
            exit 1
        fi
        cp "${PSMAC2CONFIG_GIT}" "${PSMAC2DIR}/Config"

        compile_psmac2
        PSMAC2MAKEFILE="${ANALYSISDIR}/Makefile_PSmac2"
        mv "${PSMAC2DIR}/Makefile" "${PSMAC2MAKEFILE}"
        PSMAC2CONFIG="${ANALYSISDIR}/Config_PSmac2"
        mv "${PSMAC2DIR}/Config" "${PSMAC2CONFIG}"

        PSMAC2EXECNAME=$(grep "EXEC =" "${PSMAC2MAKEFILE}" | cut -d' ' -f3)
        mv "${PSMAC2DIR}/${PSMAC2EXECNAME}" "${ANALYSISDIR}"
        PSMAC2EXEC="${ANALYSISDIR}/${PSMAC2EXECNAME}"
    else
        PSMAC2MAKEFILE="${ANALYSISDIR}/Makefile_PSmac2"
        PSMAC2CONFIG="${ANALYSISDIR}/Config_PSmac2"
        PSMAC2EXECNAME=$(grep "EXEC =" "${PSMAC2MAKEFILE}" | cut -d' ' -f3)
        PSMAC2EXEC="${ANALYSISDIR}/${PSMAC2EXECNAME}"
    fi
    set_psmac2_generic_runtime_files

    check_psmac2_generic_runtime_files

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

    set_psmac_parameterfile_snapshot_path

}

compile_psmac2() {
    echo "Compiling P-Smac2"
    cd "${PSMAC2DIR}"

    # Compile the code
    nice -n $NICE make clean
    nice -n $NICE make

    echo "... done compiling P-Smac2"
}

set_psmac2_generic_runtime_files() {
    PSMAC2PARAMETERSNAME="smac2.par"
    PSMAC2PARAMETERS_GIT="${GITHUBDIR}/${PSMAC2PARAMETERSNAME}"
    if [ ! -f "${PSMAC2PARAMETERS_GIT}" ]; then
        echo "Error: ${PSMAC2PARAMETERS_GIT} does not exist!"
        exit 1
    fi
    cp "${PSMAC2PARAMETERS_GIT}" "${ANALYSISDIR}"
    PSMAC2PARAMETERS="${ANALYSISDIR}/${PSMAC2PARAMETERSNAME}"
    PSMAC2LOGFILE="${ANALYSISDIR}/${PSMAC2LOGFILENAME}"
}

check_psmac2_generic_runtime_files() {
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
    if [ ! -f "${PSMAC2PARAMETERS}" ]; then
        echo "Error: ${PSMAC2PARAMETERS} does not exist!"
        exit 1
    fi
}

set_psmac_parameterfile_snapshot_path() {
    # Use ${#GADGETSNAPSHOTS}, which is the nr of snapshots?
    # SNAPMAX=${#GADGETSNAPSHOTS[@]}  # length of array
    # SNAPMAX=$(printf "%03d" $(( 10#$SNAPMAX-1 )))  # 10# -> force decimal
    # echo "Setting Input_File to: snapshot_000 snapshot_${SNAPMAX} 1"
    # Match line containing Input_File; set fits output name

    if [ -z "${GADGETSNAPSHOTS}" ]; then
        echo "Warning: no Gadget snapshots. Did you run the simulation?"
        echo "Assuming we want to run P-Smac2 with ICs!"
        echo -e "\nSetting Input_File to: /path/to/IC_file"
        ICFILE_ESCAPED=$(echo "${ICFILE}" | sed -e 's/[]\/$*.^|[]/\\&/g')
        perl -pi -e "s/Input_File.*/Input_File ${ICFILE_ESCAPED}/g" "${PSMAC2PARAMETERS}"
        grep -n --color=auto "Input_File" "${PSMAC2PARAMETERS}"
        return
    else
        FIRST="${GADGETSNAPSHOTS[0]}"  # Globbed, then sorted array
        LAST="${GADGETSNAPSHOTS[-1]}"

        if [ "${SYSTYPE}" == "MBP" ]; then
            FIRST="../snaps/${FIRST}"
            LAST="../snaps/${LAST}"
        fi

        # Escape forward slashes. If / not escaped we break perl :)
        FIRST=$(echo "${FIRST}" | sed -e 's/[]\/$*.^|[]/\\&/g')
        LAST=$(echo "${LAST}" | sed -e 's/[]\/$*.^|[]/\\&/g')

        echo -e "\nSetting Input_File to: /path/to/snapshot_000 /path/to/snapshot_max 1"
        perl -pi -e "s/Input_File.*/Input_File ${FIRST} ${LAST} 1/g" "${PSMAC2PARAMETERS}"
        #For IC only
        #perl -pi -e "s/Input_File.*/Input_File ${FIRST}/g" "${PSMAC2PARAMETERS}"
        grep -n --color=auto "Input_File" "${PSMAC2PARAMETERS}"
        echo
    fi
}

run_psmac2_for_given_module() {
    check_psmac2_generic_runtime_files
    for PROJECTION in x y z; do
        SMAC_PREFIX="${1}_projection-${PROJECTION}"
        EFFECT_MODULE="${2}"
        EFFECT_FLAG="${3}"

        PSMAC2PARAMETERS="${ANALYSISDIR}/${SMAC_PREFIX}_${PSMAC2PARAMETERSNAME}"
        PSMAC2LOGFILE="${ANALYSISDIR}/${SMAC_PREFIX}_${PSMAC2LOGFILENAME}"
        cp "${ANALYSISDIR}/${PSMAC2PARAMETERSNAME}" "${PSMAC2PARAMETERS}"

        if [ "${PROJECTION}" == "x" ]; then
            echo "  Running x-projection: [deg] along x : (0,90,90)"
            Euler_Angle_0="0"
            Euler_Angle_1="90"
            Euler_Angle_2="90"
        elif [ "${PROJECTION}" == "y" ]; then
            echo "   Running y-projection: [deg] along y : (0,90,0)"
            Euler_Angle_0="0"
            Euler_Angle_1="90"
            Euler_Angle_2="0"
        elif [ "${PROJECTION}" == "z" ]; then
            echo  "  Running z-projection: [deg] along z : (0,0,0)"
            Euler_Angle_0="0"
            Euler_Angle_1="0"
            Euler_Angle_2="0"
        else
            echo "Error: unknown projection angle!"
        fi

        # echo "    Setting Euler_Angle_0 to: ${Euler_Angle_0}"
        perl -pi -e 's/^Euler_Angle_0.*/Euler_Angle_0 '${Euler_Angle_0}'/g' "${PSMAC2PARAMETERS}"

        # echo "    Setting Euler_Angle_1 to: ${Euler_Angle_1}"
        perl -pi -e 's/^Euler_Angle_1.*/Euler_Angle_1 '${Euler_Angle_1}'/g' "${PSMAC2PARAMETERS}"

        # echo "    Setting Euler_Angle_2 to: ${Euler_Angle_2}"
        perl -pi -e 's/^Euler_Angle_2.*/Euler_Angle_2 '${Euler_Angle_2}'/g' "${PSMAC2PARAMETERS}"

        # Set smac2.par to run DM, change outputfile and logfile
        # Match line containing Output_File; set fits output name
        OUTPUTFILE="${SMAC_PREFIX}.fits"
        # echo "    Setting Output_File to: ${OUTPUTFILE}"
        perl -pi -e 's/Output_File.*/Output_File '${OUTPUTFILE}'/g' "${PSMAC2PARAMETERS}"

        # echo "    Setting Effect_Module to: ${EFFECT_MODULE}"
        perl -pi -e 's/Effect_Module.*/Effect_Module '${EFFECT_MODULE}'/g' "${PSMAC2PARAMETERS}"


        # echo "    Setting Effect_Flag to: ${EFFECT_FLAG}"
        perl -pi -e 's/Effect_Flag.*/Effect_Flag '${EFFECT_FLAG}'/g' "${PSMAC2PARAMETERS}"

        echo
        echo "    Analysis dir    : ${ANALYSISDIR}"
        echo "    Parameterfile   : ${PSMAC2PARAMETERS}"
        echo "    Logging to      : ${PSMAC2LOGFILE}"
        echo "    Output fits file: ${OUTPUTFILE}"
        echo "    Effect_Module   : ${EFFECT_MODULE}"
        echo "    Effect_Flag     : ${EFFECT_FLAG}"
        echo "    Projection      : ${PROJECTION}"
        echo "    Euler_Angle_0   : ${Euler_Angle_0}"
        echo "    Euler_Angle_1   : ${Euler_Angle_1}"
        echo "    Euler_Angle_2   : ${Euler_Angle_2}"

        echo
        echo "    Checking ${PSMAC2PARAMETERS}"
        echo -n "      "
        grep -n --color=auto "Output_File" "${PSMAC2PARAMETERS}"
        echo -n "      "
        grep -n --color=auto "Effect_Module" "${PSMAC2PARAMETERS}"
        echo -n "      "
        grep -n --color=auto "Effect_Flag" "${PSMAC2PARAMETERS}"
        echo -n "      "
        grep -n --color=auto "^Euler_Angle_0" "${PSMAC2PARAMETERS}"
        echo -n "      "
        grep -n --color=auto "^Euler_Angle_1" "${PSMAC2PARAMETERS}"
        echo -n "      "
        grep -n --color=auto "^Euler_Angle_2" "${PSMAC2PARAMETERS}"
        echo

        if [ ! -f "${PSMAC2PARAMETERS}" ]; then
            echo "Error: ${PSMAC2PARAMETERS} does not exist!"
            exit 1
        fi

        echo "    Running P-Smac2..."
        cd "${ANALYSISDIR}"
        SECONDS=0
        if [[ "${SYSTYPE}" == *".lisa.surfsara.nl" ]]; then
            # the pee-bee-es situation fixes nodes/threads?
            mpiexec "${PSMAC2EXEC}" "${PSMAC2PARAMETERS}" >> "${PSMAC2LOGFILE}" 2>&1
        elif [ "${SYSTYPE}" == "taurus" ]; then
            OMP_NUM_THREADS=$THREADS nice --adjustment=$NICE mpiexec.hydra -np $NODES "${PSMAC2EXEC}" "${PSMAC2PARAMETERS}" >> "${PSMAC2LOGFILE}" 2>&1 
        elif [ "${SYSTYPE}" == "MBP" ]; then
            # EXEC="./${PSMAC2EXEC##*/}"
            # PARM="${PSMAC2PARAMETERS##*/}"
            OMP_NUM_THREADS=$THREADS nice -n $NICE mpiexec -np $NODES "${PSMAC2EXEC}" "${PSMAC2PARAMETERS}" >> "${PSMAC2LOGFILE}" 2>&1
        fi
        RUNTIME=$SECONDS
        HOUR=$(($RUNTIME/3600))
        MINS=$(( ($RUNTIME%3600) / 60))
        SECS=$(( ($RUNTIME%60) ))
        printf "      Runtime = %d s, which is %02d:%02d:%02d\n" "$RUNTIME" "$HOUR" "$MINS" "$SECS"

        echo "    ... done running P-Smac2"
    done
}


# Main
# Uncomment if options are required
# if [ $# = 0 ]; then _usage && exit 2; fi

if [[ ! "${SYSTYPE}" == *".lisa.surfsara.nl" ]]; then
    parse_options $@
fi

echo -e "\nStart of program at $(date)\n"

setup_system
setup_toycluster
setup_gadget
# echo "Press enter to continue" && read enterKey
#setup_psmac2

# TODO: if directories do not exist the parameter is empty and echoing it makes no sense...

case "${EFFECT}" in
    "rho")
        echo "Running P-Smac2 for physical density."
        #  0 - Physical Density [g/cm^2]; no Flag
        run_psmac2_for_given_module "physical-density" "0" "0"
        exit 0
        ;;
    "xray")
        # 2 - X-Ray Surface Brightness; no Flag
        echo "Running P-Smac2 for X-Ray Surface Brightness."
        run_psmac2_for_given_module "xray-surface-brightness" "2" "0"
        exit 0
        ;;
    "Tem")
        # 4 - Temperature
        #     2 - Emission Weighted
        echo "Running P-Smac2 for Emission Weighted Temperature."
        run_psmac2_for_given_module "temperature-emission-weighted" "4" "2"
        exit 0
        ;;
    "Tspec")
        # 4 - Temperature
        #     3 - Spectroscopic - Chandra, XMM (Mazotta+ 04)
        echo "Running P-Smac2 for Spectroscopic Temperature (Chandra)."
        run_psmac2_for_given_module "temperature-spectroscopic" "4" "3"
        exit 0
        ;;
    "Pth")
        # 5 - Pressure
        #     0 - Thermal
        echo "Running P-Smac2 for Thermal Pressure."
        run_psmac2_for_given_module "presssure" "5" "0"
        exit 0
        ;;
    "SZy")
        # 7 - SZ Effect
        #     0 - Compton-y (=Smac1 thermal DT/T)
        echo "Running P-Smac2 for Sunyaev-Sel'dovic effect: Compton-y parameter."
        run_psmac2_for_given_module "SZ-Compton-y" "7" "0"
        exit 0
        ;;
    "SZdt")
        # 7 - SZ Effect
        #     1 - Thermal SZ DT/T
        echo "Running P-Smac2 for Sunyaev-Sel'dovic effect: Thermal, DT/T"
        run_psmac2_for_given_module "sz-thermal-dt-over-t" "7" "1"
        exit 0
        ;;
    "DMrho")
        echo "Running P-Smac2 for Dark Matter density."
        # 10 - DM Density; no Flag
        run_psmac2_for_given_module "dm-density" "10" "0"
        exit 0
        ;;
    "DMan")
        # 11 - DM Annihilation Signal (rho^2); no Flag
        echo "Running P-Smac2 for Dark Matter Annihilation."
        run_psmac2_for_given_module "dm-annihilation" "11" "0"
        exit 0
        ;;
    "all")
        echo "Running P-Smac2 for physical density."
        run_psmac2_for_given_module "physical-density" "0" "0"

        echo "Running P-Smac2 for X-Ray Surface Brightness."
        run_psmac2_for_given_module "xray-surface-brightness" "2" "0"

        echo "Running P-Smac2 for Emission Weighted Temperature."
        run_psmac2_for_given_module "temperature-emission-weighted" "4" "2"

        echo "Running P-Smac2 for Spectroscopic Temperature (Chandra)."
        run_psmac2_for_given_module "temperature-spectroscopic" "4" "3"

        echo "Running P-Smac2 for Thermal Pressure."
        run_psmac2_for_given_module "presssure" "5" "0"

        echo "Running P-Smac2 for Sunyaev-Sel'dovic effect: Compton-y parameter."
        run_psmac2_for_given_module "sz-compton-y" "7" "0"

        echo "Running P-Smac2 for Sunyaev-Sel'dovic effect: Thermal, DT/T"
        run_psmac2_for_given_module "sz-thermal-dt-over-t" "7" "1"

        echo "Running P-Smac2 for Dark Matter density."
        #run_psmac2_for_given_module "dm-density" "10" "0"

        echo "Running P-Smac2 for Dark Matter Annihilation."
        run_psmac2_for_given_module "dm-annihilation" "11" "0"
esac

if [ "$MAIL" = true ]; then
    send_mail
fi

echo -e "\nEnd of program at $(date)\n"
