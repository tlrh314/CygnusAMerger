#!/usr/bin/env bash
#PBS -lnodes=1
#PBS -lwalltime=00:30:00

# load the necessary module(s), for example:
# module load openmpi/gnu
# we need fftw, gsl, icc, openmpi, etc
# TODO: load modules and link in Makefile

# File: runGadget2.sh
# Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
# Date created: Fri Dec 04, 2015 03:44 PM
# Last modified: Wed Apr 27, 2016 11:08 PM
#
# Description: Compile Gadget-2, run simulation, copy Gadget-2 makefile/param

# "A good strategy for doing this in practice is to make a copy
# of the whole simulation source code together with its makefile in the output
# directory of each simulation run, and then use this copy to compile the code
# and to run the simulation. The code and its settings become then a logical
# part of the output generated in the simulation directory."
# --> We do not copy the entire directory but do copy the makefile/param files.

set -e

loglevel='INFO'

# Checking for indianness should not be necessary?
test_indianness() {
    gcc checkIndianness.c -o checkIndianness
    ind="$(./checkIndianness)"

    if [ $ind -eq 1 ]; then
        echo "This machine is little endian."
    elif [ $ind -eq 0 ]; then
        echo "This machine is big endian."
    else
        echo "ERROR: Unknown if big endian or little endian."
        exit 1
    fi
}

_usage() {
cat <<EOF
Usage: `basename $0` <[options]> <[filename]>
Compile Gadget-2, execute FILE with Gadget-2

Options:
  -r   --restart        set Gadget-2 restart flag to continue execution
                        implies the Gadget-2 executable will not be regenerated
                        i.e. the source will not be build.
  -l   --loglevel       set loglevel to 'ERROR', 'WARNING', 'INFO', 'DEBUG'
  -h   --help           display this help and exit
  -t   --test           test indianness of machine
  -v   --version        output version information and exit

Examples:
  `basename $0` merge_clusters.c
  `basename $0` -r merge_clusters.c
EOF
}


parse_options() {
    # It is possible to use multiple arguments for a long option.
    # Specifiy here how many are expected.
    declare -A longoptspec
    longoptspec=( [loglevel]=1 )

    optspec=":rlhtv-:"
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
            l|loglevel)
                loglevel=$OPTARG
                echo "The loglevel is $loglevel"
                ;;
            r|restart)
                restart_gadget=$OPTARG
                echo "Restart flag set"
                ;;
            h|help)
                _usage
                exit 2  # 2 means incorrect usage
                ;;
            t|test)
                test_indianness
                exit 0
                ;;
            v|version)
                grep "# Version:" $0
                grep "# Last modified:" $0
                echo "Copyright (C) 2015 Timo L. R. Halbesma, BSc."
                exit 0
                ;;
        esac
    break; done
    done

    # Not sure if this is needed...
    # shift $((OPTIND-1))
}

# 'Main'

# SYSTEM setup
# -----------------------------------------------------------------------------
# Set up logging of the output
# TIMESTAMP=$(date +"%Y%m%dT%H%M")
# If use time timestamp dir set by runToycluster.sh we can have all simulations with those ICs in the same dir
# So TIMESTAMP is a 'unique' simulation number?
TIMESTAMP="20160424T0016"
LOGFILE="runGadget.log"

echo "Start of program at `date`"
echo "Logging to: ${LOGFILE}"

# Uncomment if options are required
# if [ $# = 0 ]; then _usage; fi
parse_options $@

SYSTYPE=`hostname`
echo "Using machine ${SYSTYPE}."

# Set the OMP number of threads
# On OSX $(sysctl -n hw.ncpu)
THREADS=$(grep -c ^processor /proc/cpuinfo)
NICE=0  # default is 0
BASEDIR="$HOME"

if [ "${SYSTYPE}" == "taurus" ]; then
    echo "Taurus is also used by others. Maybe not use all resources :-)..."
    echo "  Maximum threads = ${THREADS}, but we will use 4!"
    THREADS=4
    NICE=19
    echo "  Running with nice -n ${NICE}"
    BASEDIR="/scratch/timo"
fi
echo -e "OMP_NUM_THREADS = $THREADS \n"
# -----------------------------------------------------------------------------


# Directory setup
# -----------------------------------------------------------------------------
# Set correct directory paths
WORKDIR="${BASEDIR}/Simulations/${TIMESTAMP}"
GADGETDIR="${BASEDIR}/Gadget-2.0.7/Gadget2"

PARAMETERFILE="gadget2.par"
MAKEFILE="${BASEDIR}/CygnusAMerger/Makefile_Gadget2"

if [ ! -d "${WORKDIR}" ]; then
    echo "Initial conditions directory ${TIMESTAMP} does not exist :-("
    exit 1
fi

if [ -f "${LOGFILE}" ]; then
    echo "Warning: logfile already exists!"
    exit 1
    LOGFILE="runGadget.log_$(date +"%Y%m%dT%H%M")"
    echo "Actually logging to: ${LOGFILE}"
    # TODO: if we will run simulation also rename makefile and parameterfile?
    # rm "${LOGFILE}"
fi
# -----------------------------------------------------------------------------


# Compilation
# -----------------------------------------------------------------------------
# Compile the code in the GADGETDIR, but run simulation in WORKDIR
cd "${GADGETDIR}"
# Make sure we use the latest makefile from the Git repository
cp "${MAKEFILE}" "${GADGETDIR}/Makefile"

# If needed: unpack tar archive with source code
# ( cd source ; tar xzv --strip-components=1 -f - ) < gadget-2.0.7.tar.gz

# Compile :-)
make help
exit 0
nice -n $NICE make clean
nice -n $NICE make -j8

mv -i Makefile "${WORKDIR}/Makefile_Gadget2"  # move makefile to output dir
mv -i Gadget2 "${WORKDIR}"                    # move executable to output dir
# -----------------------------------------------------------------------------


# Code execution
# -----------------------------------------------------------------------------
cd "${WORKDIR}"
# Use the parameterfile from the Git repository.
cp "$BASEDIR/CygnusAMerger/${PARAMETERFILE}" "${WORKDIR}/${PARAMETERFILE}"

# Run the code
nice --adjustment=$NICE mpiexec.hydra -np $THREADS ./Gadget2 "${PARAMETERFILE}" # >> "${LOGFILE}"


sleep 2
echo -n "Done executing code;"
# Send email to my own gmail =). TODO: change this for Lisa. TODO: make optional
# -----------------------------------------------------------------------------
msg="To: timohalbesma@gmail.com\nFrom: tlrh@${SYSTYPE}\n\
Subject: ${0} @ ${SYSTYPE} is done executing :-)!\n\n\
Dear Timo,\n\n\
\"${0} $@\" is now done executing.\n\n\
Cheers,\n${SYSTYPE}"
(echo -e $msg | sendmail -t timohalbesma@gmail.com) && echo -n " sent mail;"
# -----------------------------------------------------------------------------
