#!/usr/bin/env bash
# TODO: figure out how to make required number of nodes variable
# TODO: figure out how to hosts file for mpiexec at Lisa?
#PBS -lnodes=1
#PBS -lwalltime=00:30:00

# load the necessary module(s), for example:
# module load openmpi/gnu
# we need CFITSIO library?
# TODO: load modules and link in Makefile

# File: runPSmac2.sh
# Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
# Date created: Tue Apr 26, 2016 05:31 pm
# Last modified: Wed Apr 27, 2016 09:11 PM
#
# Description: Compile P-Smac2, obtain observables

set -e

# SYSTEM setup
# -----------------------------------------------------------------------------
SIMULATION="20160423T2219"
LOGFILE="runPSmac2.log"

echo "Start of program at `date`"
echo "Logging to: ${LOGFILE}"

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
echo -e "THREADS = $THREADS \n"
# -----------------------------------------------------------------------------


# Directory setup
# -----------------------------------------------------------------------------
# Set correct directory paths
SIMULATIONDIR="${BASEDIR}/Simulations/${SIMULATION}"
PSMAC2DIR="${BASEDIR}/P-Smac2"

PARAMETERFILE="smac2.par"
MAKEFILE="${BASEDIR}/CygnusAMerger/Makefile_PSmac2"
CONFIGFILE="${BASEDIR}/CygnusAMerger/Config_PSmac2"
OUTDIR="${SIMULATIONDIR}/${TIMESTAMP}" 

if [ ! -d "${SIMULATIONDIR}" ]; then
    echo "Initial conditions directory ${TIMESTAMP} does not exist :-("
    exit 1
fi

if [ ! -d "${OUTDIR}" ]; then
    mkdir "${OUTDIR}"
fi
# -----------------------------------------------------------------------------


# Compilation
# -----------------------------------------------------------------------------
# Compile the code in the PSMAC2DIR, but run simulation in OUTDIR
cd "${PSMAC2DIR}"
# Make sure we use the latest makefile from the Git repository
cp "${MAKEFILE}" "${PSMAC2DIR}/Makefile"
cp "${CONFIGFILE}" "${PSMAC2DIR}/Config"

# If needed: unpack tar archive with source code
# ( cd source ; tar xzv --strip-components=1 -f - ) < gadget-2.0.7.tar.gz

# Compile :-)
nice -n $NICE make clean
nice -n $NICE make -j8

mv -i Makefile "${OUTDIR}/Makefile_PSmac2"   # move makefile to output dir
mv -i Config "${OUTDIR}/Config_PSmac2"       # move config to output dir
mv -i P-Smac2 "${OUTDIR}"                    # move executable to output dir
# -----------------------------------------------------------------------------


# # MPI daemon
# # -----------------------------------------------------------------------------
# # Kill and restart mpi daemon, or start it.
# if [ -f "${BASEDIR}/.mpd.pid" ]; then
#     echo -e "\nWarning, MPI daemon is already running!"
#     echo "Continuing requires to kill mpd, thus, could kill running code."
#     echo -n "Kill mpd and continue? [y/n]: "
#     read ans
#     case ${ans:=y} in [yY]*) ;;
#         *) echo -e "\nCan't ever be too cautious, smart move...exiting\n" && exit 1 ;;
#     esac
#     mpdexit localmpd && echo "Cleanly exited mpd."
#     OUT=$?
#     if [ ! $OUT -eq 0 ]; then  # mpdexit failed: nuke it; kill it with fire!
#         PID=$(cat "${BASEDIR}/.mpd.pid")
#         if ! kill $PID > /dev/null 2>&1; then
#             echo "Could not send SIGTERM to process $PID" >&2
#             echo "MPD was not running? But now it is :-)"
#         else
#             echo "Successfully send SIGTERM to process $PID" >&2
#             echo "MPD was killed, then started."
#         fi
#     fi
# else
#     echo "No mpi daemon pid file found. Starting it."
# fi
# # TODO: do we always have to restart the mpi daemon?
# nice -n $NICE mpd --daemon --ncpus=$THREADS --pid="${BASEDIR}/.mpd.pid" < /dev/null
# -----------------------------------------------------------------------------

# Code execution
# -----------------------------------------------------------------------------
cd "${OUTDIR}"
# Use the parameterfile from the Git repository.
cp "$BASEDIR/CygnusAMerger/${PARAMETERFILE}" "${OUTDIR}/${PARAMETERFILE}"

# Run the code
OMP_NUM_THREADS=$THREADS nice --adjustment=$NICE mpiexec.hydra -np 1 ./P-Smac2 "${PARAMETERFILE}" # >> "${LOGFILE}"
# -----------------------------------------------------------------------------


sleep 2
echo -n "Done executing code;"

rm -f "${BASEDIR}/.mpd.pid" && echo -n " lock removed;"
mpdexit localmpd && echo " cleanly exited mpd."
exit 0
# Send email to my own gmail =). TODO: change this for Lisa. TODO: make option
# -----------------------------------------------------------------------------
msg="To: timohalbesma@gmail.com\nFrom: tlrh@${SYSTYPE}\n\
Subject: ${0} @ ${SYSTYPE} is done executing :-)!\n\n\
Dear Timo,\n\n\
\"${0} $@\" is now done executing.\n\n\
Cheers,\n${SYSTYPE}"
(echo -e $msg | sendmail -t timohalbesma@gmail.com) && echo "Mail Sent."
# -----------------------------------------------------------------------------
