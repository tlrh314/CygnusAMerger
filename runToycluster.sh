#!/bin/bash
#PBS -lnodes=1
#PBS -lwalltime=00:30:00

# load the necessary module(s), for example:
# module load openmpi/gnu

set -e

# Set up logging of the output
TIMESTAMP=$(date +"%Y%m%dT%H%M")
LOGFILE="runToycluster.log"

#exec > >(tee "${LOGFILE}")
# Pipe stderr to the logfile too.
#exec 2>&1

echo "Start of program at `date`"
echo "Logging to: ${LOGFILE}"


# Settings for machine
SYSTYPE=`hostname`
echo "Using machine ${SYSTYPE}."
# Set correct SYSTYPE in Makefile. Now detected in Makefile based on hostname!
#SYSTYPE="Taurus"
#SYSTYPE="Lisa"
#echo "Running on $SYSTYPE"
#if [ "$SYSTYPE" == "Taurus" ]; then
#    # Uncomment Taurus; Comment out Lisa
#    perl -pi -e 's/.*SYSTYPE="Lisa"/#SYSTYPE="Lisa"/g' Makefile
#    perl -pi -e 's/.*SYSTYPE="Taurus"/SYSTYPE="Taurus"/g' Makefile
#elif [ "$SYSTYPE" == "Lisa" ]; then
#    # Uncomment Lisa; comment out Taurus
#    perl -pi -e 's/.*SYSTYPE="Lisa"/SYSTYPE="Lisa"/g' Makefile
#    perl -pi -e 's/.*SYSTYPE="Taurus"/#SYSTYPE="Taurus"/g' Makefile
#fi

# Set the OMP number of threads
# On OSX $(sysctl -n hw.ncpu)
THREADS=$(grep -c ^processor /proc/cpuinfo)
NICE=0  # default is 0
BASEDIR="$HOME/Code"

if [ "${SYSTYPE}" == "taurus" ]; then
    echo "Taurus is also used by others. Maybe not use all resources :-)..."
    echo "  Maximum threads = ${THREADS}, but we will use 16!"
    THREADS=16
    NICE=19
    echo "  Running with nice -n ${NICE}"
    BASEDIR="/scratch/timo/Code"
fi
echo -e "OMP_NUM_THREADS = $OMP_NUM_THREADS \n"


# Set correct working directory
cd "${BASEDIR}/Toycluster"

# Make sure we use the latest makefile from the Git repository
cp "$BASEDIR/CygnusAMerger/Makefile_Toycluster" Makefile

PARAMETERFILE="cluster.par"
# Also use the parameterfile from the Git repository.
cp "$BASEDIR/CygnusAMerger/${PARAMETERFILE}" "${PARAMETERFILE}"

OUTDIR="../workdir/ToyclusterICs/${TIMESTAMP}"
if [ ! -d "${OUTDIR}" ]; then
    mkdir "${OUTDIR}"
fi

OUT=$(grep "Output_file" "${PARAMETERFILE}" |  tr ' ' '\n' | head -2 | tail -1)

# Compile the code
nice -n $NICE make clean
nice -n $NICE make -j8

# Run the code
OMP_NUM_THREADS=$THREADS nice --adjustment=$NICE ./Toycluster "${PARAMETERFILE}" >> "${LOGFILE}"

# Copy output to the output directory
mv Makefile "${OUTDIR}"
mv cluster.par "${OUTDIR}"
mv "${OUT}" "${OUTDIR}"
mv "${LOGFILE}" "${OUTDIR}"

echo "End of program at `date`"
