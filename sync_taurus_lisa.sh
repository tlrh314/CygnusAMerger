#!/usr/bin/env bash
#
# File: sync_taurus_lisa.sh
# Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
# Date created: Sat Apr 23, 2016 11:09 PM
# Last modified:

set -e

SYSTYPE=`hostname`
echo "Using machine ${SYSTYPE}."

TAURUSDIR="/scratch/timo/Code/Simulations/"
LISADIR="/home/timoh/Code/Simulations/"


# Rsynch both ways :-)
# -a archive mode, implies -rlptgoD 
# -u update mode; skip files that are newer on the receiver
# -H preserve hard links
# -x don't cross filesystem boundaries
# -z compress file data during the transfer
# -r recurse into directories
# -l copy symlinks as symlinks
# -p preserve permissions
# -t preserve times
# -g preserve group
# -o preserve owner, (only if run as root?)
# -D transfer deviced and sockets (only if run as root?)

# On Taurus
if [ "${SYSTYPE}" == "taurus" ]; then
    NICE=19
    nice -n $NICE rsync -auHxz --progress "${TAURUSDIR}" lisa:"${LISADIR}"
    nice -n $NICE rsync -auHxz --progress lisa:"${LISADIR}" "${TAURUSDIR}"
elif [ "${SYSTYPE}" == "*.lisa.surfsara.nl" ]; then
    # On Lisa
    NICE=0  # default is 0
    nice -n $NICE rsync -auHxz --progress taurus:"${TAURUSDIR}" "${LISADIR}"
    nice -n $NICE rsync -auHxz --progress "${LISADIR}" taurus:"${TAURUSDIR}"
else
    echo "Unknown system: ${SYSTYPE}. Exiting"
    exit 1
fi

