#!/usr/bin/env bash

set -e

mv="mv -i"
cp="cp -i"
rm="rm -i"

BASEDIR="/Users/timohalbesma/Documents/Educatie/UvA/Master of Science Astronomy and Astrophysics/Jaar 3 (20152016)/Masterproject MScProj/Code"
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

    # echo "snapshot_run2_0${snapnr}" "snapshot_${snapnr_new}"
    mv "snapshot_run2_0${snapnr}" "snapshot_${snapnr_new}"
done
