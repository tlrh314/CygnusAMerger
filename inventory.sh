#!/usr/bin/env bash

# File: ./inventory.sh
# Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
# Date created: Thu May 19, 2016 10:31 am
# Last modified: Fri May 27, 2016 08:34 pm
#
# Description: Make table of all run simulation parameters

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
    INVENTORY="${DATADIR}/inventory.csv"

    if [ -f "${INVENTORY}" ]; then
        echo "Removing ${INVENTORY}..."
        rm "${INVENTORY}"
    fi

    if [ "${LOGLEVEL}" == "DEBUG" ]; then
        echo "System settings"
        echo "Using machine   : ${SYSTYPE}"
        echo "Base directory  : ${BASEDIR}"
        echo "Github directory: ${GITHUBDIR}"
        echo "Sim output dir  : ${DATADIR}"
        echo "Inventory csv   : ${INVENTORY}"
        echo
    fi
}

get_toycluster_parm() {
    if [ ! -f "${TOYCLUSTERPARAMETERS}" ]; then
        echo "Error: Toycluster parameterfile does not exist!"
    fi
    PARM=$(grep "$1" "${TOYCLUSTERPARAMETERS}" | head -n 1 | awk '{print $2}')
    #echo "    $1: ${PARM}"
    echo "${PARM}"
}
get_gadget_parm() {
    if [ ! -f "${GADGETPARAMETERS}" ]; then
        echo "Error: Gadget parameterfile does not exist!"
    fi
    PARM=$(grep "$1" "${GADGETPARAMETERS}" | head -n 1 | awk '{print $2}')
    #echo "    $1: ${PARM}"
    echo "${PARM}"
}


setup_system

echo -n "SimID, " >> "${INVENTORY}"
echo -n "Ntotal, Mtotal, Mass_Ratio, ImpactParam, Redshift, c_nfw_0, v_com_0, rc_0, c_nfw_1, v_com_1, rc_1, " >> "${INVENTORY}"
echo "TimeBegin, TimeMax, TimeBetSnapshot, BoxSize" >> "${INVENTORY}"

for RUN in "${DATADIR}"/*
do
    if [ -d "${RUN}" ]; then
        DIRNAME="${RUN##*/}"
        TOYCLUSTERPARAMETERS="${RUN}/ICs/toycluster.par"
        GADGETPARAMETERS="${RUN}/snaps/gadget2.par"

        echo "${DIRNAME}"
        echo -n "${DIRNAME}, " >> "${INVENTORY}"

        if [ -d "${RUN}/ICs" ]; then
            echo "  ICs generated   : yes"
            Ntotal=$(get_toycluster_parm "Ntotal")
            echo "    Ntotal          : $Ntotal"
            Mtotal=$(get_toycluster_parm "Mtotal")
            echo "    Mtotal          : $Mtotal"
            Mass_Ratio=$(get_toycluster_parm "Mass_Ratio")
            echo "    Mass_Ratio      : $Mass_Ratio"
            ImpactParam=$(get_toycluster_parm "ImpactParam")
            echo "    ImpactParam     : $ImpactParam"
            #Cuspy
            Redshift=$(get_toycluster_parm "Redshift")
            echo "    Redshift        : $Redshift"
            #Bfld_Norm
            #Bfld_Eta
            #Bfld_Scale
            #bf
            #h_100
            #UnitLength_in_cm
            #UnitMass_in_g
            #UnitVelocity_in_cm_per_s
            c_nfw_0=$(get_toycluster_parm "c_nfw_0")
            echo "    c_nfw_0         : $c_nfw_0"
            v_com_0=$(get_toycluster_parm "v_com_0")
            echo "    v_com_0         : $v_com_0"
            rc_0=$(get_toycluster_parm "rc_0")
            echo "    rc_0            : $rc_0"
            c_nfw_1=$(get_toycluster_parm "c_nfw_1")
            echo "    c_nfw_1         : $c_nfw_1"
            v_com_1=$(get_toycluster_parm "v_com_1")
            echo "    v_com_1         : $v_com_1"
            rc_1=$(get_toycluster_parm "rc_1")
            echo "    rc_1            : $rc_1"

            echo -n "${Ntotal}, ${Mtotal}, ${Mass_Ratio}, ${ImpactParam}, ${Redshift}, ${c_nfw_0}, ${v_com_0}, ${rc_0}, ${c_nfw_1}, ${v_com_1}, ${rc_1}, " >> "${INVENTORY}"
        else
            echo "  ICs generated   : no"
            echo -n ", , , , , , , , , , , " >> "${INVENTORY}"
        fi

        if [ -d "${RUN}/snaps" ]; then
            echo "  snaps generated : yes"
            TimeBegin=$(get_gadget_parm "TimeBegin")
            echo "    TimeBegin       : $TimeBegin"
            TimeMax=$(get_gadget_parm "TimeMax")
            echo "    TimeMax         : $TimeMax"
            TimeBetSnapshot=$(get_gadget_parm "TimeBetSnapshot")
            echo "    TimeBetSnapshot : $TimeBetSnapshot"
            BoxSize=$(get_gadget_parm "BoxSize")
            echo "    BoxSize         : $BoxSize"
            echo "${TimeBegin}, ${TimeMax}, ${TimeBetSnapshot}, ${BoxSize}" >> "${INVENTORY}"
        else
            echo "  snaps generated : no"
            echo ", , , " >> "${INVENTORY}"
        fi

        if [ -d "${RUN}/analysis" ]; then
            echo "  smac generated  : yes"
        else
            echo "  smac generated  : no"
        fi
        echo
    fi
done

open -a /Applications/Microsoft\ Office\ 2008/Microsoft\ Excel.app/ "${INVENTORY}"
