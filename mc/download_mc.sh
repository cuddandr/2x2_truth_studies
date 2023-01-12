#!/bin/bash

FILEPREFIX="https://portal.nersc.gov/project/dune/data/2x2/simulation/edepsim/NuMI_FHC_CHERRY/Merged2x2MINERvA_noRock_NuMI_FHC_CHERRY_5E17"
FILESUFFIX="EDEPSIM.root"
LOG_FILE="download_mc.log"
DOWNLOAD_FILE=""

if [ -z ${1} ]; then
    echo "Need start file number"
    exit 121
fi

if [ -z ${2} ]; then
    echo "Need end file number"
    exit 121
fi

START=${1}
NUMFILES=${2}
for (( c=$START; c<=$NUMFILES; c++ )); do
    printf -v DOWNLOAD_FILE "${FILEPREFIX}_%03d_${FILESUFFIX}\n" "$c"
    echo "Downloading $DOWNLOAD_FILE"

    wget -nv -a $LOG_FILE $DOWNLOAD_FILE
done
