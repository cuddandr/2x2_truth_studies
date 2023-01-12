#!/bin/bash

FILEPREFIX="kaon_2x2_5E17"
FILESUFFIX="EDEPSIM.root"
# OUTDIR="kaon_mc"
OUTDIR="all_kaons"

i=0
for IN_FILE in *.root; do
    printf -v OUT_FILE "${FILEPREFIX}_%03d_${FILESUFFIX}" "$i"

    if [ ! -f "$OUTDIR/$OUT_FILE" ]; then
        echo "Reading $IN_FILE to $OUTDIR/$OUT_FILE"
        python3 kaon_picker.py -o $OUTDIR/$OUT_FILE -i $IN_FILE
    fi
    
    i=$((i+1))
done
