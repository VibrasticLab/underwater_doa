#!/bin/bash

# ls -1 input/Resample1per10/*.wav > fwavin.list

export TERMGUI=mate-terminal

export INDIR='input/Resample9per10'
export OUTDIR='output/konvolusi9'
export FLIST=fwavin.list

for i in `cat list/$FLIST`;do
    FNAME=$(basename $i .wav)
    FWAVIN=$INDIR/${FNAME}.wav
    FWAVOUT=$OUTDIR/${FNAME}_konv.wav
    #echo "Processing $FWAVIN to $FWAVOUT"
    $TERMGUI -e "matlab -nodisplay -nosplash -nodesktop -r matcliconv(\'$FWAVIN\',\'$FWAVOUT\')"
done
