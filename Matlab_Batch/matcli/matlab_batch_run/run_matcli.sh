#!/bin/bash

# ls -1 data_rekaman/*.wav > fwavin.list

export TERMGUI=mate-terminal
export DATADIR=data_rekaman

export FLIST=5p10_3.list
export OUTDIR='Resample/resample6per10'

export SAMPLOW=5000
export SAMPHIGH=10000

for i in `cat list/$FLIST`;do
    FNAME=$(basename $i .wav)
    FWAVIN=$DATADIR/${FNAME}.wav
    FWAVOUT=$OUTDIR/${FNAME}_re.wav
    $TERMGUI -e "matlab -nodisplay -nosplash -nodesktop -r matclirun(\'$FWAVIN\',\'$FWAVOUT\',\'$SAMPLOW\',\'$SAMPHIGH\')"
done
