#!/bin/bash

# ls -1 input/Resample1per10/*.wav > fwavin.list

export TERMGUI=mate-terminal

export INDIR='input/konvolusi7'
export OUTDIR='output/konvolusi7'

mkdir -p $OUTDIR

for i in `cat dest_angle.txt`;do
    WAVIN1=$(ls $INDIR/${i}_Track_1_*)
    WAVIN2=$(ls $INDIR/${i}_Track_2_*)
    WAVIN3=$(ls $INDIR/${i}_Track_3_*)
    WAVIN4=$(ls $INDIR/${i}_Track_4_*)

    CSVNAME=$OUTDIR/hasil_rata_${i}.csv

    #echo "file $WAVIN1, $WAVIN2, $WAVIN3, $WAVIN4ls"

    $TERMGUI -e "matlab -nodisplay -nosplash -nodesktop -r estangle(\'$WAVIN1\',\'$WAVIN2\',\'$WAVIN3\',\'$WAVIN4\',\'$CSVNAME\')"
done
