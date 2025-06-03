#!/bin/bash

# Re-runs the quasilove calculations for all waveforms listed in the logfile.
# This is (presumably) to get the new (format of?) results after an update in the quasilove code.

dir=$1

logfile="${dir}/results/logfile.txt"

source codes/quasilove_fns.sh

while read line; do
    fstem=$(echo $line | awk '{print $1}')
    fstem="${dir}/waveforms/$fstem"
    # echo $fstem
    qlmain_no_plotting $fstem "$dir/results"
done < $logfile