#!/bin/bash

source ./codes/quasilove_fns.sh

dir=$1
logfile1="$1/logfile.txt"
logfile2="${dir}/results/logfile.txt"

if [ -z $2 ]; then
	start='2100'
else
	start=$2
fi

start="$start"0000000000
echo start $start

for f in $(echo ${dir}/waveforms/??????????????_*.??Z | sed $'s/ /\\\n/g' | sort -r); do
	fstem=${f%.*}
	evstnm=${fstem##*/}
	date=$(echo "$evstnm" | cut -f 1 -d'_')
	if (( $date > $start )); then
		continue
	fi
	if grep -q "$evstnm" "$logfile1"; then
		echo "done"
		continue
	fi
	echo $evstnm
	# ./codes/quasilove.sh $fstem "$dir/results"
	qlmain $fstem "$dir/results"

	echo $evstnm >> $logfile1
done