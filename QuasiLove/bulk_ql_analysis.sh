#!/bin/bash

dir=$1
logfile1="$1/logfile.txt"
logfile2="${dir}/results/logfile.txt"

if [ -z $2 ]; then
	start='1900'
else
	start=$2
fi

start="$start"0000000000
echo start $start

for f in $(echo ${dir}/waveforms/??????????????_*.??Z | sed $'s/ /\\\n/g' | sort -r); do
	fstem=${f%.*}
	evnm=${fstem##*/}
	date=$(echo "$evnm" | cut -f 1 -d'_')
	if (( $date > $start )); then
		continue
	fi
	if grep -q "$evnm" "$logfile1"; then
		echo "done"
		continue
	fi
	echo $evnm
	./codes/quasilove.sh $fstem "$dir/results"
	
	open "$dir/results/$evnm.pdf"
	
	# Here's a weird command to try to bring the Terminal window to the front. It seems to work [sometimes]
	osascript -e 'tell application "Terminal" to activate'
	
	# Now decide how good the split is
	rating=""
	echo "Accept? 1 = yes, 0 = no, 2 = ...maybe? "
	read -s -n1 rating
	
	until [ $rating == "0" ] || [ $rating == "1" ] || [ $rating == "2" ]; do
		echo "Frightfully sorry, you must choose 0, 1 or 2. Do try again: "
		read rating
	done
	echo $evnm >> $logfile1
	echo $evnm $rating | awk '{printf "%s\t%s\n", $1, $2}' >> $logfile2
	
	cont=""
	echo "Continue? [0 to quit] "
	read -s -n1 cont
	if [ $cont == "0" ]; then
		echo "Okay, see you next time!"
		exit
	fi
	
	pid=$(pgrep -n Preview) # Try to close the Preview window we've opened
	kill $pid
done