#!/usr/bin/bash

# Script to inspect and rate QL analysis results.
# Usage: ./scripts/view_ql_results.sh <directory> [start_date]
# This script is expecting the <directory> to contain a 'waveforms' subdirectory,
# in which the SAC waveforms are stored with filenames formatted as 'YYYYMMDDHHMMSS_[NET].[STA].[CHA]'.
# There should be Z, T and R channels (i.e. already rotated).

# This looks at the waveforms in *reverse* chronological order, optionally starting from a given (upper bound) date

# This is written for a Mac, so it inspects the figures in Preview,
# and also uses osascript to bring the Terminal window to the front when rating the results.
# This is for speed purposes so you don't have to click. 
# Note that it will (try to) kill the Preview windows.

# Figures & results are saved in the 'results' subdirectory of the given directory.

# Ratings are saved in a logfile, with a rating of 0, 1 or 2, meant to mean
# 0 = no, 1 = yes, 2 = maybe (perhaps to have another look later).
# If the filename already appears in the logfile, it is skipped.

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
	if grep -q "$evstnm" "$logfile2"; then
		echo "done"
		continue
	fi
	if grep -q "$evstnm" "$logfile1"; then
		echo "doing"
	else
		continue
	fi
	echo $evstnm
	
	snr=$( awk 'NR==2 {print $10}' "$dir/results/$evstnm.results")
	
	if (( $(echo "$snr < 5" | bc -l) )); then
		echo "SNR = $snr too small"
		echo $evstnm 0 | awk '{printf "%s\t%s\n", $1, $2}' >> $logfile2
		continue
	fi
	
	# open -g "$dir/results/$evstnm.pdf"
	
	qlmanage -p "$dir/results/$evstnm.pdf" &>/dev/null &
	QL_PID=$!
	
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
	echo $evstnm $rating | awk '{printf "%s\t%s\n", $1, $2}' >> $logfile2
	
	cont=""
	echo "Continue? [0 to quit] "
	read -s -n1 cont
	if [ $cont == "0" ]; then
		echo "Okay, see you next time!"
		exit
	fi
		#
	# pid=$(pgrep -n Preview) # Try to close the Preview window we've opened
	# kill $pid
	kill $QL_PID
done