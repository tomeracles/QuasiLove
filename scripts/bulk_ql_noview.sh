#!/usr/bin/bash
# Script to run the QuasiLove detection on a set of waveforms.
# Recommended to run from parent directory of the QuasiLove directory.
# Usage: ./scripts/bulk_ql_noview.sh <directory> [start_date]
# This script is expecting the <directory> to contain a 'waveforms' subdirectory,
# in which the SAC waveforms are stored with filenames formatted as 'YYYYMMDDHHMMSS_[NET].[STA].[CHA]'.
# There should be Z, T and R channels (i.e. already rotated).

# This goes through the waveforms in *reverse* chronological order, optionally starting from a given (upper bound) date

# Figures & results are saved in the 'results' subdirectory of the given directory, to be inspected later (e.g. with the view_ql_results.sh script.)

# After analysis is done, each filename is saved in the logfile in the given directory (<directory>/logfile.txt, not <directory>/results/logfile.txt).
# If the filename already appears in the logfile, it is skipped.

source ./QuasiLove/quasilove_fns.sh # This is where the qlmain function is defined. If you're not in the parent directory, change this path accordingly.

dir=$1
logfile1="$1/logfile.txt"
# logfile2="${dir}/results/logfile.txt"

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
		# skip file if already in logfile
		echo "done"
		continue
	fi
	echo $evstnm

	# Do analysis
	qlmain $fstem "$dir/results"

	echo $evstnm >> $logfile1
done