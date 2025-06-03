#!/usr/bin/bash

source ./QuasiLove/quasilove_fns.sh

if [ -z "$1" ]; then
	echo "Usage: $0 <SAC file>"
	exit 1
fi

inf=$1 # One of the SAC files, doesn't matter which

infdir=$(dirname "$inf")
parentdir=$(dirname "$infdir")
outdir="${parentdir}/results"
mkdir -p ${outdir}


fstem=${inf%.*}
evstnm=${fstem##*/}

echo $evstnm
# ./codes/quasilove.sh $fstem "$dir/results"

# sac <<-EOF 
# r $inf
# lh
# q
# EOF

qlmain $fstem "$outdir"

# open "$dir/results/$evstnm.pdf"




