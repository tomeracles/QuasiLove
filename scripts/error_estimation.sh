#!/bin/bash

results=$1

fstem="${results%.*}"

out="${fstem}_err.txt"

echo "EVTIME EVLO EVLA EVDP MAG STLO STLA GCARC BAZ SNR G1 R1 DTQL DX SCLO SCLA AMPL H/V(R) H/V(QL) CORRAMP DXERR" > "$out"

while read line; do
	qldt=$(echo $line | awk '{print $13}')
	qldist=$(echo $line | awk '{print $14}')
	qlerror=$(echo "25 * $qldist / $qldt" | bc -l)
	qlerror=$(printf "%.1f" $qlerror)
	echo $qldt $qldist $qlerror
	# echo "$line $qlerror"
	echo "$line $qlerror" >> "$out"
done < $results

