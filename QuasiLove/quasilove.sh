#!/bin/bash
# * Look for potential Quasi-Love waves
# * USAGE: quasilove.sh filestem outdir 
# * where the components are expected to
# * be in the file format [filestem].?H[R,T,Z]
# *
# * waveforms should be filtered with a fmax of 0.01 Hz
# *
# * Love and Rayleigh waves are picked based on windows covering
# * approx max and min group velocities at periods >100s based on
# * model GDM52 (Ekstrom 2011); i.e.,
# * 3.8 km/s < V(G) < 4.7 km/s
# * 3.4 km/s < V(R) < 4.1 km/s
# *

sac=/usr/local/LINUX_sac/bin/sacinit.sh

fstem="$1"
outdir="$2"

if [ -z $outdir ]; then
	echo "Expected quasilove.sh [filestem] [outdir]"
	echo "Try again"
	exit
else
	echo "Running..."
fi

def_variables() {
	fstem="$1"
	outdir="$2"
}

fdir="${fstem%/*}"
fstem="${fstem##*/}"

R="${fdir}/$fstem.?HR"
T="${fdir}/$fstem.?HT"
Z="${fdir}/$fstem.?HZ"

evnm=$(echo $fstem | cut -d'_' -f 1)
evdate="${evnm:0:4}-${evnm:4:2}-${evnm:6:2}"
evtime="${evnm:8:2}:${evnm:10:2}:${evnm:12:2}"

outstem="${outdir}/$fstem"

main() {
	# echo "=============================================================================="
	# echo "-                           QUASI LOVE WAVE ANALYSIS                         -"
	# echo "------------------------------------------------------------------------------"
	# echo "-                                Tom Merry 2023                              -"
	# echo "=============================================================================="
	echo $fstem
	mkdir tmp
	cp $R $T $Z ./tmp/
	fstem=./tmp/$fstem
	R=$fstem.?HR
	T=$fstem.?HT
	Z=$fstem.?HZ
	write_sac_macro
	execute
	write_results
	plot_results
	clean
	
}



write_sac_macro() {
	cat <<- EOF > ql.m
	setbb fstem \$1
	setbb R %fstem%.?HR T %fstem%.?HT Z %fstem%.?HZ
	r %T
	setbb dist &1,dist baz &1,baz
	setbb stlo &1,stlo stla &1,stla
	setbb evlo &1,evlo evla &1,evla
	setbb evdp &1,evdp mag &1,mag
	setbb stnm &1,kstnm net &1,knetwk
	setbb gcarc &1,gcarc
	
	setbb tGmin (div &1,dist 4.7)
	setbb tGmax (div &1,dist 3.8)
	ch user1 %tGmin
	ch kuser1 tGmin
	ch user2 %tGmax
	ch kuser2 tGmax
	mtw %tGmin %tGmax
	markptp l 250 to t0
	setbb Gwind1 ((min &1,t0 &1,t1 ) - 100)
	setbb Gwind2 ((max &1,t0 &1,t1 ) + 100)
	wh
	
	* Find Love arrival
	cut %Gwind1 %Gwind2
	r %T
	cut off
	abs
	mtw b e
	markptp l 250 to t2
	setbb tl &1,t3 
	setbb Gamp &1,user0
	r %T
	div %Gamp
	w append _norm
	
	* Find background noise
	cut 300 600
	r %T
	cut off
	mtw b e
	abs
	markptp l 250 to t2
	setbb noise &1,user0
	
	setbb snr ( %Gamp / %noise )

	* find Rayleigh arrival
	r %Z
	setbb tRmin (max ( div &1,dist 4.1 ) ( %tl + 100 ) )
	setbb tRmax (div &1,dist 3.3)
	ch user1 %tRmin
	ch kuser1 tRmin
	ch user2 %tRmax
	ch kuser2 tRmax
	mtw %tRmin %tRmax
	markptp l 250 to t0
	setbb Rwind1 ((min &1,t0 &1,t1 ) - 100)
	setbb Rwind2 ((max &1,t0 &1,t1 ) + 100)
	wh
	cut %Rwind1 %Rwind2
	r %Z
	cut off
	abs
	mtw b e
	markptp l 250 to t2
	setbb tr &1,t3
	setbb Ramp &1,user0
	
	* Calculate H/V for Rayleigh
	cut ( %tr - 100 ) ( %tr + 100 )
	r %R %Z
	cut off
	setbb HVR ( (&1,depmax  + ( &1,depmin * -1 )) / (&2,depmax  + ( &2,depmin * -1 )) )

	* normalise R by Rayleigh amplitude
	r %R
	* setbb max (max &1,depmax (mul &1,depmin -1))
	* div %max
	div %Ramp
	w append _norm

	* calc gradient of R
	r %R
	dif
	w append _diff
	setbb max (max &1,depmax (mul &1,depmin -1))
	div %max
	mul -1
	w append _diff_norm
	
	* calc hilbert of R
	r %R
	hilbert
	w append _hilb
	* setbb max (max &1,depmax (mul &1,depmin -1))
	* div %max
	div %Ramp
	w append _hilb_norm
	
	* calc envelope of R
	r %R
	setbb max (max &1,depmax (mul &1,depmin -1))
	div %max
	envelope
	w append _env
	

	* Normalise Z
	r %Z
	div %Ramp
	* div %Gamp
	w append _norm
	
	* calc envelope of Z
	r %Z
	setbb max (max &1,depmax (mul &1,depmin -1))
	div %max
	* div %Gamp
	envelope
	w append _env
	
	* sum env
	do renv wild %R%_env
	addf \$renv
	enddo
	w append _R_env_sum
	
	
	* Get Rayleigh amplitude
	r %Z
	sqr
	w append _sqr
	r %R
	sqr
	do zsqr wild %Z%_sqr
	addf \$zsqr
	enddo
	sqrt
	setbb max (max &1,depmax (mul &1,depmin -1))
	div %max
	w append _Z_mod

	* Sum the differentiated R and the normalised Z
	r %R%_diff_norm
	do znorm wild %Z%_norm
	binoperr npts warning
	addf \$znorm
	enddo
	w append _Z_norm_sum
	
	* Sum the Hilberted R and the normalised Z
	r %R%_hilb_norm
	do znorm wild %Z%_norm
	binoperr npts warning
	addf \$znorm
	enddo
	w append _Z_norm_sum


	* window the T about the Love arrival and correlate with the stack
	cut %Gwind1 %Gwind2
	r more %T
	cut off
	div 1 (max &2,depmax (mul &2,depmin -1))
	correlate master 2
	abs

	* Now pick the maximum value, but it has to be between 0 and the Rayleigh pick
	* Let's say it has to be at least 50% way before the Rayleigh pick?
	setbb Rdelt ( %tR - %tL )
	ch t2 %Rdelt kt2 R1
	mtw 0 ( %Rdelt * 0.5 )
	markptp l 200 to t8
	ch t8
	setbb qldelt &1,t9
	setbb tql ( %qldelt + %tl )
	setbb corrmax &1,depmax
	w append _corr
	
	* Is this negative or positive correlation?
	r %R%_hilb_norm_Z_norm_sum
	cut %Gwind1 %Gwind2
	r more %T
	cut off
	div 1 (max &2,depmax (mul &2,depmin -1))
	correlate master 2
	dc 2
	ch t0 0
	cutim t0 ( %qldelt - 1 ) ( %qldelt + 1 )
	if &1,depmax GE 0
	setbb corramp &1,depmax
	else
	setbb corramp &1,depmin
	endif
	
	
	* Calc QL amplitude and H/V
	cut ( %tql - 100 ) ( %tql + 100 )
	r %R %Z
	cut off
	setbb qlamp ( max &1,depmax ( &1,depmin * -1 ) &2,depmax ( &2,depmin * -1 ) )
	setbb qlrelamp ( %qlamp / %Gamp )
	setbb HVQL ( (&1,depmax  + ( &1,depmin * -1 )) / (&2,depmax  + ( &2,depmin * -1 )) )
	
	* and distance
	setbb dx ( %qldelt * %dist  / ( %tR - %tL ) )
	setbb dxx ( %dx * ( sine ( %baz * ( pi ) / 180 ) ) )
	setbb dxy ( %dx * ( cosine ( %baz * ( pi ) / 180 ) ) )
	
	getbb
	getbb to %fstem%.bbinfo

	* xlim -1000 1000
	* p1
	EOF
}

execute() {
	$sac <<-EOF
	m ql.m $fstem
	echo off
	q
	EOF
}

clean() {
	cp ./tmp/tmp.pdf "${outstem}.pdf"
	rm -r ./tmp/
	rm ql.m
}

write_results() {
	echo "writing results"
	# cat "$fstem".bbinfo
	bb="$fstem".bbinfo
	stlo=$(awk '$1=="STLO" {print $3*1}' $bb)
	stla=$(awk '$1=="STLA" {print $3*1}' $bb)
	stnm=$(awk '$1=="STNM" {print $3}' $bb)
	net=$(awk '$1=="NET" {print $3}' $bb)
	evlo=$(awk '$1=="EVLO" {print $3*1}' $bb)
	evla=$(awk '$1=="EVLA" {print $3*1}' $bb)
	evdp=$(awk '$1=="EVDP" {print $3*1}' $bb)
	mag=$(awk '$1=="MAG" {print $3*1}' $bb)
	gcarc=$(awk '$1=="GCARC" {print $3*1}' $bb)
	baz=$(awk '$1=="BAZ" {print $3*1}' $bb)
	
	tR=$(awk '$1=="TR" {print $3*1}' $bb)
	tL=$(awk '$1=="TL" {print $3*1}' $bb)
	tQ=$(awk '$1=="TQL" {print $3*1}' $bb)
	Gamp=$(awk '$1=="GAMP" {print $3*1}' $bb)
	snr=$(awk '$1=="SNR" {print $3*1}' $bb)
	Ramp=$(awk '$1=="RAMP" {print $3*1}' $bb)
	QLamp=$(awk '$1=="QLAMP" {print $3*1}' $bb)
	qlrelamp=$(awk '$1=="QLRELAMP" {print $3*1}' $bb)
	HVR=$(awk '$1=="HVR" {print $3*1}' $bb)
	HVQL=$(awk '$1=="HVQL" {print $3*1}' $bb)
	
	
	corrmax=$(awk '$1=="CORRMAX" {print $3*1}' $bb)
	corramp=$(awk '$1=="CORRAMP" {print $3*1}' $bb)
	rwind1=$(grep RWIND1 $bb | awk '{print $3}')
	rwind2=$(grep RWIND2 $bb | awk '{print $3}')
	gwind1=$(grep GWIND1 $bb | awk '{print $3}')
	gwind2=$(grep GWIND2 $bb | awk '{print $3}')
	trmin=$(grep TRMIN $bb | awk '{print $3}')
	trmax=$(grep TRMAX $bb | awk '{print $3}')
	tgmin=$(grep TGMIN $bb | awk '{print $3}')
	tgmax=$(grep TGMAX $bb | awk '{print $3}')
	qldelt=$(grep QLDELT $bb | awk '{print $3*1}')
	
	dx=$(awk '$1=="DX" {print $3*1}' $bb)
	dxx=$(awk '$1=="DXX" {print $3*1}' $bb)
	dxy=$(awk '$1=="DXY" {print $3*1}' $bb)
	
	sclo=$(echo $dxx $dxy | gmt mapproject -JE${stlo}/${stla}/10 -Rg -I -C -Fk | awk '{print $1}')
	scla=$(echo $dxx $dxy | gmt mapproject -JE${stlo}/${stla}/10 -Rg -I -C -Fk | awk '{print $2}')
	
	res=${outstem}.results
	
	cat <<- EOF > $res
	EVTIME EVLO EVLA EVDP MAG STLO STLA GCARC BAZ SNR G1 R1 DTQL DX SCLO SCLA AMPL H/V(R) H/V(QL) CORRAMP
	${evdate}T${evtime} $evlo $evla $evdp $mag $stlo $stla $gcarc $baz $snr $tL $tR $qldelt $dx $sclo $scla $qlrelamp $HVR $HVQL $corramp
	EOF
	
}



plot_window() {
	min=$1; max=$2; color=$3
	gmt plot -G$color -L <<-EOF
	$min -10
	$max -10
	$max 10
	$min 10
	EOF
}

plot_tl_tr() {
	#statements
	gmt plot -W1p,darkgreen,- <<-EOF
	$tL -10
	$tL 10
	EOF
	gmt plot -W1p,blue,- <<-EOF
	$tR -10
	$tR 10
	EOF
	gmt plot -W1p,red,. <<-EOF
	$tQ -10
	$tQ 10
	EOF
}




plot_results() {
	tmin=$(echo "$tL - 1500" | bc -l)
	tmax=$(echo "$tR + 1500" | bc -l)
	
	gmt begin tmp/tmp
	gmt basemap -JX20/3 -R$tmin/$tmax/-1.1/1.1 -BWesn -Bya1f.5 -Bxa1000f500 -Y30c
	gmt text -N -F+cTC -D0/1.7 <<-EOF
	$evdate $evtime   EVLO: $(printf "%.2f" $evlo)@~\260@~   EVLA: $(printf "%.2f" $evla)@~\260@~   EVDP: $(printf "%.0f" $evdp) km   M: ${mag}
	EOF
	gmt text -N -F+cTC -D0/1.2 <<-EOF
	${net}.${stnm}   STLO: $(printf "%.2f" $stlo)@~\260@~   STLA: $(printf "%.2f" $stla)@~\260@~   \
		GCARC: $(printf "%.2f" $gcarc)@~\260@~   BAZ: $(printf "%.1f" $baz)@~\260@~   SNR: $(printf "%.0f" $snr) 
	EOF
	gmt text -N -F+cTC -D0/0.7 <<-EOF
	QL/G = $( printf "%.2f" $qlrelamp )  H/V@-R@- = $( printf "%.2f" $HVR )  H/V@-QL@- = $( printf "%.2f" $HVQL )
	EOF

	# sac2xy ${T}_norm /dev/stdout | head
	# sac2xy ${Z}_norm /dev/stdout | head
	
	
	# sac2xy ${T}_norm /dev/stdout | gmt plot -h1 -W1p,darkgreen
	gmt sac ${T}_norm -W1p,darkgreen
	plot_window $tgmin $tgmax green@80
	echo "T" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold,darkgreen
	plot_tl_tr
	
	# sac2xy ${Z}_norm /dev/stdout | gmt plot -h1 -W1p,blue -Y-3.2c -BWesn -Bya1f.5 -Bxa1000f500
	gmt sac ${Z}_norm -W1p,blue -Y-3.2c -BWesn -Bya1f.5 -Bxa1000f500
	plot_window $trmin $trmax blue@80
	# sac2xy ${R}_Z_mod /dev/stdout | gmt plot -h1 -Wthinnest,black
	# gmt sac ${R}_Z_mod -Wthinnest,black
	echo "Z" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold,blue
	plot_tl_tr
	
	gmt sac ${R}_norm -W1p,red -Y-3.2c -BWesn -Bya1f.5 -Bxa1000f500
	echo "R" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold,red
	plot_tl_tr
	
	gmt sac ${R}_hilb_norm -W1p,red,- -Gp+ggrey@80 -Y-3.2c -BWesn -Bya1f.5 -Bxa1000f500 
	gmt sac ${Z}_norm -W1p,blue -Gp+ggrey@80
	# sac2xy ${R}_env /dev/stdout | gmt plot -h1 -W1p,black
	# sac2xy ${Z}_env /dev/stdout | gmt plot -h1 -W1p,orange
	
	echo "Z" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold,blue
	echo "H(R)" | gmt text -D0.2/0.2 -F+cBL+f12p,Helvetica-Bold,red
	plot_tl_tr
	
	gmt basemap -R$tmin/$tmax/-2.1/2.1 -BWeSn -Bya2f1 -Bxa1000f500+l"Time (s) since earthquake" -Y-3.2c
	echo "${R}_hilb_norm_Z_norm_sum"
	gmt sac ${R}_hilb_norm_Z_norm_sum -W1p,magenta
	# gmt sac ${Z}_R_env_sum -W1p,orange
	if (( $(echo "$corramp > 0" | bc -l) )); then
		echo $corramp bigger than zero
		Tscale=1
	else
		echo $corramp smaller than zero
		Tscale=-1
	fi
	gmt sac ${T}_norm -T+s${qldelt} -M${Tscale}/0 -W1p,darkgreen 
	echo "Z @~\053@~ H(R) stack" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold,magenta
	echo "T time shifted by @~\144@~t (${qldelt} s)" | gmt text -D0.2/0.2 -F+cBL+f12p,Helvetica-Bold,darkgreen
	plot_tl_tr
	# sac2xy ${R}_norm /dev/stdout | gmt plot tmp.tmp -h1 -W1p,blue -Y-3.2c -BWesn -Bya1f.5 -Bxa1000f500
	
	gmt basemap -R$(echo "$tmin - $tL" | bc -l)/$(echo "$tmax - $tL" | bc -l)/0/$corrmax -BWeSn -Byaf -Bxa1000f500+l"Time lag (s)" -Y-5c
	gmt sac ${R}_hilb_norm_corr -W1p,black
	gmt basemap -R$tmin/$tmax/-2.1/2.1 -A > /dev/null
	plot_tl_tr
	
	r=$(echo "$Gamp * 1.1" | bc -l)
	gmt basemap -Y16c -X22c -JX4 -R-$r/$r/-$r/$r -BNWes -Byg100000+l"Radial" -Bxg100000+l"Transverse"
	paste <(sac2xy $T /dev/stdout ) <(sac2xy $R /dev/stdout) | awk -v t=$tL '$1 > (t -75) && $1 < (t +75)' | gmt plot -i1,3 -W1p -S~n4:+sc1c+d
	echo "G1" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold
	
	r=$(echo "$Ramp * 1.1" | bc -l)
	gmt basemap -Y-5.5 -R-$r/$r/-$r/$r -BNWes -Byg100000+l"Vertical" -Bxg100000+l"Radial"
	paste <(sac2xy $R /dev/stdout ) <(sac2xy $Z /dev/stdout) | awk -v t=$tR '$1 > (t -75) && $1 < (t +75)' | gmt plot -i1,3 -W1p -S~n2:+st0.5+an
	echo "R1" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold
	
	r=$(echo "$QLamp * 1.1" | bc -l)
	gmt basemap -Y-5.5 -R-$r/$r/-$r/$r -BNWes -Byg100000+l"Vertical" -Bxg100000+l"Radial"
	paste <(sac2xy $R /dev/stdout ) <(sac2xy $Z /dev/stdout) | awk -v t=$tQ '$1 > (t -75) && $1 < (t +75)' | gmt plot -i1,3 -W1p
	echo "QL" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold
	
	
	gmt basemap -Y-5.5 -R-$r/$r/-$r/$r -BNWes -Byg100000+l"Vertical" -Bxg100000+l"H(Radial)"
	paste <(sac2xy ${R}_hilb /dev/stdout ) <(sac2xy $Z /dev/stdout) | awk -v t=$tQ '$1 > (t -75) && $1 < (t +75)' | gmt plot -i1,3 -W1p
	echo "QL" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold
	
	
	
	
	
	gmt coast -JE${stlo}/${stla}/10 -Rg -Wthinnest -Ggrey@50 -A100000 -B -X5 -Y5
	echo $stlo $stla | gmt plot -St0.2c -Gred -Wthinnest
	echo $evlo $evla | gmt plot -Sa0.2c -Gred -Wthinnest
	gmt plot -W0.5p,red <<-EOF
	$stlo $stla
	$evlo $evla
	EOF
	echo $sclo $scla | gmt plot -Sc0.3c -W1p,purple
	
	gmt end
	
	echo "QL amplitude ${qlrelamp}"
}



main
