#!/usr/bin/bash
# * Look for potential Quasi-Love waves
#
# * USAGE: 
# * In a shell script try:
# * source quasilove_fns.sh
# * qlmain filestem outdir [plotoff]
# * where the components are expected to
# * be in the file format [filestem].?H[R,T,Z]
# * The optional 'plotoff' will stop it from plotting the results, to save time if not needed
# *
# * waveforms should be filtered with a fmax of 0.01 Hz
# *
# * Love and Rayleigh waves are picked based on windows covering
# * approx max and min group velocities at periods >100s based on
# * model GDM52 (Ekstrom 2011); i.e.,
# * 3.8 km/s < V(G) < 4.7 km/s
# * 3.4 km/s < V(R) < 4.1 km/s
# *

# If you need to point to a particular SAC installation, do it here.
# sac=/usr/local/LINUX_sac/bin/sacinit.sh
sac=sac

def_variables() {
	fstem="$1"
	outdir="$2"
    fdir="${fstem%/*}"
    fstem="${fstem##*/}"

    R="${fdir}/$fstem.?HR"
    T="${fdir}/$fstem.?HT"
    Z="${fdir}/$fstem.?HZ"

    evnm=$(echo $fstem | cut -d'_' -f 1)
    evdate="${evnm:0:4}-${evnm:4:2}-${evnm:6:2}"
    evtime="${evnm:8:2}:${evnm:10:2}:${evnm:12:2}"

    outstem="${outdir}/$fstem"
}


qlmain() {
    if [ -z $2 ]; then
        echo "Expected quasilove.sh [filestem] [outdir]"
        echo "Try again"
        exit
    fi
	plotting=1
	if [[ ! -z $3 ]]; then
		if [[ $3 == "plotoff" ]]; then
			plotting=0
		else
			echo "Unknown option $3, did you mean to write 'plotoff' to disable plotting?"
			exit
		fi
	fi
    def_variables "$1" "$2"
	echo $fstem
	mkdir -p tmp
	cp $R $T $Z ./tmp/
	fstem=./tmp/$fstem
	R=$fstem.?HR
	T=$fstem.?HT
	Z=$fstem.?HZ
	write_sac_macro
    polarisation_code
	echo hello sac is $sac
	execute
	check_it_ran
	write_results
	if (( plotting )); then
		plot_results
	fi
	clean
}

polarisation_code() {
	# Very simple python code that calculates & spits out the polarisation angle from the covariance matrix
	# it takes the TT, RR and TR covariance values as input
    cat << EOF > ./tmp/pol.py
import sys, math
tt, rr, tr = [float(x) for x in sys.argv[1:4]]
l1 = (tt + rr + ((tt + rr)**2 - 4 * (tt * rr - tr**2))**0.5 ) / 2
if tr != 0:
	v1 = tr
	v2 = l1 - tt
else:
	v1 = 1
	v2 = 0
angle = math.atan(v2/v1) * 180 / math.pi
print(angle)
EOF
}

write_sac_macro() {
	cat <<-EOF > ql.m
	setbb fstem \$1
	setbb R %fstem%.?HR T %fstem%.?HT Z %fstem%.?HZ

	* Read some variables
	r %T
	setbb dist &1,dist baz &1,baz
	setbb stlo &1,stlo stla &1,stla
	setbb evlo &1,evlo evla &1,evla
	setbb evdp &1,evdp mag &1,mag
	setbb stnm &1,kstnm net &1,knetwk
	setbb gcarc &1,gcarc
	
	* Define search window for Love
	setbb tGmin (%dist% / 4.7)
	setbb tGmax (%dist% / 3.8)
	mtw %tGmin %tGmax
	markptp l 250 to t0
	setbb Gwind1 ((min &1,t0 &1,t1 ) - 100)
	setbb Gwind2 ((max &1,t0 &1,t1 ) + 100)

    * Find Love arrival
    abs
	mtw %Gwind1 %Gwind2
	markptp l 250 to t2
	setbb tl &1,t3
	setbb Gamp &1,user0

    * Calculate SNR
    mtw 300 600
	markptp l 250 to t2
	setbb noise &1,user0
    setbb snr ( %Gamp / %noise )

    * calc G1 polarisation
    * First save G1 window
    echo on
    cut %Gwind1 %Gwind2
    r %T %R
    cut off
    w append _G1
    echo off

    * Calc covariance matrix
    * TR
    r %T%_G1
    do RG wild %R%_G1
	binoperr npts warning
	mulf \$RG
	enddo
    int
    cutim e -0.5 0
    setbb covTR &1,depmax
    * TT
    r %T%_G1
    do TG wild %T%_G1
	binoperr npts warning
	mulf \$TG
	enddo
    int
    cutim e -0.5 0
    setbb covTT &1,depmax
    * RR
    r %R%_G1
    do RG wild %R%_G1
	binoperr npts warning
	mulf \$RG
	enddo
    int
    cutim e -0.5 0
    setbb covRR &1,depmax

    sc python ./tmp/pol.py %covTT %covRR %covTR > ./tmp/G1.pol
	
    readtable ./tmp/G1.pol
    setbb G1pol &1,depmax

    * Normalise T
	r %T
	div %Gamp
	w append _norm
	
    * Define search window for Rayleigh
	r %Z
	setbb tRmin (max ( %dist% / 4.1 ) ( %tl + 100 ) )
	setbb tRmax (%dist% / 3.3)
	mtw %tRmin %tRmax
	markptp l 250 to t0
	setbb Rwind1 ((min &1,t0 &1,t1 ) - 100)
	setbb Rwind2 ((max &1,t0 &1,t1 ) + 100)

	* find Rayleigh arrival
	abs
	mtw %tRmin %tRmax
	markptp l 250 to t2
	setbb tr &1,t3
	setbb Ramp &1,user0

    * calc R1 polarisation
    * First save R1 window
    cut %Rwind1 %Rwind2
    r %T %R
    cut off
    w append _R1

    * Calc covariance matrix
    * TR
    r %T%_R1
    do RR wild %R%_R1
	binoperr npts warning
	mulf \$RR
	enddo
    int
    cutim e -0.5 0
    setbb covTR &1,depmax
    * TT
    r %T%_R1
    do TR wild %T%_R1
	binoperr npts warning
	mulf \$TR
	enddo
    int
    cutim e -0.5 0
    setbb covTT &1,depmax
    * RR
    r %R%_R1
    do RR wild %R%_R1
	binoperr npts warning
	mulf \$RR
	enddo
    int
    cutim e -0.5 0
    setbb covRR &1,depmax

    sc python ./tmp/pol.py %covTT %covRR %covTR > ./tmp/R1.pol
    readtable ./tmp/R1.pol
    setbb R1pol ( 90 + &1,depmax )

    * Normalise Z
	r %Z
	div %Ramp
	* div %Gamp
	w append _norm
	
	* Calculate H/V for Rayleigh
	cut ( %tr - 100 ) ( %tr + 100 )
	r %R %Z
	cut off
	setbb HVR ( (&1,depmax  + ( &1,depmin * -1 )) / (&2,depmax  + ( &2,depmin * -1 )) )
	

	* normalise R by Rayleigh amplitude
	r %R
	div %Ramp
	w append _norm
	
	* calc hilbert of R
	r %R
	hilbert
	w append _hilb
	div %Ramp
	w append _hilb_norm
	
	
	* Sum the Hilberted R and the normalised Z
	r %R%_hilb_norm
	do znorm wild %Z%_norm
	binoperr npts warning
	addf \$znorm
	enddo
	w append _Z_norm_sum


	* window the T about the Love arrival and correlate with the stack
	cut ( %Gwind1 - 50 ) ( %Gwind2 + 50 )
	r more %T
	cut off
	div 1 (max &2,depmax ( &2,depmin * -1 ))
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
    cut ( %Gwind1 - 50 ) ( %Gwind2 + 50 )
	r more %T
	cut off
	div 1 (max &2,depmax ( &2,depmin * -1 ))
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
	
	getbb to %fstem%.bbinfo
	EOF
}

execute() {
	echo $sac
	$sac <<-EOF > /dev/null
	m ql.m $fstem
	echo off
	q
	EOF
	
}

check_it_ran(){
	if [ ! -f "${fstem}.bbinfo" ]; then
		echo "Something went wrong, no sac bbinfo file created"
		exit 1
	fi
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
	stnm=$(awk '$1=="STNM" {print $3}' $bb | tr -d \')
	net=$(awk '$1=="NET" {print $3}' $bb | tr -d \')
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
    G1pol=$(awk '$1=="G1POL" {print $3*1}' $bb)
    R1pol=$(awk '$1=="R1POL" {print $3*1}' $bb)
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

	writedir=${res%/*}
	if [ ! -d "$writedir" ]; then
		mkdir -p "$writedir"
	fi
	
	cat <<- EOF > $res
	EVTIME EVLO EVLA EVDP MAG STLO STLA GCARC BAZ SNR G1 R1 DTQL DX SCLO SCLA AMPL H/V(R) H/V(QL) CORRAMP G1POL R1POL
	${evdate}T${evtime} $evlo $evla $evdp $mag $stlo $stla $gcarc $baz $snr $tL $tR $qldelt $dx $sclo $scla $qlrelamp $HVR $HVQL $corramp $G1pol $R1pol
	EOF
}


sac2xy(){
	sacf=$1
	npts=$(gmt convert $sacf -bi1i+l -q79:79)
	delta=$(gmt convert $sacf -bi1f+l -q0:0)
	tb=$(gmt convert $sacf -bi1f+l -q5:5)
	gmt convert $1 -bi1f+l -hi632 | awk -v npts=$npts -v delta=$delta -v tb=$tb 'NR<=npts {
		t = tb + (NR - 1) * delta
		print t, $1
	}'
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
	R1Pol = $( printf "%.2f" $R1pol ) G1Pol = $( printf "%.2f" $G1pol )  QL/G = $( printf "%.2f" $qlrelamp )  H/V@-R@- = $( printf "%.2f" $HVR )  H/V@-QL@- = $( printf "%.2f" $HVQL )
	EOF

	# sac2xy ${T}_norm /dev/stdout | head
	# sac2xy ${Z}_norm /dev/stdout | head
	
	plot_window $tgmin $tgmax green@80
	# sac2xy ${T}_norm /dev/stdout | gmt plot -h1 -W1p,darkgreen
	gmt sac ${T}_norm -W1p,darkgreen
	
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
	paste <(sac2xy $T /dev/stdout ) <(sac2xy $R /dev/stdout) | 
		awk -v t=$tL '$1 > (t -75) && $1 < (t +75)' | 
		gmt plot -i1,3 -W1p -S~n4:+sc1c+d
	echo "G1" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold
	
	r=$(echo "$Ramp * 1.1" | bc -l)
	gmt basemap -Y-5.5 -R-$r/$r/-$r/$r -BNWes -Byg100000+l"Vertical" -Bxg100000+l"Radial"
	paste <(sac2xy $R /dev/stdout ) <(sac2xy $Z /dev/stdout) | 
		awk -v t=$tR '$1 > (t -75) && $1 < (t +75)' | 
		gmt plot -i1,3 -W1p -S~n2:+st0.5+an
	echo "R1" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold
	
	r=$(echo "$QLamp * 1.1" | bc -l)
	gmt basemap -Y-5.5 -R-$r/$r/-$r/$r -BNWes -Byg100000+l"Vertical" -Bxg100000+l"Radial"
	paste <(sac2xy $R /dev/stdout ) <(sac2xy $Z /dev/stdout) | 
		awk -v t=$tQ '$1 > (t -75) && $1 < (t +75)' | 
		gmt plot -i1,3 -W1p
	echo "QL" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold
	
	
	gmt basemap -Y-5.5 -R-$r/$r/-$r/$r -BNWes -Byg100000+l"Vertical" -Bxg100000+l"H(Radial)"
	paste <(sac2xy ${R}_hilb /dev/stdout ) <(sac2xy $Z /dev/stdout) | 
		awk -v t=$tQ '$1 > (t -75) && $1 < (t +75)' | 
		gmt plot -i1,3 -W1p
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

