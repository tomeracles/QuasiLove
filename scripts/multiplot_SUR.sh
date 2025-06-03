#!/bin/bash

mkdir tmp 2> /dev/null
sac=/usr/local/LINUX_sac/bin/sacinit.sh

stn=$1


results=./Africa/results/accepted_results.txt

normalise() {
    Gamp=$(sac2xy Africa/waveforms/$f.?HT /dev/stdout | 
            awk -v t=$tG1 '$1 > (t - 0.5) && $1 < (t + 0.5) {print $2}' | 
            awk 'NR==1')
    $sac <<- EOF
    r Africa/waveforms/$f.?H[R,T,Z]
    div $Gamp
    w append _norm
    q
EOF
}

stack_stack() {
    cp Africa/waveforms/$f.?HR Africa/waveforms/$f.?HZ tmp/
    zf=$(echo tmp/$f.?HZ)
    $sac <<- EOF
*    echo on
    r tmp/$f.?HR
    hilbert
*    do zcmp wild tmp/$f.?HZ
*	binoperr npts warning
	addf $zf
*    enddo
    div $Gamp
    w Africa/waveforms/$f.qlstack
	
    echo off
    q
EOF
    rm tmp/$f.?HR tmp/$f.?HZ
}

get_stuff() {
    # grep _G.TAM. ./Africa/results/logfile.txt | 
    #     awk '$2==1 {print"Africa/waveforms/"$1".?HT"}' > tmp/sac.list

    alllist=$(grep _$stn$ ./Africa/logfile.txt | awk '{print $1}')
    glist=$(grep _$stn. ./Africa/results/logfile.txt | awk '$2==1{print $1}')

    grep _$stn. ./Africa/results/logfile.txt | awk '{print $1}' | head

    echo "" > tmp/allT.saclst
    echo "" > tmp/allR.saclst
    echo "" > tmp/allZ.saclst
    echo "" > tmp/goodT.saclst
    echo "" > tmp/goodR.saclst
    echo "" > tmp/goodZ.saclst
    echo "" > tmp/allqlstack.saclst
    echo "" > tmp/goodqlstack.saclst

    for f in $alllist; do
        # echo $f
        snr=$(awk 'NR==2{print $10}' Africa/results/$f.results)
        qlamp=$(awk 'NR==2{print $17}' Africa/results/$f.results)
        if (( $(echo "$snr > 10" | bc -l) )) && (( $(echo "$qlamp < 0.5" | bc -l) )); then
            baz=$(awk 'NR==2{print $9}' Africa/results/$f.results)
            tG1=$(awk 'NR==2{print $11}' Africa/results/$f.results)
            

            
            
            normalise

            stack_stack

            echo Africa/waveforms/$f.?HT_norm -$tG1 $baz >> tmp/allT.saclst
            echo Africa/waveforms/$f.?HR_norm -$tG1 $baz >> tmp/allR.saclst
            echo Africa/waveforms/$f.?HZ_norm -$tG1 $baz >> tmp/allZ.saclst
            echo Africa/waveforms/$f.qlstack -$tG1 $baz >> tmp/allqlstack.saclst
            # sac2xy Africa/waveforms/$f.?HT /dev/stdout | 
            #     awk -v t=$tG1 a=$Gamp 'NR>1{print ($1 - t), $2/a}' >  Africa/waveforms/$f.T.xy
            # sac2xy Africa/waveforms/$f.?HR /dev/stdout | 
            #     awk -v t=$tG1 a=$Gamp 'NR>1{print ($1 - t), $2/a}' >  Africa/waveforms/$f.R.xy
            # sac2xy Africa/waveforms/$f.?HZ /dev/stdout | 
            #     awk -v t=$tG1 a=$Gamp 'NR>1{print ($1 - t), $2/a}' >  Africa/waveforms/$f.Z.xy
        fi
    done
    for f in $glist; do
        # echo $f
        baz=$(awk 'NR==2{print $9}' Africa/results/$f.results)
        tG1=$(awk 'NR==2{print $11}' Africa/results/$f.results)
        snr=$(awk 'NR==2{print $10}' Africa/results/$f.results)
        echo Africa/waveforms/$f.?HT_norm -$tG1 $baz >> tmp/goodT.saclst
        echo Africa/waveforms/$f.?HR_norm -$tG1 $baz >> tmp/goodR.saclst
        echo Africa/waveforms/$f.?HZ_norm -$tG1 $baz >> tmp/goodZ.saclst
        echo Africa/waveforms/$f.qlstack -$tG1 $baz >> tmp/goodqlstack.saclst
    done
}

plot_sac() {
    gmt begin tmp/tmp
        gmt basemap -JX10/25 -R-200/550/20/110 -BWrnS \
            -Byaf+l"Backazimuth" -Bxaf+l"Time after G1 (s)"
        gmt sac tmp/allT.saclst -Wthinnest,grey -M1/0 -Gp+gred@80 -Gn+gblue@80
        gmt sac tmp/goodT.saclst -W0.5p,blue -M1/0
        gmt plot -W0.5p,black <<-EOF
        0 -1000
        0 1000
EOF
        echo "Transverse" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold

        gmt basemap -X10 -BlrnS -Bx
        gmt sac tmp/allR.saclst -Wthinnest,grey -M1/0 -Gp+gred@80 -Gn+gblue@80
        gmt sac tmp/goodR.saclst -W0.5p,blue -M1/0
        gmt plot -W0.5p,black <<-EOF
        0 -1000
        0 1000
EOF
        echo "Radial" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold

        gmt basemap -X10 -BlrnS -Bx
        gmt sac tmp/allZ.saclst -Wthinnest,grey -M2/0 -Gp+gred@80 -Gn+gblue@80
        gmt sac tmp/goodZ.saclst -W0.5p,blue -M2/0
        gmt plot -W0.5p,black <<-EOF
        0 -1000
        0 1000
EOF
        echo "Vertical (x2)" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold

        gmt basemap -X10 -BlEnS -Bx+l"Time after G1(s)" -By+l"Backazimuth"
        gmt sac tmp/allqlstack.saclst -Wthinnest,grey -M1/0 -Gp+gred@80 -Gn+gblue@80
        gmt sac tmp/goodqlstack.saclst -W0.5p,blue -M1/0
        gmt plot -W0.5p,black <<-EOF
        0 -1000
        0 1000
EOF
        echo "H(R) + Z" | gmt text -D0.2/-0.2 -F+cTL+f12p,Helvetica-Bold
		
		awk -v stlo=$stlo -v stla=$stla '$6==stlo && $7==stla {print $13, $9}' $results | head
		awk -v stlo=$stlo -v stla=$stla '$6==stlo && $7==stla {print $13, $9}' $results | gmt plot -Sy0.3c -W2p,red
		
    gmt end show
}

plot_map() {
	gmt plot -JE
}




main(){
    # get_stuff
    plot_sac
}

a=$(grep "_II.SUR\t" ./Africa/results/logfile.txt | awk '$2==1 {print $1; exit;}')
echo $a
stlo=$(tail -n 1 ./Africa/results/$a.results | awk '{print $6}')
stla=$(tail -n 1 ./Africa/results/$a.results | awk '{print $7}')
echo $stlo $stla

main