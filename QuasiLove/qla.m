* Look for potential Quasi-Love waves
* USAGE: qla filestem, where the components are expected to
* be in the file format [filestem].?H[R,T,Z]
*
* waveforms should be filtered with a fmax of 0.01 Hz
*
* Love and Rayleigh waves are picked based on windows covering
* approx max and min group velocities at periods >100s based on
* model GDM52 (Ekstrom 2011); i.e.,
* 3.8 km/s < V(G) < 4.7 km/s
* 3.4 km/s < V(R) < 4.1 km/s
*

setbb fstem $1
setbb R %fstem%.?HR T %fstem%.?HT Z %fstem%.?HZ


message "=============================================================================="
message "-                           QUASI LOVE WAVE ANALYSIS                         -"
message "------------------------------------------------------------------------------"
message "-                                Tom Merry 2023                              -"
message "=============================================================================="

* find Love arrival
r %T
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
setbb max (max &1,depmax (mul &1,depmin -1))
div %max
w append _norm

* find Rayleigh arrival
r %Z
setbb tRmin (div &1,dist 4.1)
setbb tRmax (div &1,dist 3.4)
ch user1 %tRmin
ch kuser1 tRmin
ch user2 %tRmax
ch kuser2 tRmax
mtw %tRmin %tRmax
markptp l 250 to t0
setbb Rwind1 ((min &1,t0 &1,t1 ) - 100)
setbb Rwind2 ((max &1,t0 &1,t1 ) + 100)
wh

* calc gradient of R
r %R
dif
setbb max (max &1,depmax (mul &1,depmin -1))
div %max
w append _norm
mul -1
w change _norm _diff

* Normalise Z
r %Z
setbb max (max &1,depmax (mul &1,depmin -1))
div %max
w append _norm

* Sum the differentiated R and the normalised Z
r %R%_diff
do znorm wild %Z%_norm
binoperr npts warning
addf $znorm
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
* Let's say it has to be at least 150 s before the Rayleigh pick??
setbb Rdelt ( %Rwind1 - %Gwind1 )
ch t2 %Rdelt kt2 R1
mtw 0 ( %Rdelt - 150 )
markptp l 200 to t8
ch t8
setbb qldelt &1,t9

xlim -1000 1000
p1



! sc echo "0.5 0.5" | gmt plot -JX5 -R-1/1/-1/1 -Sa0.5c -pdf tmp2


sc gmt begin tmp
sc gmt basemap -JX20/3 -R500/6000/-1.1/1.1 -BWesn -Bya1f.5 -Bxa1000f500 -Y30c

sc sac2xy %T%_norm /dev/stdout > tmp.tmp
sc gmt plot tmp.tmp -h1 -W1p,darkgreen
sc sac2xy %Z%_norm /dev/stdout > tmp.tmp
sc gmt plot tmp.tmp -h1 -W1p,blue -Y-3.2c -BWesn -Bya1f.5 -Bxa1000f500
sc sac2xy %R%_norm /dev/stdout > tmp.tmp
sc gmt plot tmp.tmp -h1 -W1p,blue -Y-3.2c -BWesn -Bya1f.5 -Bxa1000f500




sc gmt end show