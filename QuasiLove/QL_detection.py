#!/usr/bin/env python

# import obspy

networks = { 'IRIS': [ 'II', 'IU', 'IM', 'TT', 'WM', 'MN', 'AF', 'DZ', 'G', 'GT', 'IP', 'NJ', 'GH',  'NR', 'BX'],
            'GFZ': ['GE'],
            'ORFEUS' : ['DZ'],
            'RESIF' : ['QM']
}

class SAC:
    def __init__(self, launch_path):
        self.launch = launch_path
    def run(self):
        import subprocess
        output = subprocess.run(self.launch, shell=True, input=self.input, capture_output=True, universal_newlines=True)
        return output

def download(station,event,client,inv, outdir):
    from obspy.geodetics import locations2degrees
    from obspy.io.sac.sactrace import SACTrace
    stlo = station.longitude
    stla = station.latitude
    evo = event.origins[0].time
    evnm = f'{evo.year}{evo.month:02d}{evo.day:02d}{evo.hour:02d}{evo.minute:02d}{evo.second:02d}'
    evlo = event.origins[0].longitude
    evla = event.origins[0].latitude
    evdp = event.origins[0].depth
    mag = event.magnitudes[0].mag
    
    dist = locations2degrees(evla, evlo, stla, stlo)
    print(f'downloading for {event.origins[0].time}')
    try:
        st = client.get_waveforms(network=station.network, station = station.code, 
                                  location = '*', channel = '?H?',
                         starttime = event.origins[0].time,
                         endtime = event.origins[0].time + 6000
                        )
    except Exception as e:
        print(e)
        print('no data available')
        return
    
    # Get only the longest traces from each channel
    traces = []
    ids = set([tr.id for tr in st])
    for trid in ids:
        if len(st.select(id=trid)) == 1:
            traces.append(st.select(id=trid)[0])
        else:
            tlen = 0
            for tr in st.select(id=trid):
                if tr.stats.npts > tlen:
                    longest = tr
                    tlen = tr.stats.npts
            traces.append(longest)
    st = obspy.Stream(traces=traces)
    
    # Choose a channel type (based on priority list) and a location (based on longest trace)
    moveon=0
    chs = set([tr.stats.channel[0] for tr in st])
    for ch in ['L', 'M', 'B', 'H', 'V']:
        if ch in chs:
            sttmp = st.select(channel = f'{ch}H?')
            tlen=0
            for loc in set([tr.stats.location for tr in sttmp]):
                if len(st.select(channel = f'{ch}H?', location = loc)) == 3:
                    tr = st.select(channel = f'{ch}H?', location = loc)[0]
                    if tr.stats.npts > tlen:
                        longest = st.select(channel = f'{ch}H?', location = loc)
                        tlen = tr.stats.npts
                    st1 = longest
                    moveon=1          
        if moveon==1:
            break
    if moveon == 0:
        print('no suitable data')
        return
    
    from obspy.geodetics import gps2dist_azimuth
    st1.detrend(type='linear')
    st1.taper(0.1)
    st1.resample(1)
    st1.filter('bandpass',freqmax=0.01,freqmin=0.002,corners=3,zerophase=True)

    evd = f'{evo.year}{evo.month:02d}{evo.day:02d}'
    mag = event.magnitudes[0].mag
    evdp = event.origins[0].depth / 1000
    _, _, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
    
    try:
        st1.rotate('->ZNE', inventory = inv)
    except Exception as e:
        print(e)
        return
    try:
        st1.rotate('NE->RT', back_azimuth = baz)
    except Exception as e:
        print(e)
        return
    
    for tr in st1:
        sac = SACTrace.from_obspy_trace(tr)
        sac.lcalda = 1
        sac.o = otime
        sac.iztype = "io"
        sac.evlo = evlo
        sac.evla = evla
        sac.evdp = evdp
        sac.stlo = stlo
        sac.stla = stla
        sac.mag = mag
        dest = f'{outdir}/{evnm}_{sac.knetwk}.{sac.kstnm}.{sac.kcmpnm}'
        sac.write(dest)


def ql_analysis_from_sac(filestem, outdir, logf, sac):
    # from obspy import read as obspyread
    import time
    from obspy.io.sac.sactrace import SACTrace
    import numpy as np
    from glob import glob

    t0 = time.time()
    tsac = SACTrace.read(f'{filestem}T', byteorder = 'little')
    rsac = SACTrace.read(f'{filestem}R', byteorder = 'little')
    zsac = SACTrace.read(f'{filestem}Z', byteorder = 'little')
    # print(tsac)
    
    #some metadata
    stlo = tsac.stlo
    stla = tsac.stla
    stnm = tsac.kstnm
    net = tsac.knetwk
    delta = tsac.delta
    evo = tsac.o
    evnm = filestem.split("/")[-1].split("_")[0]
    
    evlo = tsac.evlo
    evla = tsac.evla
    evdp = tsac.evdp
    mag = tsac.mag
    gcarc = tsac.gcarc
    baz = tsac.baz
    print(evnm, stnm, net, evlo, evla, evdp, mag)
    
    
    to = tsac.data
    ro = rsac.data
    zo = zsac.data

    t1 = time.time()

    # normalise
    tn = to/max(abs(to))
    rn = ro/max(abs(ro))
    zn = zo/max(abs(zo))
    
    
    print(delta)
    
    print(t1 - t0, 'seconds elapsed')
    
    # Find max Love amplitude
    
    #get search window
    minGv = 3.8 # km/s, approx min from GDM52 (Ekstrom 2011)
    maxGv = 4.7
    minGt = tsac.dist / maxGv
    maxGt = tsac.dist / minGv
    minGi = int((minGt - tsac.b)/delta)
    maxGi = int((maxGt - tsac.b)/delta)
    
    Gwind = tn[minGi:maxGi]
    G1max=max(abs(Gwind))
    G1pos = np.argmax(abs(Gwind)) + minGi
    G1maxt = G1pos * delta
    
    
    # Repeat for Rayleigh
    minRv = 3.4 # km/s, approx min from GDM52 (Ekstrom 2011)
    maxRv = 4.1
    minRt = tsac.dist / maxRv
    maxRt = tsac.dist / minRv
    minRi = int((minRt - rsac.b)/delta)
    maxRi = int((maxRt - rsac.b)/delta)
    
    Rwind = zn[minRi:maxRi]
    R1max = max(abs(Rwind))
    R1pos = np.argmax(abs(Rwind)) + minRi
    R1maxt = R1pos * delta


    #  calculate SNR 
    Tnoise=max(abs(tn[int(300/delta):int(600/delta)]))
    Tsig=G1max
    snr=Tsig/Tnoise
    
    print('SNR = ', snr)

    # get times
    tt = delta * np.arange(len(tn)) + tsac.b
    
    # calc gradient of R
    rgrad = np.gradient(rn) / max(abs(np.gradient(rn)))
    # calc hilbert of R
    from scipy.signal import hilbert
    rhilb = hilbert(rn).imag
    
    # plotting
    
    import pygmt
    fig = pygmt.Figure()
    fig.shift_origin(yshift= 45)

    xsize = 20
    ysize = 3
    yshift = ysize + 0.3

    fig.plot(
        projection = f'X{xsize}/{ysize}',
        region = [500,6000,-1.05,1.05],
        x = tt,
        y = tn,
        pen = '1p,darkgreen',
        frame = ['nWes', 'xa500f100', 'ya1']
    )
    
    # plot Love picking window
    fig.plot(
        x = [minGt, minGt, maxGt, maxGt],
        y = [-10, 10, 10, -10],
        close = True,
        color = 'lightgreen@70',
        pen = '0.5p,lightgreen'
    )
    
    def plot_tl_tr():
        fig.plot(x = [G1maxt, G1maxt],
            y = [-10, 10],
            pen = '1p,darkgreen,-')
        fig.plot(x = [R1maxt, R1maxt],
            y = [-10, 10],
            pen = '1p,orange,-')
    
    plot_tl_tr()
    fig.text(position = 'TL',
             text = 'T',
             font= '12p,OpenSans-Bold,darkgreen',
            offset = '0.2c/-0.2c'
            )
    fig.text(x = G1maxt,
             y = 1.3,
             text = 't@-L',
             font = '10p,Helvetica-Bold,darkgreen',
             no_clip = True
            )
    fig.text(x = R1maxt,
             y = 1.3,
             text = 't@-R',
             font = '10p,Helvetica-Bold,orange',
             no_clip = True
            )
    
    evtxt = f'{evnm[0:4]}-{evnm[4:6]}-{evnm[6:8]} {evnm[8:10]}:{evnm[10:12]}:{evnm[12:14]}'
    hdrtxt1 = f'{evtxt}   EVLO: {evlo:.3f}   EVLA: {evla:.3f}   EVDP: {evdp:.1f} km   M: {mag:.1f}'
    hdrtxt2 = f'{net}.{stnm}   STLO: {stlo:.3f}   STLA: {stla:.3f}   ' + \
            f'GCARC: {gcarc:.1f}@+o@+   BAZ: {baz:.1f}@+o@+   SNR: {snr:.0f} '

    fig.text(position = 'TC',
             text = hdrtxt1,
             offset = '1.2c',
             no_clip = True
            )
    fig.text(position = 'TC',
             text = hdrtxt2,
             offset = '0.8c',
             no_clip = True
            )
            

    fig.shift_origin(yshift= -yshift)

    fig.plot(
        x = tt,
        y = zn,
        pen = '1p,blue',
        frame = ['nWes', 'xa500f100', 'ya1']
    )
    # plot Rayleigh picking window
    fig.plot(
        x = [minRt, minRt, maxRt, maxRt],
        y = [-10, 10, 10, -10],
        close = True,
        color = 'lightblue@70',
        pen = '0.5p,lightblue'
    )
    plot_tl_tr()
    fig.text(position = 'TL',
             text = 'Z',
             font= '12p,blue',
            offset = '0.2c/-0.2c'
            )
    


    fig.shift_origin(yshift= -yshift)

    fig.plot(
        x = tt,
        y = rn,
        pen = '1p,red',
        frame = ['nWes', 'xa500f100', 'ya1+l"Normalised Amplitude"']
    )
    plot_tl_tr()
    fig.text(position = 'TL',
             text = 'R',
             font= '12p,red',
            offset = '0.2c/-0.2c'
            )


    fig.shift_origin(yshift= -yshift)

    fig.plot(
        x = tt,
        y = zn,
        pen = '1p,blue'
    )
    fig.plot(
        x = tt,
        y = rgrad * -1,
        pen = '1p,red,-',
        frame =  ['nWes', 'xa500f100', 'ya1']
    )
    fig.plot(
        x = tt,
        y = rhilb,
        pen = '0.5p,magenta,-',
        frame =  ['nWes', 'xa500f100', 'ya1']
    )
    plot_tl_tr()
    fig.text(position = 'TL',
             text = 'dR/dt',
             font= '12p,red',
            offset = '0.2c/-0.2c'
            )
    fig.text(position = 'BL',
             text = 'Z',
             font= '12p,blue',
            offset = '0.2c/0.2c'
            )



    
    fig.shift_origin(yshift= -yshift)

    zrg = zn - rgrad

    fig.plot(
        x = tt,
        y = zrg,
        region = [500,6000,-2.05,2.05],
        pen = '1p,magenta',
        frame = ['nWeS', 'xa500f100+l"Time after event origin (s)"', 'ya2']
    )
    plot_tl_tr()
    with pygmt.config(PS_CHAR_ENCODING = 'ISOLatin1+'):
        fig.text(position = 'TL',
                 text = 'Z + dR/dt',
                 font= '12p,magenta',
                offset = '0.2c/-0.2c'
                )


    # x-corr

    r = np.correlate(zrg, tn, mode='full')
    gr = np.gradient(abs(r))
    ind = []
    ii=0
    for jr in range(len(gr) - 1):
        if gr[jr]>0 and gr[jr + 1] < 0:
            ii +=1
            ind.append(jr)

    ind2 = []
    ii=0
    for jr in range(1,len(ind) - 1):
        if (abs(r[ind[jr]]) > abs(r[ind[jr-1]])) and (abs(r[ind[jr]]) > abs(r[ind[jr+1]])):
            ii+=1
            ind2.append(jr)

    inda = np.array(ind)
    qli=None
    for i in inda[ind2]:
        if i > len(zrg):
            qli = i
            break
    if qli is None:
        print('no QL wave')
        # with open(logf, 'a') as f:
        #     f.write(f'{evnm}_{station.code} NO\n')
        return
        
    qldp = qli - len(zrg)
    qldpt = qldp * delta
    qldpt
    qlt = qldpt + tt[G1pos]
    
    qlpos = qldp + G1pos
    print(qli, qldp,qldpt, qlpos)

    if r[qli] > 0:
        # positive corr
        fig.plot(x = tt + qldpt, 
                 y = tn,
                 pen = '1p,darkgreen'
                )
    else:
        fig.plot(x = tt + qldpt, 
                 y = -tn,
                 pen = '1p,darkgreen'
                )

    fig.plot(x = [qlt - 100, qlt - 100, qlt + 100, qlt + 100],
            y = [-10,10,10,-10],
            close = True,
            color = 'grey@70',
            pen = 'thinnest,grey'
    )
    fig.plot(x = [tt[G1pos] + qldpt, tt[G1pos] + qldpt],
            y = [-2.05, 2.05],
            pen = '0.5p,darkgreen,-')

    fig.text(position = 'BL',
             text = f'T time shifted by @~\144@~t ({qldpt} s)',
             font= '12p,darkgreen',
            offset = '0.2c/0.2c'
            )

    fig.shift_origin(yshift= -yshift - 1)
    lags = (np.arange(len(r)) - len(zrg)) * delta 
    fig.plot( 
        x = lags, 
        y = abs(r),
        pen = '1p,black',
        region = [500 - tt[G1pos], 6000 - tt[G1pos], 0, max(abs(r)) * 1.05],
        frame = ['nWeS', 'xa500f100+l"Time lag relative to t@-L@- (s)"', 'ya2000f1000+l"X-corr ampl."']
    )
    fig.plot(
        x = lags[ind],
        y = abs(r[ind]),
        pen = 'thinnest,blue'
    )
    fig.plot(
        x = lags[inda[ind2]],
        y = abs(r[inda[ind2]]),
        style = 'c0.1c',
        pen='0.5p,lightblue'
    )
    fig.plot(
        x = [0, 0],
        y = [0, 100000],
        pen = '1p,darkgreen,-'
    )
    fig.plot(
        x = [lags[qli], lags[qli]],
        y = [0, 100000],
        pen = '0.5p,darkgreen,-'
    )
    
    
    

    fig.plot(
        x = rn[int(qlpos-100):int(qlpos+100)],
        y = zn[int(qlpos-100):int(qlpos+100)],
        projection = 'X4/4',
        region = [-1,1,-1,1],
        frame = True,
        xshift = xsize + .5,
        pen = '1p,black'
    )
    
    zwind = zo[int(qlpos-100):int(qlpos+100)]
    rwind = ro[int(qlpos-100):int(qlpos+100)]
    minx = min([min(rwind), min(zwind)])
    maxx = max([max(rwind), max(zwind)])
    fig.plot(
        x = rwind,
        y = zwind,
        projection = 'X4/4',
        region = [minx,maxx,minx,maxx],
        frame = True,
        yshift = 5,
        pen = '1p,black'
    )
    
    fig.show(method = 'external')
    
    return

    fig.shift_origin(yshift = 6)

    # calc position

    G1time = time[pos]
    R1time = time[posR]

    dx = qldpt * dist / (R1time - G1time)

    fig.coast(projection = f'E{stlo}/{stla}/10', frame = True,
                  region = 'g',
                  shorelines = 'thinnest',
                 land = 'grey@50',
                  resolution = 'l', area_thresh = 100000
                 )
    fig.plot(x = stlo,
             y = stla,
             style = 't0.2c',
             color = 'red',
             pen = 'thinnest'
            )
    fig.plot(x = evlo,
             y = evla,
             style = 'a0.2c',
             color = 'red',
             pen = 'thinnest'
            )
    fig.plot(x = [stlo, evlo],
             y = [stla, evla],
             pen = '0.5p,red'
            )

    from obspy.geodetics import degrees2kilometers
    G1time = time[pos]
    R1time = time[posR]
    
    if R1time < G1time or R1time == G1time:
        print('no good')
        with open(logf, 'a') as f:
            f.write(f'{evnm}_{station.code} ERROR\n')
        return

    dx = qldpt * dist / (R1time - G1time)
    dxk = degrees2kilometers(dx)

    dxx = dxk * np.sin(baz * np.pi/180)
    dxy = dxk * np.cos(baz * np.pi/180)
    print(dx, dxk, dxx, dxy)

    with open('tmp/tmpf.tmp', 'w') as tmpf:
        tmpf.write(f'{dxx} {dxy}\n')
        
    # better fix this!
    # !gmt mapproject tmp/tmpf.tmp -JE{stlo}/{stla}/10 -Rg -I -C -Fk > tmp/lola.tmp

    with open('tmp/lola.tmp', 'r') as f:
        loc = f.read().split()
        sclo = float(loc[0])  
        scla = float(loc[1])
        print(sclo, scla)

    fig.plot('tmp/lola.tmp',
                style = 'c0.3c',
                pen= '1p,purple',
                projection = f'E{stlo}/{stla}/10',
                region = 'g'
               )
    fig.show(width = 1000)

    qual = input('Keep? [1 for yes, 0 for no, m for maybe]\n')
    print(f'Quality is {qual}')

    while qual not in ['0','1', 'm']:
        qual = input('Try again, Keep? [1 for yes, 0 for no, m for maybe]\n')
        print(f'Quality is {qual}')
    
    if qual == 'm':
        qual = 'MAYBE'
    # save results/figure

    resf = f'{outdir}/{evnm}_{station.code}.result'

    hstr = 'EVTIME EVLO EVLA EVDP MAG STLO STLA DIST BAZ SNR G1 R1 DTQL DX SCLO SCLA AMPL QUAL\n'
    outstr = f'{evo} {evlo:.3f} {evla:.3f} {evdp:.1f} {mag:.1f} ' + \
            f'{stlo:.3f} {stla:.3f} {dist:.1f} {baz:.1f} {snr:.0f} ' + \
            f'{G1time:.1f} {R1time:.1f} {qldpt:.1f} {dx:.2f} {sclo:.3f} {scla:.3f} {r[qli]:.1f} {qual}\n'

    with open(resf, 'w') as f:
        f.write(hstr)
        f.write(outstr)
    with open(logf, 'a') as f:
        f.write(f'{evnm}_{station.code} {qual}\n')

    if qual == '1':
        fig.text(position = 'TC',
                offset = '0c/3c',
                 text = 'ACCEPTED',
                 font = '24p,Helvetica-Bold,green',
                 no_clip=True
                )
    elif qual == '0':
        fig.text(position = 'TC',
                offset = '0c/3c',
                 text = 'REJECTED',
                 font = '24p,Helvetica-Bold,red',
                 no_clip=True
                )
    elif qual == 'MAYBE':
        fig.text(position = 'TC',
                offset = '0c/3c',
                 text = 'MAYBE',
                 font = '24p,Helvetica-Bold,orange',
                 no_clip=True
                )

    fig.savefig(f'{outdir}/{evnm}_{station.code}.pdf')


def ql_analysis(station,event,inv, outdir, logf):
    from obspy.geodetics import locations2degrees
    stlo = station.longitude
    stla = station.latitude
    evo = event.origins[0].time
    evnm = f'{evo.year}{evo.month:02d}{evo.day:02d}{evo.hour:02d}{evo.minute:02d}{evo.second:02d}'
    evlo = event.origins[0].longitude
    evla = event.origins[0].latitude
    
    dist = locations2degrees(evla, evlo, stla, stlo)
    print(f'downloading for {event.origins[0].time}')
    try:
        st = client.get_waveforms(network=net, station = sta, location = '*', channel = '?H?',
                         starttime = event.origins[0].time,
                         endtime = event.origins[0].time + 6000
                        )
    except:
        print('no data available')
        with open(logf, 'a') as f:
            f.write(f'{evnm}_{station.code} NO_DATA\n')
        return
            
    print(st)
    # st1 = st.select(channel='BH?', location = '00')
    # st1 = st.select(channel='VH?')
    chs = set([tr.stats.channel[0] for tr in st])
    print(chs)
    moveon=0
    if 'L' in chs:
        print('try L')
        sttmp = st.select(channel = 'LH?')
        for loc in set([tr.stats.location for tr in sttmp]):
            if len(st.select(channel = 'LH?', location = loc)) > 2:
                st1 = st.select(channel = 'LH?', location = loc)
                moveon=1
                break
    if 'B' in chs and moveon == 0:
        print('try B')
        sttmp = st.select(channel = 'BH?')
        for loc in set([tr.stats.location for tr in sttmp]):
            if len(st.select(channel = 'BH?', location = loc)) > 2:
                st1 = st.select(channel = 'BH?', location = loc)
                moveon=1
                break
    if 'H' in chs and moveon == 0:
        print('try H')
        sttmp = st.select(channel = 'HH?')
        for loc in set([tr.stats.location for tr in sttmp]):
            if len(st.select(channel = 'HH?', location = loc)) > 2:
                st1 = st.select(channel = 'HH?', location = loc)
                moveon=1
                break      
    if moveon == 0:
        print('try anything else')
        for ch in chs:
            st1 = st.select(channel = f'{ch}H?', location = '*')
            if len(st1) == 3:
                moveon = 1
                break
            elif len(st1) > 3:
                locations = set([tr.stats.location for tr in st1])
                if len(locations) > 1:
                    for loc in locations:
                        st2 = st1.select(location = loc)
                        if len(st2) > 2:
                            st1 = st2
                            moveon = 1
                            break
            if moveon == 1:
                break
        if moveon == 0:
            print('3 component data not available')
            with open(logf, 'a') as f:
                f.write(f'{evnm}_{station.code} BAD_DATA\n')
            return
            
#         st1 = st.select(channel = f'{st[0].stats.channel[0:1]}H?', location = st[0].stats.location)
    
    from obspy.geodetics import gps2dist_azimuth
    st1.detrend(type='linear')
    st1.taper(0.1)
    st1.resample(1)
    st1.filter('bandpass',freqmax=0.01,freqmin=0.002,corners=3,zerophase=True)

    evlo = event.origins[0].longitude
    evla = event.origins[0].latitude
    stlo = station.longitude
    stla = station.latitude

    evo = event.origins[0].time
    evnm = f'{evo.year}{evo.month:02d}{evo.day:02d}{evo.hour:02d}{evo.minute:02d}{evo.second:02d}'
    evd = f'{evo.year}{evo.month:02d}{evo.day:02d}'
    mag = event.magnitudes[0].mag
    evdp = event.origins[0].depth / 1000

    _, _, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
    
#     lens = [tr.stats.npts for tr in st1]
    zne= True
    for tr in st1:
        if tr.stats.channel[2] not in ['Z', 'N', 'E']:
            zne = False
    # if zne is False:
    try:
        st1.rotate('->ZNE', inventory = inv)
    except:
        print('error in components')
        with open(logf, 'a') as f:
            f.write(f'{evnm}_{station.code} ERROR\n')
        return
    try:
        st1.rotate('NE->RT', back_azimuth = baz)
    except:
        print('error rotating')
        with open(logf, 'a') as f:
            f.write(f'{evnm}_{station.code} ERROR\n')
        return
        
    
    to = st1.select(channel='?HT')[0].data
    ro = st1.select(channel='?HR')[0].data
    zo = st1.select(channel='?HZ')[0].data

    # normalise
    tn = to/max(abs(to))
    rn = ro/max(abs(ro))
    zn = zo/max(abs(zo))

    import numpy as np
    delta = st1[0].stats.delta
    # plot relative to max Love arrival
    G1max=max(abs(tn))
    pos = np.argmax(abs(tn))
    maxt = pos * delta

    R1max=max(abs(zn))
    posR=np.argmax(abs(zn))

    #  calculate SNR 
    Tnoise=max(abs(tn[int(300/delta):int(600/delta)]))
    Tsig=G1max;
    snr=Tsig/Tnoise;
    # get times
    time = delta * np.arange(len(tn))
    # calc gradient of R
    rgrad = np.gradient(rn) / max(abs(np.gradient(rn)))
    
    # plotting
    
    import pygmt
    fig = pygmt.Figure()
    fig.shift_origin(yshift= 45)

    xsize = 20
    ysize = 3
    yshift = ysize + 0.3

    fig.plot(
        projection = f'X{xsize}/{ysize}',
        region = [0,6000,-1.05,1.05],
        x = time,
        y = tn,
        pen = '2p,darkgreen',
        frame = ['nWes', 'xa500f100', 'ya1']
    )
    
    def plot_tl_tr():
        fig.plot(x = [time[pos], time[pos]],
            y = [-10, 10],
            pen = '1p,darkgreen,-')
        fig.plot(x = [time[posR], time[posR]],
            y = [-10, 10],
            pen = '1p,orange,-')
    
    plot_tl_tr()
    fig.text(position = 'TL',
             text = 'T',
             font= '12p,OpenSans-Bold,darkgreen',
            offset = '0.2c/-0.2c'
            )
    fig.text(x = time[pos],
             y = 1.3,
             text = 't@-L',
             font = '10p,Helvetica-Bold,darkgreen',
             no_clip = True
            )
    fig.text(x = time[posR],
             y = 1.3,
             text = 't@-R',
             font = '10p,Helvetica-Bold,orange',
             no_clip = True
            )

    hdrtxt1 = f'{evo}   EVLO: {evlo:.3f}   EVLA: {evla:.3f}   EVDP: {evdp:.1f} km   M: {mag:.1f}'
    hdrtxt2 = f'{net}.{sta}   STLO: {stlo:.3f}   STLA: {stla:.3f}   ' + \
            f'DIST: {dist:.1f}@+o@+   BAZ: {baz:.1f}@+o@+   SNR: {snr:.0f} '

    fig.text(position = 'TC',
             text = hdrtxt1,
             offset = '1.2c',
             no_clip = True
            )
    fig.text(position = 'TC',
             text = hdrtxt2,
             offset = '0.8c',
             no_clip = True
            )

    fig.shift_origin(yshift= -yshift)

    fig.plot(
        x = time,
        y = zn,
        pen = '1p,blue',
        frame = ['nWes', 'xa500f100', 'ya1']
    )
    plot_tl_tr()
    fig.text(position = 'TL',
             text = 'Z',
             font= '12p,blue',
            offset = '0.2c/-0.2c'
            )


    fig.shift_origin(yshift= -yshift)

    fig.plot(
        x = time,
        y = rn,
        pen = '1p,red',
        frame = ['nWes', 'xa500f100', 'ya1+l"Normalised Amplitude"']
    )
    plot_tl_tr()
    fig.text(position = 'TL',
             text = 'R',
             font= '12p,red',
            offset = '0.2c/-0.2c'
            )


    fig.shift_origin(yshift= -yshift)

    fig.plot(
        x = time,
        y = zn,
        pen = '1p,blue'
    )
    fig.plot(
        x = time,
        y = rgrad * -1,
        pen = '1p,red,-',
        frame =  ['nWes', 'xa500f100', 'ya1']
    )
    plot_tl_tr()
    fig.text(position = 'TL',
             text = 'dR/dt',
             font= '12p,red',
            offset = '0.2c/-0.2c'
            )
    fig.text(position = 'BL',
             text = 'Z',
             font= '12p,blue',
            offset = '0.2c/0.2c'
            )

    fig.shift_origin(yshift= -yshift)

    zrg = zn - rgrad

    fig.plot(
        x = time,
        y = zrg,
        region = [500,6000,-2.05,2.05],
        pen = '1p,magenta',
        frame = ['nWeS', 'xa500f100+l"Time after event origin (s)"', 'ya2']
    )
    plot_tl_tr()
    with pygmt.config(PS_CHAR_ENCODING = 'ISOLatin1+'):
        fig.text(position = 'TL',
                 text = 'Z + dR/dt',
                 font= '12p,magenta',
                offset = '0.2c/-0.2c'
                )


    # x-corr

    r = np.correlate(zrg, tn, mode='full')
    gr = np.gradient(abs(r))
    ind = []
    ii=0
    for jr in range(len(gr) - 1):
        if gr[jr]>0 and gr[jr + 1] < 0:
            ii +=1
            ind.append(jr)

    ind2 = []
    ii=0
    for jr in range(1,len(ind) - 1):
        if (abs(r[ind[jr]]) > abs(r[ind[jr-1]])) and (abs(r[ind[jr]]) > abs(r[ind[jr+1]])):
            ii+=1
            ind2.append(jr)

    inda = np.array(ind)
    qli=None
    for i in inda[ind2]:
        if i > len(zrg):
            qli = i
            break
    if qli is None:
        print('no QL wave')
        with open(logf, 'a') as f:
            f.write(f'{evnm}_{station.code} NO\n')
        return
        
    qldp = qli - len(zrg)
    qldpt = qldp * delta
    qldpt

    if r[qli] > 0:
        # positive corr
        fig.plot(x = time + qldpt, 
                 y = tn,
                 pen = '1p,darkgreen'
                )
    else:
        fig.plot(x = time + qldpt, 
                 y = -tn,
                 pen = '1p,darkgreen'
                )

    fig.plot(x = [time[pos] + qldpt, time[pos] + qldpt],
            y = [-2.05, 2.05],
            pen = '0.5p,darkgreen,-')

    fig.text(position = 'BL',
             text = f'T time shifted by @~\144@~t ({qldpt} s)',
             font= '12p,darkgreen',
            offset = '0.2c/0.2c'
            )

    fig.shift_origin(yshift= -yshift - 1)
    lags = (np.arange(len(r)) - len(zrg)) * delta 
    fig.plot( 
        x = lags, 
        y = abs(r),
        pen = '1p,black',
        region = [500 - time[pos], 6000 - time[pos], 0, max(abs(r)) * 1.05],
        frame = ['nWeS', 'xa500f100+l"Time lag relative to t@-L@- (s)"', 'ya2000f1000+l"X-corr ampl."']
    )
    fig.plot(
        x = lags[ind],
        y = abs(r[ind]),
        pen = 'thinnest,blue'
    )
    fig.plot(
        x = lags[inda[ind2]],
        y = abs(r[inda[ind2]]),
        style = 'c0.1c',
        pen='0.5p,lightblue'
    )
    fig.plot(
        x = [0, 0],
        y = [0, 100000],
        pen = '1p,darkgreen,-'
    )
    fig.plot(
        x = [lags[qli], lags[qli]],
        y = [0, 100000],
        pen = '0.5p,darkgreen,-'
    )

    fig.shift_origin(yshift = 5, xshift = xsize + .5)

    # calc position

    G1time = time[pos]
    R1time = time[posR]

    dx = qldpt * dist / (R1time - G1time)

    fig.coast(projection = f'E{stlo}/{stla}/10', frame = True,
                  region = 'g',
                  shorelines = 'thinnest',
                 land = 'grey@50',
                  resolution = 'l', area_thresh = 100000
                 )
    fig.plot(x = stlo,
             y = stla,
             style = 't0.2c',
             color = 'red',
             pen = 'thinnest'
            )
    fig.plot(x = evlo,
             y = evla,
             style = 'a0.2c',
             color = 'red',
             pen = 'thinnest'
            )
    fig.plot(x = [stlo, evlo],
             y = [stla, evla],
             pen = '0.5p,red'
            )

    from obspy.geodetics import degrees2kilometers
    G1time = time[pos]
    R1time = time[posR]
    
    if R1time < G1time or R1time == G1time:
        print('no good')
        with open(logf, 'a') as f:
            f.write(f'{evnm}_{station.code} ERROR\n')
        return

    dx = qldpt * dist / (R1time - G1time)
    dxk = degrees2kilometers(dx)

    dxx = dxk * np.sin(baz * np.pi/180)
    dxy = dxk * np.cos(baz * np.pi/180)
    print(dx, dxk, dxx, dxy)

    with open('tmp/tmpf.tmp', 'w') as tmpf:
        tmpf.write(f'{dxx} {dxy}\n')
        
    # better fix this!
    # !gmt mapproject tmp/tmpf.tmp -JE{stlo}/{stla}/10 -Rg -I -C -Fk > tmp/lola.tmp

    with open('tmp/lola.tmp', 'r') as f:
        loc = f.read().split()
        sclo = float(loc[0])  
        scla = float(loc[1])
        print(sclo, scla)

    fig.plot('tmp/lola.tmp',
                style = 'c0.3c',
                pen= '1p,purple',
                projection = f'E{stlo}/{stla}/10',
                region = 'g'
               )
    fig.show(width = 1000)

    qual = input('Keep? [1 for yes, 0 for no, m for maybe]\n')
    print(f'Quality is {qual}')

    while qual not in ['0','1', 'm']:
        qual = input('Try again, Keep? [1 for yes, 0 for no, m for maybe]\n')
        print(f'Quality is {qual}')
    
    if qual == 'm':
        qual = 'MAYBE'
    # save results/figure

    resf = f'{outdir}/{evnm}_{station.code}.result'

    hstr = 'EVTIME EVLO EVLA EVDP MAG STLO STLA DIST BAZ SNR G1 R1 DTQL DX SCLO SCLA AMPL QUAL\n'
    outstr = f'{evo} {evlo:.3f} {evla:.3f} {evdp:.1f} {mag:.1f} ' + \
            f'{stlo:.3f} {stla:.3f} {dist:.1f} {baz:.1f} {snr:.0f} ' + \
            f'{G1time:.1f} {R1time:.1f} {qldpt:.1f} {dx:.2f} {sclo:.3f} {scla:.3f} {r[qli]:.1f} {qual}\n'

    with open(resf, 'w') as f:
        f.write(hstr)
        f.write(outstr)
    with open(logf, 'a') as f:
        f.write(f'{evnm}_{station.code} {qual}\n')

    if qual == '1':
        fig.text(position = 'TC',
                offset = '0c/3c',
                 text = 'ACCEPTED',
                 font = '24p,Helvetica-Bold,green',
                 no_clip=True
                )
    elif qual == '0':
        fig.text(position = 'TC',
                offset = '0c/3c',
                 text = 'REJECTED',
                 font = '24p,Helvetica-Bold,red',
                 no_clip=True
                )
    elif qual == 'MAYBE':
        fig.text(position = 'TC',
                offset = '0c/3c',
                 text = 'MAYBE',
                 font = '24p,Helvetica-Bold,orange',
                 no_clip=True
                )

    fig.savefig(f'{outdir}/{evnm}_{station.code}.pdf')