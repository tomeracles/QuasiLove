#!/usr/bin/env python
# Script to download waveforms for QL analysis
# Will download waveforms for M > 6.5, depth <50km, dist > 70deg events from NEIC PDE catalogue
# and for a given network and station(s).
#
# Usage: ./scripts/download_waveforms.py --network <net> [--stations <sta>] --outdir ./<dir>/waveforms [--minbaz <minbaz>] [--maxbaz <maxbaz>]
#

import obspy
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
iris = Client("IRIS")

# These are the newtworks I was using and their associated data clients.
networks = { 'IRIS': [ 'II', 'IU', 'IM', 'TT',  'AF', 'G', 'GT', 'IP', 'NJ', 'GH',  'NR', 'BX', '4J', 'WM', 'AU'],
            'GFZ': ['GE', ],
            'ORFEUS' : ['DZ'],
            'RESIF' : ['QM'],
            'ICGC' : ['ES'],
            'INGV' : ['MN'],
            'NOA' : ['HL']
}

eqcatfile = './data/NEIC_PDE_M6.5_1980-20230206.xml'
# get earthquake catalogue: M > 6.5
def get_earthquake_catalogue(out = eqcatfile):
    # This was done once to download the catalogue
    starttime = UTCDateTime('1980-01-01T')
    endtime = UTCDateTime('2024-12-12T')
    cat = iris.get_events(starttime=starttime, endtime=endtime,
                            minmagnitude=6.5, catalog="NEIC PDE", maxdepth=50)
    cat.write(out, format='QUAKEML')

cat = obspy.read_events(eqcatfile)

def get_inv(net,sta):
    from obspy import Inventory
    from obspy.clients.fdsn import Client
    inv = Inventory()
    dmc = [d for d in networks if net in networks[d]][0]
    print(dmc)
    client = Client(dmc)
    inv.extend(client.get_stations(network=net, station=sta, level = 'channel'))
    return inv, client


def download(station,event,client,inv, outdir = './Africa/waveforms'):
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
    st1.filter('bandpass',freqmax=0.01,freqmin=0.001,corners=3,zerophase=True)

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
        sac.o = evo
        sac.iztype = "io"
        sac.evlo = evlo
        sac.evla = evla
        sac.evdp = evdp
        sac.stlo = stlo
        sac.stla = stla
        sac.mag = mag
        dest = f'{outdir}/{evnm}_{sac.knetwk}.{sac.kstnm}.{sac.kcmpnm}'
        sac.write(dest)

def main(outdir = './Africa/waveforms', network = 'GE', stations = 'CSS,KARP',
    minbaz = 0, maxbaz = 360):
    from obspy.geodetics import locations2degrees
    
    inv, client = get_inv(network, stations)
    for sta in inv[0]:
        sta.network = inv[0].code
        stlo = sta.longitude
        stla = sta.latitude
        print(sta.code)
        for event in cat:
            if sta.start_date > event.origins[0].time or sta.end_date < event.origins[0].time:
                continue
            evo = event.origins[0].time
            evnm = f'{evo.year}{evo.month:02d}{evo.day:02d}{evo.hour:02d}{evo.minute:02d}{evo.second:02d}'
            evlo = event.origins[0].longitude
            evla = event.origins[0].latitude
            dist = locations2degrees(evla, evlo, stla, stlo)
            _, az, baz = gps2dist_azimuth(evla, evlo, stla, stlo)

            if dist > 70 and baz > minbaz and baz < maxbaz:
                print(event.origins[0].time)
                download(sta,event,client,inv,outdir = outdir)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Download waveforms for QL analysis.')
    parser.add_argument('--outdir', type=str, default='./Africa/waveforms',
                        help='Output directory for waveforms.')
    parser.add_argument('--network', type=str, 
                        help='Network code to download from.')
    parser.add_argument('--stations', type=str, default='*',
                        help='Station code, or comma-separated list of station codes, or * for all available stations, to download.')
    parser.add_argument('--minbaz', type=float, default=0,
                        help='Minimum back azimuth for station-event pairs.')
    parser.add_argument('--maxbaz', type=float, default=360,
                        help='Maximum back azimuth for station-event pairs.')
    args = parser.parse_args()
    main(outdir=args.outdir, network=args.network, stations=args.stations)