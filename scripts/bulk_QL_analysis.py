#!/usr/bin/env python
import sys
sys.path.append('./QuasiLove')
from QL_detection import ql_analysis_from_sac
from QL_detection import SAC
from glob import glob

sac = SAC('/usr/local/LINUX_sac/bin/sacinit.sh') # Here you should put the path to your favourite SAC executable or initialisation script

indir = './Africa/waveforms'
outdir = './Africa/results'
logf = outdir + '/logfile.txt'

flist = sorted(list(set([f[0:-1] for f in glob(f'{indir}/??????????????_*.*.???')])))
flist.reverse()

fstem = flist[0]
ql_analysis_from_sac(fstem, outdir, logf, sac)