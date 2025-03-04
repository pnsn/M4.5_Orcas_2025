import os
from pathlib import Path
from obspy import *
from obspy.clients.fdsn import Client

from obsplus import WaveBank

from pyrocko import obspy_compat
from pyrocko.gui.snuffler.marker import save_markers, load_markers

obspy_compat.plant()

ROOT = Path(__file__).parent.parent
SAVEPATH = ROOT/'data'/'snuffler'
ALL_MARKERS = ROOT/'data'/'snuffler'/'_markers_working.dat'
try:
    os.makedirs(str(SAVEPATH), exist_ok=False)
except:
    pass

IRIS = Client('IRIS')

# t0 = UTCDateTime('2025-03-03T13:02:37') + 3600*5.5
t0 = UTCDateTime('2025-03-04T15:00:00')
tnow = UTCDateTime()
STAS = 'OLGA,MCW,ORCA,GUEM,LUMI,TURTL,LOPEZ,OAKH,SAXON,DONK'
CHAS='HH?,BH?,EH?,HN?'
NETS = 'UW'
inv = IRIS.get_stations(network=NETS, station=STAS)
# cat = IRIS.get_events(starttime = t0 - 10, endtime=t0 + 10, latitude=48.607, longitude=-122.804, maxradius=0.1/111.2)
st = IRIS.get_waveforms(network=NETS, station=STAS, location='*', channel=CHAS, starttime=t0, endtime=tnow)

# TODO - Load running marker file
tmp_mrkr = load_markers(str(ALL_MARKERS))
return_tag, markers = st.snuffle(ntracks=len(st), markers=tmp_mrkr)


# breakpoint()
# Only save new incremental file if there are new markers
if len(markers) > 0:
    savename = SAVEPATH/f'markers_{tnow}.dat'
    save_markers(markers, str(savename))

# Always save running marker file
save_markers(markers, str(ALL_MARKERS))


