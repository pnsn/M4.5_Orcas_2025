"""
:module: M4.5_Orcas_2025/src/snuffle_aftershocks.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: This script queries a specified period of waveform data from
    stations near the M4.5 earthquake at Orcas Island on March 3rd 2025
    and launches a Pyrocko Snuffler instance to display those waveforms
    for analyst picking. Picks are automatically saved to the `data/snuffler`
    directory and in subsequent sessions called by this script, the
    `data/snuffler/_markers_working.dat` file is re-loaded into snuffler
    to load prior picking progress.
"""


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
ALL_MARKERS = ROOT/'data'/'snuffler'/'_markers_working_62079456.dat'
try:
    os.makedirs(str(SAVEPATH), exist_ok=False)
except:
    pass

IRIS = Client('IRIS')

# t0 = UTCDateTime('2025-03-03T13:02:37') + 3600*5.5
t0 = UTCDateTime('2025-03-05T20:00:00')
tnow = UTCDateTime()
STAS = 'CHIMA,DOSE,LRIV,OSQM'
#'OLGA,MCW,ORCA,GUEM,LUMI,TURTL,LOPEZ,OAKH,SAXON,DONK'
CHAS='HH?,BH?,EH?,HN?'
NETS = 'UW'
inv = IRIS.get_stations(network=NETS, station=STAS)
# cat = IRIS.get_events(starttime = t0 - 10, endtime=t0 + 10, latitude=48.607, longitude=-122.804, maxradius=0.1/111.2)
st = IRIS.get_waveforms(network=NETS, station=STAS, location='*', channel=CHAS, starttime=t0, endtime=tnow)

# TODO - Load running marker file
try:
    tmp_mrkr = load_markers(str(ALL_MARKERS))
except:
    tmp_mrkr = None
return_tag, markers = st.snuffle(ntracks=len(st), markers=tmp_mrkr)


# breakpoint()
# Only save new incremental file if there are new markers
if len(markers) > 0:
    savename = SAVEPATH/f'markers_{tnow}.dat'
    save_markers(markers, str(savename))

# Always save running marker file
save_markers(markers, str(ALL_MARKERS))


