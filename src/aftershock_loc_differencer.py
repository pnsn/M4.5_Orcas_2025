"""
:module: M4.5_Orcas_2025/src/aftershock_loc_differencer.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: A set of scripts for parsing a plain-text paste of a Jiggle origin/magnitude solution string
  and getting time and distance differences between a given event (ideally an aftershock) and the 
  M4.5 mainshock at Orcas Island on March 3rd 2025.

"""


from pathlib import Path
import pandas as pd
from obspy import UTCDateTime
from obspy.geodetics import locations2degrees

eg_hdr = '61503963 UW 2025-03-03 21:14:40.17   48.6057 -122.8103  17.70  0.92 Ml  eq  L amyw     UW 01  H   2  -  H P5   17.52  0.18'


LAT, LON, MZ = 48.611, -122.805, 15.97
T0 = UTCDateTime('2025-03-03T13:02:37.850')
def parse_jiggle_origin_header(hdr):
    # Split & remove extra whitespace
    parts = hdr.split(' ')
    parts = [_e for _e in parts if _e != '']

    # Parse lines
    evid = parts[0]
    t0 = UTCDateTime(f'{parts[2]}T{parts[3]}')
    lon = float(parts[4])
    lat = float(parts[5])
    mz = float(parts[6])
    mag = float(parts[7])
    magtype = parts[8]
    etype = parts[9]
    output = [evid, t0, lon, lat, mz, mag, magtype, etype]
    output = dict(zip(['evid','t0','lon','lat','mz','mag','magtype','etype'], output))
    return output

def get_distances_hdr(hdr):
    after = parse_jiggle_origin_header(hdr)
    dist = locations2degrees(LAT,LON, after['lat'], after['lon'])
    dist_km = 111.2/dist
    dist_3d_km = (dist_km**2 + (after['mz'] - MZ)**2)**0.5
    dt = after['t0'] - T0
    print(f'Epicentral distance: {dist_km: .3f} km')
    print(f'Hypocentral distance {dist_3d_km:.3f} km')
    print(f'Origin time {dt/3600:.3f} hrs after mainshock')
    print(f'Origin time {dt/60:.3f} min after mainshock')


def get_distances(df):
    df = df.sort_values('MAG', ascending=False)
    ser_main = df.loc[62078906]
    holder = []
    for evid, row in df.iterrows():
        # Get epicentral offsets
        dh_deg = locations2degrees(ser_main.LAT, ser_main.LON, row.LAT, row.LON)
        dh_km = 111.2*dh_deg
        dz_km = row.MZ - ser_main.MZ
        dx_km = (dz_km**2 + dh_km**2)**0.5
        dt_sec = (row.DATETIME - ser_main.DATETIME).total_seconds()
        holder.append([dx_km, dh_km, dt_sec])
    df_out = pd.DataFrame(data=holder, index=df.index, columns=['hyp_off_km','epi_off_km','orig_off_sec'])
    return df_out


if __name__ == '__main__':
    ROOT = Path(__file__).parent.parent
    EFILE = ROOT / 'data' / 'jiggle' / 'Event_Table_Output_4MAR2025_1900UTC.csv'
    # Load events with datetimes
    df = pd.read_csv(EFILE, index_col=[0], parse_dates=['DATETIME'])
    # Sort by descending magnitude
    df = df.sort_values(by='MAG', ascending=False)
    # Calculate time/location differences
    _df = get_distances(df)
    # Attach time/location differences
    df = pd.concat([df,_df], axis=1, ignore_index=False)