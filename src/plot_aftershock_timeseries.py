import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.offsetbox import AnchoredText

# from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees
from obspy.imaging.beachball import beach
import pandas as pd


import cartopy.crs as ccrs
from cartopy.io.img_tiles import OSM
import cartopy.feature as cfeature

# Define absolute path to repository root
ROOT = Path(__file__).parent.parent
# Save path for writing figure files
FIGPATH = ROOT / 'figures'
# Path to PNSN Logo PNG
LOGO_PNG = ROOT / 'data' / 'resources' / 'PNSN_Small_RGB.png'
# Path to AQMS Event Table CSV output from Jiggle
AQMS_CSV = ROOT / 'data' / 'jiggle' / 'Event_Table_Output_4MAR2025_1900UTC.csv'

# Set figure saving/resolution controls
issave = True
DPI = 250
FMT = 'png'


# Set if figure should be plt.show()'d
isshow = True

UTM10N = ccrs.UTM(zone=10, southern_hemisphere=False)
WGS84 = ccrs.PlateCarree()
ZOOM = 11
RAD_KM = 12.5

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

def rad2llur(rlat, rlon, rad_m=50000.):
    
    # Convert reference location to northing & easting
    mEo, mNo = UTM10N.transform_point(rlon,rlat,WGS84)
    llE = mEo - rad_m
    llN = mNo - rad_m
    urE = mEo + rad_m
    urN = mNo + rad_m
    lowerleft = WGS84.transform_point(llE, llN, UTM10N)
    upperright = WGS84.transform_point(urE, urN, UTM10N)
    
    return [lowerleft[0], upperright[0], lowerleft[1], upperright[1]]


# Load event data
df = pd.read_csv(AQMS_CSV, index_col=[0], parse_dates=['DATETIME'])
# Add distances
df = pd.concat([df, get_distances(df)], axis=1, ignore_index=False)
# Get series of mainshock statistics
ser_main = df.loc[62078906]



# Create connection to OSM tiles
imagery = OSM(cache=True)

# Specify moment tensor for beachball
np1 = [15., 55., 90.]
# xy = UTM10N.transform_point(x=ser_main.LON, y=ser_main.LAT, src_crs=WGS84)
x, y = imagery.crs.transform_point(x=ser_main.LON, y=ser_main.LAT, src_crs=ccrs.Geodetic())
bb_main = beach(np1, xy=(x,y), width=6000, zorder=1)

# Initialize Figure
fig = plt.figure(figsize=(5.6,7.7))
gs = fig.add_gridspec(ncols=1, nrows=3, hspace=0)
# Create axis for map
axmap = fig.add_subplot(gs[:2], projection=imagery.crs)
axmap.set_extent(rad2llur(df.LAT.median(), df.LON.median(), rad_m=RAD_KM*1e3), ccrs.PlateCarree())
# Add basemap imagery
axmap.add_image(imagery, ZOOM)
text = AnchoredText('©OpenStreetMap contributors',
                    loc=4, prop={'size': 6}, frameon=True)
# Add attribution
axmap.add_artist(text)

# Plotmainshock Beachball
axmap.add_collection(bb_main)

# Plot Events
adf = df[df.index.values != ser_main.name]
axmap.scatter(adf.LON, adf.LAT, s=3**(2 + adf.MAG),
              c='none', edgecolors='red',
              alpha=0.75, transform=ccrs.PlateCarree(),
              linewidths=2)

# Add legend
# x, y = imagery.crs.transform_point(x=48.67, y=-122.75, src_crs=ccrs.Geodetic())
# axmap.add_patch(mpatches.Rectangle(xy=[x,y], width=0.1, height=0.03,transform=ccrs.PlateCarree()))

# Plot Stations
# axmap.plot(df_inv.lon, df_inv.lat, )

# Add coordinates
gl = axmap.gridlines(draw_labels=True, zorder=1)
gl.bottom_labels=False
gl.left_labels=False
gl.xlines=False
gl.ylines=False


# MAGNITUDE TIME-SERIES
axts = fig.add_subplot(gs[-1])
# Do Mainshock stem
ml, sl, bl = axts.stem(0, ser_main.MAG, bottom=-1)
# Format marker lines
ml.set_markerfacecolor('white')
ml.set_markeredgecolor('blue')
ml.set_linewidth(2)
sl.set_color('blue')
sl.set_linewidth(2)
bl.set_color('none')

# Do Stem plot
ml, sl, bl = axts.stem(adf.orig_off_sec/3600, adf.MAG, bottom = -1)
# Format marker lines
ml.set_markerfacecolor('white')
ml.set_markeredgecolor('red')
ml.set_linewidth(2)
sl.set_color('red')
sl.set_linewidth(2)
bl.set_color('none')
# Add gridlines
axts.grid(linestyle=':')
# Scale Figure
axts.set_ylim([-1, 5])
# Label axes
axts.set_xlabel(f'Hours Since {df.DATETIME.min().strftime("%Y-%m-%d %H:%M:%S")} (UTC)')
axts.set_ylabel('Magnitude')


# ADD PNSN LOGO
logoax = fig.add_axes([0.01, 0.9, 0.3, 0.3], anchor='SE', zorder=-1)
im = plt.imread(str(LOGO_PNG))
logoax.imshow(im)
logoax.axis('off')

# SAVE FIGURE (IF SWITCH IS TURNED ON)
if issave:
    try:
        os.makedirs(str(FIGPATH), exist_ok=False)
    except:
        pass
    plt.savefig(str(FIGPATH/f'Aftershock_Timeseries_{pd.Timestamp.now()}_{DPI}dpi.{FMT}'), format=FMT, dpi=DPI)

# DISPLAY FIGURE (IF SWITCH IS TURNED ON)
if isshow:
    plt.show()