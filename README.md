# M4.5_ORCAS_2025  
Visual review tools for identifying aftershocks of the M4.5 10 km east of Orcas Island on March 3rd 2025.  

# Description  
Scripts for pulling waveform data from the IRIS Data Management Center webservices
and viewing/reviewing them in Pyrocko's "Snuffler" application. Data are saved as
timestamped snuffler marker files in `data/snuffler_markers/markers_{timestamp}` for each review session
launched by the `src/snuffle_aftershocks.py` script and all markers are saved in a working
data file `data/snuffler_markers_working.dat` that is loaded at the start of a review session
and overwritten at the end of each review session.

# Author
Nathan T. Stevens (ntsteven@uw.edu)

# License
GNU General Public License v3 (GPLv3)