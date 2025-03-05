import os, logging
from pathlib import Path

from obspy import UTCDateTime
from obspy.clients.fdsn import Client

from eqcutil import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger


if __name__ == '__main__':
    # Create logger
    Logger = setup_terminal_logger(name='run_match_filter', level=logging.INFO)

    ROOT = Path(__file__).parent.parent.parent
    CTR_FILE = ROOT / 'processed_data' / 'templates' / 'templates_for_match_filter.tgz'

    ## TEMPLATE MATCH PARAMETERIZATION SECTION ##
    T0 = UTCDateTime('2025-01-01T00:00:00')
    T1 = UTCDateTime('2025-02-01T00:00:00') #T0 + 6.*3600.
    THRESH_TYPE = 'MAD'
    THRESH = 8.     # From Shelly & Beroza (2007)
    TRIG_INT = 1.   # 1 second gap between detections minimum
    DAYLONG = False  # Processing 24 hour chunks?
    PARPROC = True  # Processing in parallel
    SAVEPROGRESS = True
    NCORES = 12       
    RETURN_STREAM = True

    ## PROCESSING SECTION ##

    # Connect to client
    IRIS = Client('IRIS')
    Logger.info(f'Connected to client')
    # Load templates
    ctr = ClusteringTribe().read(str(CTR_FILE))
    Logger.info(f'Loaded {len(ctr)} templates')
    # Run template matching
    outs = ctr.client_detect(
        client = IRIS,
        starttime = T0,
        endtime = T1,
        threshold = THRESH,
        threshold_type = THRESH_TYPE,
        trig_int = TRIG_INT,
        daylong = DAYLONG,
        concurrent_processing=False,
        parallel_process=PARPROC,
        save_progress=SAVEPROGRESS,
        process_cores=NCORES,
        return_stream=RETURN_STREAM
    )

    if RETURN_STREAM:
        party, full_st = outs

    else:
        party = outs
     
    
