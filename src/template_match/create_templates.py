from pathlib import Path

import pandas as pd

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.event import (
    Catalog, Event, Origin, Arrival, Magnitude,
    ResourceIdentifier, QuantityError, OriginQuality,
    OriginUncertainty)

from eqcorrscan import Tribe

from eqcutil import ClusteringTribe
from eqcutil.catalog.model_phases import model_picks

def aqms2cat(df, inv, phases=['P'], velocity_model='P4'):
    """
    Convert a dataframe representation of an AQMS event table
    exported from Jiggle into a :class:`~obspy.core.event.Catalog`
    object and populate events with modeled pick times for stations
    included in **inv**
    """
    cat = Catalog()
    # Iterate across events
    for evid, row in df.iterrows():
        #Convert row into event, origin, and magnitude
        event = Event(
            resource_id = ResourceIdentifier(prefix='smi:local/Event/'),
            event_type=row.ETYPE[1:-1])
        origin = Origin(
            resource_id=ResourceIdentifier(id=f'quakeml:uw.anss.org/Origin/UW/{evid}'),
            time=UTCDateTime(str(row.DATETIME)),
            longitude=row.LON,
            latitude=row.LAT,
            depth=row.MZ,
            depth_type = {'0':'from location','1':'operator assigned'}[row.ZF[1:-1]],
            depth_errors=QuantityError(uncertainty=row.ERR_Z*1e3),
            epicenter_fixed=row.HF=="'1'",
            time_fixed=row.TF=="'1'",
            quality=OriginQuality(associated_phase_count=row.OBS,
                                  used_phase_count=row.USED,
                                  minimum_distance=row.DIST/111.2,
                                  azimuthal_gap=row.GAP),
            origin_uncertainty=OriginUncertainty(horizontal_uncertainty=row.ERR_H*1e3),
            evaluation_mode={"'F'":'manual',"'I'":'manual', "'A'":'automatic', "'H'":'automatic'}[row.ST],
            evaluation_status={"'F'":'final',"'I'":'reviewed',"'H'":'confirmed',"'A'":'preliminary'}[row.ST],
        )

        magnitude = Magnitude(mag=row.MAG,
                              magnitude_type=row.MTYP,
                              origin_id=origin.resource_id,
                              mag_errors=QuantityError(uncertainty=row.MERR))

        # Model Picks
        picks = model_picks(origin, inv, model_name=velocity_model, phases=phases)
        # Associate picks to origin & attach picks to event
        for _p in picks:
            arr = Arrival(
                pick_id=_p.resource_id,
                phase=_p.phase_hint,
                earth_model_id=ResourceIdentifier(
                    id=f'quakeml:uw.anss.org/VelocityModel/UW/{velocity_model}'
                    )
                )
            event.picks.append(_p)
            origin.arrivals.append(arr)
        
        # Attach origin to event
        event.origins.append(origin)
        # Attach magnitude to event
        event.magnitudes.append(magnitude)
        # Set origin as preferred origin
        event.preferred_origin = origin.resource_id
        # Set magnitude as preferred magnitude
        event.preferred_magnitude = magnitude.resource_id

        # Attach to catalog
        cat.events.append(event)
    return cat
        
if __name__ == '__main__':
    ROOT = Path(__file__).parent.parent.parent
    AQMS_DATA = ROOT / 'data' / 'jiggle' / 'Event_Table_Output_4MAR2025_1900UTC.csv'

    STAS = 'OLGA,TURTL'
    NETS = 'UW'
    CHANS = 'HHZ'

    tckwargs = {
        'method': 'from_client',
        'lowcut': 5.,
        'highcut': None,
        'filt_order': 4,
        'samp_rate': 100.,
        'length': 6.,
        'prepick': 1.,
        'process_length': 3600.,
        'min_snr': 2.,
        'num_cores': 4,
        'save_progress': False
    }

    xcckwargs = {
        'method': 'xcc',
        'replace_nan_distances_with': 1,
        'shift_len': 2,
        'corr_thresh': 0.6,
        'allow_individual_trace_shifts': True,
        'cores': 'all'
    }


    # PROCESSING SECTION #

    # Connect to webservices
    IRIS = Client('IRIS')
    # Attach client to construction kwargs
    tckwargs.update({'client_id': IRIS})

    # Read event table
    df = pd.read_csv(AQMS_DATA, index_col=[0], parse_dates=['DATETIME'])
    # Strip off mainshock
    adf = df #.iloc[1:]
    # Get station inventory
    inv = IRIS.get_stations(
        station=STAS,
        network=NETS,
        channel=CHANS,
        level='channel')
    # Convert event table into catalog 
    cat = aqms2cat(adf, inv)

    # Attach catalog to tckwargs
    tckwargs.update({'catalog': cat})

    # Construct templates
    tribe = Tribe().construct(**tckwargs)
    # Rename templates
    for tmp in tribe:
        tmp.name = ''.join(str(tmp.event.preferred_origin).split('/')[-2:])
    # Convert tribe into clusteringtribe
    ctr = ClusteringTribe(templates=tribe.templates)
    # Run template xcorr clustering
    ctr.cluster(**xcckwargs)



