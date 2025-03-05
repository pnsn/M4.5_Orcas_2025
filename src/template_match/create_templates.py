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

### SUPPORTING METHOD FOR CONVERTING AQMS EVENT CSV INTO OBSPY CATALOG

def aqms2cat(df, inv, phases=['P','p'], velocity_model='P4', pick_preference='earliest'):
    """
    Convert a dataframe representation of an AQMS event table
    exported from Jiggle into a :class:`~obspy.core.event.Catalog`
    object and populate events with modeled pick times for stations
    included in **inv**.
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
        
        nslc_pick_times = {}
        for _p in picks:
            nslc = _p.waveform_id.id
            if nslc in nslc_pick_times.keys():
                if pick_preference == 'earliest':
                    if nslc_pick_times[nslc] < _p.time:
                        continue
                    else:
                        nslc_pick_times[nslc] = _p.time
                elif pick_preference == 'latest':
                    if nslc_pick_times[nslc] > _p.time:
                        continue
                    else:
                        nslc_pick_times[nslc] = _p.time
                elif pick_preference in ['all', None]:
                    pass
            else:
                nslc_pick_times.update({nslc: _p.time})

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
        
### MAIN CODE ###

if __name__ == '__main__':
    ROOT = Path(__file__).parent.parent.parent
    AQMS_DATA = ROOT / 'data' / 'jiggle' / 'Event_Table_Output_4MAR2025_1900UTC.csv'

    OUTPUT_DIR = ROOT / 'processed_data' / 'templates'

    # DEFINE INVENTORY
    STAS = 'OLGA,TURTL'
    NETS = 'UW'
    CHANS = 'HHZ'

    # Production station delays from the 2023 P5 station delay set (P5.del)
    OFFICIAL_STATION_DELAYS = {'UW.OLGA..HHZ': 0.06, 'UW.TURTL..HHZ': 0.01}
    # Additional delays applied based on visual review of earlier template generation
    AD_HOC_STATION_DELAYS = {'UW.OLGA..HHZ': 2.15, 'UW.TURTL..HHZ': 0.5 }
    # Minimum mean Signal to Noise Ratio for template matching
    TEMPLATE_SNR_MIN = 5.

    tckwargs = {
        'method': 'from_client',
        'lowcut': 5.,
        'highcut': None,
        'filt_order': 4,
        'samp_rate': 100.,
        'length': 8.,
        'prepick': 1.5,
        'process_length': 3600.,
        'min_snr': 1.3,
        'num_cores': 4,
        'save_progress': False
    }

    xcckwargs = {
        'method': 'xcc',
        'replace_nan_distances_with': 'mean',
        'shift_len': 4,
        'corr_thresh': 0.45,
        'allow_individual_trace_shifts': False,
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
    
    # Convert event table into catalog & model arrival times
    cat = aqms2cat(adf, inv, pick_preference='earliest')
    # Manually apply station delays based on model P5 2023 P-wave station delays

    # Apply station delays
    for event in cat.events:
        for pick in event.picks:
            # Get composite station delay correction
            dt = OFFICIAL_STATION_DELAYS[pick.waveform_id.id] +\
                 AD_HOC_STATION_DELAYS[pick.waveform_id.id]
            # Apply station delay to pick
            pick.time = pick.time + dt

    # Attach catalog to tckwargs
    tckwargs.update({'catalog': cat})

    # Construct templates
    tribe = Tribe().construct(**tckwargs)
    # Rename templates & merge multiple traces
    for tmp in tribe:
        tmp.name = ''.join(str(tmp.event.preferred_origin).split('/')[-2:])
        if len(tmp.st) > 2:
            breakpoint()
        tmp.st.merge(method=1)
    # Convert tribe into clusteringtribe
    ctr = ClusteringTribe(templates=tribe.templates)
    # Run template xcorr clustering
    ctr.cluster(**xcckwargs)
    # Populate clusters dataframe
    ctr.populate_event_metadata()

    # Save whole clustering tribe
    ctr.write(str(OUTPUT_DIR/'aqms_event_templates.tgz'))

    # Create subset of highest SNR templates per xcorr cluster
    pref_names = []
    for _gn in ctr._c.xcc.unique():
        # Subset by group
        _df = ctr._c[ctr._c.xcc == _gn]
        # Get highest snr
        snr_max = _df.mean_snr_dB.max()
        # If the max snr meets minimum requirements
        if snr_max >= TEMPLATE_SNR_MIN:
            name = _df[_df.mean_snr_dB == _df.mean_snr_dB.max()].index.values[0]
            pref_names.append(name)
        else:
            continue

    tctr = ctr.get_subset(names=pref_names)
    tctr.write(str(OUTPUT_DIR/'templates_for_match_filter.tgz'))







