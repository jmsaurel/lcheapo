import warnings
from pathlib import Path
from obspy.core.inventory import read_inventory, Inventory, Station
from copy import deepcopy

# SPOBS1 channel map is flipped from lcplotw, to make SISMANTILLES2 data work
# Although it does make sense that the SPOBS1 and SPOBS2 would have the same
# hydrophone address.
chan_maps = {'SPOBS1': ['BDH:00', 'SH3:00'],
             'SPOBS2': ['BDH:00', 'SH2:00', 'SH1:00', 'SH3:00'],
             'BBOBS1': ['BH2:00', 'BH1:00', 'BHZ:00', 'BDG:00'],
             'HYDROCT1': ['BDH:01', 'BDH:02', 'BDH:03', 'BDH:04']}

obs_types = list(chan_maps.keys())


def load_station(obs_type, sample_rate, **kwargs):
    """
    Load station and channels corresponding to OBS type

    Args:
        obs_type (str): obs type (must be in channel_maps)
        sample_rate (float): sampling rate
        kwargs (dict): any argument for Inventory.select (channel, location,
            starttime, endtime, ...)

    Returns:
        resp (:class:`~obspy.core.station.Station`): station info
    """
    if not obs_type in obs_types:
        raise Exception(f'{obs_type=} not in {obs_types=}')
    try:
        inv_file = Path(__file__).parent / "data" / f'{obs_type}.INSU-IPGP.station.xml'
        inv = read_inventory(inv_file)
    except Exception:
        warnings.warn(f'Could not read inventory file {inv_file}')
        return []
    assert isinstance(inv, Inventory), f"{obs_type=}: {inv=} is not an Inventory"
    assert len(inv.networks) == 1, 'Error: inventory has more than one network'
    inv_select = inv.select(**kwargs)
    if len(inv_select) == 0:
        raise Exception(f'{obs_type} has no channels matching {kwargs=}')
    station = None
    for s in inv_select[0].stations:
        if s.channels[0].sample_rate == sample_rate:
            if station is None:
                station = s
            else:
                raise Exception(f'{obs_type}: {sample_rate=} and {kwargs=} matches more than one instance')
    if station is None:
        raise Exception(f'{obs_type}: {sample_rate=} and {kwargs=} matched no instances')
    assert isinstance(station, Station), f"{obs_type=}: {sta=} is not a Station"
    return station
