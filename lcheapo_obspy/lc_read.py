#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read LCHEAPO data into an obspy stream
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

from lcheapo.lcheapo import (LCDataBlock, LCDiskHeader)
from obspy.core import UTCDateTime, Stream, Trace
import argparse
import warnings
import struct
import time
import numpy as np

chan_maps = {'SPOBS1': ['SH3', 'BDH'],
             'SPOBS2': ['BDH', 'SH2', 'SH1', 'SH3'],
             'BBOBS': ['SH2', 'SH1', 'SHZ', 'BDH'],
             'HYDROCT': ['BDH:00', 'BDH:01', 'BDH:02', 'BDH:03']}


def read(filename, starttime=None, endtime=None, network='XX', station='SSSSS',
         chan_map=None, debug=False):
    """
    Read LCHEAPO data into an obspy stream

    :param filename: LCHEAPO filename
    :type  filename: str
    :param starttime: Start time as a ISO8601 string, a UTCDateTime object,
        or a number.  In the latter case, is seconds since the file start
    :type  endtime: :class:`~obspy.core.utcdatetime.UTCDateTime`
    :param endtime: End time as a ISO8601 string, a UTCDateTime object,
        or a number.  In the latter case, is seconds after starttime
    :type  network: str
    :param network: Network code (up to two characters)
    :type  station: str
    :param station: FDSN station name (up to five characters)
    :type  chan_map: list of str, or str
    :param chan_map: Channel mappings: list of names or str of known type
    :return: stream
    :rtype:  class `obspy.stream.Stream`

    .. rubric:: Example

    >>> from lcheapo_obspy import read
    >>> st = read("/path/to/four_channels.lch")
    >>> print(st)  # doctest: +ELLIPSIS
    2 Trace(s) in Stream:
    BW.UH3..EHE | 2010-06-20T00:00:00.279999Z - ... | 200.0 Hz, 386 samples
    BW.UH3..EHZ | 2010-06-20T00:00:00.279999Z - ... | 200.0 Hz, 386 samples
    """
    # Check/complete input variables
    if len(network) > 2:
        network = network[:2]
    if len(station) > 5:
        network = network[:5]
    chan_map = _fill_chan_map(chan_map)

    if debug:
        print('Reading LCHEAPO file {}'.format(filename), end='', flush=True)
        start_time = time.time()
    with open(filename, 'rb') as fp:
        data = _read_data(starttime, endtime, fp)
    if debug:
        print('took {:.2f} seconds'.format(time.time() - start_time))
    data = _stuff_info(data, network, station, chan_map)
    return data


def get_data_timelimits(lcheapo_object):
    """
    Return data start and endtimes

    :param lcheapo_object: Filename or open file like object that contains the
        binary Mini-SEED data.  Any object that provides a read() method will
        be
        considered to be a file like object.
    """

    if hasattr(lcheapo_object, "read"):
        fp = lcheapo_object
    else:
        fp = open(lcheapo_object, 'rb')

    lcHeader = LCDiskHeader()
    block = LCDataBlock()
    lcHeader.readHeader(fp)

    # read starttime
    startBlock = lcHeader.dataStart
    block.seekBlock(fp, startBlock)
    block.readBlock(fp)
    starttime = block.getDateTime()

    # read endtime
    lastBlock = block.determineLastBlock(fp)
    block.seekBlock(fp, lastBlock)
    block.readBlock(fp)
    endtime = block.getDateTime()

    return UTCDateTime(starttime), UTCDateTime(endtime)


def _read_data_old(starttime, endtime, fp):
    """
    Return data

    Reads block by block

    :param starttime: start time
    :type  starttime: :class: `~obspy.UTCDateTime`
    :param endtime: end time
    :type  endtime: :class: `~obspy.UTCDateTime`
    :param fp: file pointer
    :type  fp: class `file`
    """
    starttime, endtime = _get_time_limits(starttime, endtime, fp)

    lcHeader = LCDiskHeader()
    block = LCDataBlock()
    lcHeader.readHeader(fp)
    samp_rate = lcHeader.realSampleRate
    n_start_block = _get_block_number(starttime, fp)
    n_end_block = _get_block_number(endtime, fp)
    n_chans = lcHeader.numberOfChannels
    stream = Stream()

    block.seekBlock(fp, n_start_block)
    arrays = []
    stats = []
    # Stuff first record for each channel
    for i in range(0, n_chans):
        block.readBlock(fp)
        arrays.append(np.array(block.convertDataTo24BitValues()))
        stats.append({'sampling_rate': samp_rate,
                      'starttime': UTCDateTime(block.getDateTime())})
    # Read the rest
    for n_block in range(n_start_block+n_chans, n_end_block + n_chans):
        block.readBlock(fp)
        i_channel = block.muxChannel
        arrays[i_channel] = np.append(arrays[i_channel],
            np.array(block.convertDataTo24BitValues()))
    i = 0
    for array in arrays:
        stream.append(Trace(data=array.astype('int32'), header=stats[i]))
        i += 1
    return stream.slice(starttime, endtime)


def _read_data(starttime, endtime, fp, debug=False):
    """
    Return data

    Reads all blocks at once and extracts as slices

    :param starttime: start time
    :type  starttime: :class: `~obspy.UTCDateTime`
    :param endtime: end time
    :type  endtime: :class: `~obspy.UTCDateTime`
    :param fp: file pointer
    :type  fp: class `file`
    """
    if debug:
        start_read = time.time()
    starttime, endtime = _get_time_limits(starttime, endtime, fp)

    lcHeader = LCDiskHeader()
    block = LCDataBlock()
    lcHeader.readHeader(fp)
    samp_rate = lcHeader.realSampleRate
    n_start_block = _get_block_number(starttime, fp)
    n_end_block = _get_block_number(endtime, fp)
    n_chans = lcHeader.numberOfChannels
    stream = Stream()

    chan_blocks = 1 + int(((n_end_block - n_start_block) / n_chans))
    read_blocks = chan_blocks * n_chans
    block.seekBlock(fp, n_start_block)
    buf = fp.read(read_blocks * 512)
    if debug:
        print('{:.3f}s for fp.read()'.format(time.time()-start_read))
        last_read = time.time()
    dt = np.dtype('b')
    all = np.frombuffer(buf, dtype=dt)
    if debug:
        print('{:.3f}s to .frombuffer'.format(time.time() - last_read))
        last_read = time.time()
    a = np.reshape(all, (read_blocks, 512)) # keep data contiguous
    if debug:
        print('{:.3f}s to reshape'.format(time.time() - last_read))
        last_read = time.time()
    # The following seems to take the same time ...
    # a = all.view()
    # a.shape = (read_blocks,512)
    headers = a[:, :14]
    data = a[:, 14:]
    if debug:
        print('{:.3f}s to extract headers and data'.format(time.time() - last_read))
        last_read = time.time()

    # Extract data
    s_per_block = 166 / samp_rate
    stats = {'sampling_rate': samp_rate}
    for i in range(0, n_chans):
        # Get header information and whether data are contiguous
        first_time = _get_header_time(headers[0, :])
        last_time = _get_header_time(headers[-1, :])
        if debug:
            print('{:d}: {:.3f}s to read times'.format(i, time.time() - last_read))
            last_read = time.time()
        offset = (last_time - first_time) - (chan_blocks - 1) * s_per_block
        if offset > 0:
            warnings.warn('Last block {:g}s ({:g} blocks) late!'.format(
                    offset, offset/s_per_block))
        elif offset < 0:
            warnings.warn('Last block {:g}s ({:g} blocks) early!'.format(
                    -offset, -offset/s_per_block))
        stats['starttime'] = first_time
        # Get data
        chan_data = data[i:read_blocks:n_chans, :].flatten()
        if debug:
            print('   {:.3f}s to slice data'.format(time.time() - last_read))
            last_read = time.time()
        # could be quicker using the np.flat() iterator ?
        t32 = np.array(chan_data[0:498*chan_blocks:3]*(1<<16) +\
                       chan_data[1:498*chan_blocks:3].astype('B')*(1<<8) +\
                       chan_data[2:498*chan_blocks:3].astype('B'),
                       dtype='int32')
        # t32 = np.array([x * (1 << 16) + y * (1 << 8) + z
        #                 for x, y, z in 
        #                 [struct.unpack(">bBB", chan_data[x:x + 3])
        #                  for x in range(0, 498 * chan_blocks, 3)]],
        #                 dtype='int32')
        if debug:
            print('   {:.3f}s to convert to 24-bit'.format(time.time() - last_read))
            last_read = time.time()
        stream.append(Trace(data=t32, header=stats))
        if debug:
            print('   {:.3f}s to stream.append'.format(time.time() - last_read))
            last_read = time.time()
    if debug:
        print('{:.3f}s total read time '.format(time.time() - start_read))
    return stream


def _get_header_time(header):
    """
    Return time from an LCHEAPO header

    header = 14-byte LCHEAPO HEADER
    """
    (msec, second, minute, hour, day, month, year) = struct.unpack(
        '>HBBBBBB', header.data[:8])
    return UTCDateTime(year + 2000, month, day, hour, minute, second,
                       msec + 1000)


def _get_time_limits(starttime, endtime, fp):
    """
    Return starttime and endtime as UTCDateTimes

    startime = UTCDateTime: no change
             = string: convert to UTCDateTime
             = numeric: interpret as seconds from file start
    endtime = UTCDateTime: no change
            = string: convert to UTCDateTime
            = numeric: interpret as seconds starttime

    Also makes sure that starttime and endtime are within the data bounds
    """
    data_start, data_end = get_data_timelimits(fp)
    if not starttime:
        starttime = 0
    if endtime == 0:
        endtime = data_end
    if not endtime:
        endtime = 3600
    if isinstance(starttime, (float, int)):
        starttime = data_start + starttime
    else:
        starttime = UTCDateTime(starttime)
    if isinstance(endtime, (float, int)):
        endtime = starttime + endtime
    # Handle bad time ranges
    assert endtime > starttime, "endtime is before starttime"
    assert starttime < data_end, f"starttime is after data end ({data_end})"
    assert endtime > data_start, f"endtime is before data start ({data_start})"
    if starttime < data_start:
        starttime = data_start
    if endtime > data_end:
        endtime = data_end
    return starttime, endtime


def _get_block_number(time, fp):
    """
    Return block number containing first channel containing the given time

    :param time: the time
    :type  time: class `~obspy.UTCDateTime`
    :param fp: file pointer
    :type  fp: class `file`
    """
    lcHeader = LCDiskHeader()
    block = LCDataBlock()
    lcHeader.readHeader(fp)
    data_start_block = lcHeader.dataStart
    block.seekBlock(fp, data_start_block)
    block.readBlock(fp)

    starttime, _ = get_data_timelimits(fp)
    block_len_s = block.numberOfSamples / lcHeader.realSampleRate
    num_chans = lcHeader.numberOfChannels
    record_offset = int((time-starttime) / block_len_s)

    return data_start_block + record_offset * num_chans


def _stuff_info(stream, network, station, chan_map):
    """
    Put network, station and channel information into stream

    :param stream: obspy stream
    :type  stream: class `~obspy.stream.Stream`
    :param network: network code
    :type  network: str
    :param station: station code
    :type  station: str
    :param chan_map: channel mapping
    :type  chan_map: list

    :return data: informed data
    :rtype  data: class `~obspy.stream.Stream`
    """
    chan_map_list = list(chan_map)
    for trace in stream:
        trace.stats.network = network
        trace.stats.station = station
        if len(chan_map) >= 1:
            chan_loc = chan_map_list[0].split(':')
            trace.stats.channel = chan_loc[0]
            if len(chan_loc) > 1:
                trace.stats.location = chan_loc[1]
            chan_map_list.pop(0)
    return stream


def _fill_chan_map(inp):
    """
    Fill and/or verify channel mappings

    inp: either a known LCHEAPO INSTRUMENT MODEL or a list of channel names
    """
    if isinstance(inp, str):
        if inp in chan_maps:
            return chan_maps[inp]
        else:
            print('Unknown LCHEAPO instrument: {}, choose from {}'.format(
                inp, [s for s in chan_maps]))
    elif isinstance(inp, list):
        if _valid_chan_map(inp):
            return(inp)
        else:
            print("Bad channel map")
    else:
        print("Input is neither an inst type nor a list of channel names")
    return None


def _valid_chan_map(chan_map):
    """
    Verify that chan_map is valid
    """
    if not isinstance(chan_map, list):
        return False
    for elem in chan_map:
        if not isinstance(elem, str):
            return False
        if len(elem) <= 3:
            return True
        else:
            if ':' in elem:
                elems = elem.split(':')
                if len(elems) == 2:
                    if len(elems[0]) <= 3 and len(elems[1]) <= 2:
                        return True
    return False


def _plot_command():
    """
    Command-line plotting interface
    """
    parser = argparse.ArgumentParser(
        description=__doc__)
    parser.add_argument("infile", help="Input filename(s)")
    parser.add_argument(
        "starttime", default=None,
        help="start time (ISO8601, or seconds from file start)")
    parser.add_argument("endtime", default=None,
                        help="end time (ISO8601, or seconds from start time)")
    parser.add_argument("-t", "--obs_type", default='SPOBS2', help="obs type",
                        choices=[s for s in chan_maps])
    parser.add_argument("-n", "--network", default='XX', help="network code")
    parser.add_argument("-s", "--station", default='SSSSS',
                        help="station code")
    args = parser.parse_args()

    endtime = _normalize_time_arg(args.endtime)
    if endtime == 0:
        endtime = 3600.
    stream = read(args.infile,
                  _normalize_time_arg(args.starttime),
                  _normalize_time_arg(args.endtime),
                  network=args.network,
                  station=args.station,
                  chan_map=args.obs_type)
    print(stream)
    stream.plot(size=(800, 600), equal_scale=False, method='full')


def _to_mseed_command():
    """
    Command-line conversion to miniSEED
    """
    parser = argparse.ArgumentParser(
        description=__doc__)
    parser.add_argument("infile", help="Input filename(s)")
    parser.add_argument(
        "starttime", default=None,
        help="start time (ISO8601, or seconds from file start)")
    parser.add_argument("endtime", default=None,
                        help="end time (ISO8601, or seconds from starttime)")
    parser.add_argument("-g", "--granularity", type=int, default=86400,
                        help="granularity for reading (seconds)")
    parser.add_argument("-t", "--obs_type", default='SPOBS2', help="obs type",
                        choices=[s for s in chan_maps])
    parser.add_argument("-n", "--network", default='XX', help="network code")
    parser.add_argument("-s", "--station", default='SSSSS',
                        help="station code")
    args = parser.parse_args()

    stream = read(args.infile,
                  _normalize_time_arg(args.starttime),
                  _normalize_time_arg(args.endtime),
                  network=args.network,
                  station=args.station,
                  chan_map=args.obs_type)

    for tr in stream:
        fname = tr.stats.starttime.strftime('%Y-%m-%dT%H%M%S') +\
                "{}.{}.{}.{}.mseed".format(tr.stats.network, tr.stats.station,
                                           tr.stats.channel, tr.stats.location)
        tr.write(fname, format='MSEED', encoding=3, reclen=4096)


def _normalize_time_arg(a):
    """
    Convert time from string to float if it is numeric
    """
    try:
        temp = float(a)
    except ValueError:
        return a
    else:
        return temp


# ---------------------------------------------------------------------------
# Run 'main' if the script is not imported as a module
# ---------------------------------------------------------------------------
# if __name__ == '__main__':
#     main()
