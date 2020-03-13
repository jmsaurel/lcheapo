#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read LCHEAPO data into an obspy stream
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

import argparse
import warnings
import struct
import time
import os
import sys
import inspect

import numpy as np
from obspy.core import UTCDateTime, Stream, Trace
from obspy import read_inventory

from lcheapo.lcheapo import (LCDataBlock, LCDiskHeader)

chan_maps = {'SPOBS1': ['SH3', 'BDH'],
             'SPOBS2': ['BDH', 'SH2', 'SH1', 'SH3'],
             'BBOBS1': ['BH2', 'BH1', 'BHZ', 'BDH'],
             'HYDROCT': ['BDH:00', 'BDH:01', 'BDH:02', 'BDH:03']}


def read(filename, starttime=None, endtime=None, network='XX', station='SSSSS',
         obs_type=None):
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
    :type  obs_type: str
    :param obs_type: OBS type (must match a key in chan_maps)
    :return: stream
    :rtype:  :class:`~obspy.core.stream.Stream`

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
    if not obs_type:
        warnings.warn('No obs_type provided, assuming SPOBS2')
        obs_type = 'SPOBS2'

    with open(filename, 'rb') as fp:
        data = _read_data(starttime, endtime, fp)
    data = _stuff_info(data, network, station, obs_type)
    return data


def get_data_timelimits(lcheapo_object):
    """
    Return data start and endtimes

    :param lcheapo_object: Filename or open file-like object that contains the
        binary Mini-SEED data.  Any object that provides a read() method will
        be considered a file-like object.
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


def _read_data(starttime, endtime, fp):
    """
    Return data

    Reads all blocks at once and extracts as slices

    :param starttime: start time
    :type  starttime: :class:`~obspy.UTCDateTime`
    :param endtime: end time
    :type  endtime: :class:`~obspy.UTCDateTime`
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

    chan_blocks = 1 + int(((n_end_block - n_start_block) / n_chans))
    read_blocks = chan_blocks * n_chans
    block.seekBlock(fp, n_start_block)
    buf = fp.read(read_blocks * 512)
    dt = np.dtype('b')
    all = np.frombuffer(buf, dtype=dt)
    a = np.reshape(all, (read_blocks, 512))  # keep data contiguous
    # The following seems to take the same time ...
    # a = all.view()
    # a.shape = (read_blocks,512)
    headers = a[:, :14]
    data = a[:, 14:]

    # Extract data
    s_per_block = 166 / samp_rate
    stats = {'sampling_rate': samp_rate}
    for i in range(0, n_chans):
        # Get header information and whether data are contiguous
        first_time = _get_header_time(headers[0, :])
        last_time = _get_header_time(headers[-1, :])
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
        # could be quicker using the np.flat() iterator ?
        t32 = np.array(chan_data[0: 498*chan_blocks: 3]*(1 << 16)
                       + chan_data[1: 498*chan_blocks: 3].astype('B')
                       * (1 << 8)
                       + chan_data[2: 498*chan_blocks: 3].astype('B'),
                       dtype='int32')
        stream.append(Trace(data=t32, header=stats))
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
    :type  time: :class:`~obspy.UTCDateTime`
    :param fp: file pointer
    :type  fp: :class:`~file`
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


def _stuff_info(stream, network, station, obs_type):
    """
    Put network, station and channel information into station stream

    :param stream: obspy stream
    :type  stream: :class:`~obspy.stream.Stream`
    :param network: network code
    :type  network: str
    :param station: station code
    :type  station: str
    :param obs_type: type of obs (must be a key in channel_maps)
    :type  obs_type: str

    :return data: informed data
    :rtype  data: :class:`~obspy.stream.Stream`
    """
    assert obs_type in chan_maps
    chan_map = chan_maps[obs_type]
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
        trace.stats.response = _load_response(obs_type, trace.stats.channel,
                                              trace.stats.starttime)
    return stream

def _load_response(obs_type, channel, datetime):
    """
    Load response corresponding to OBS type and component

    :param obs_type: obs type
    :type  obs_type: str
    :param channel: trace channel code
    :param datetime: time for which to get response
    :type datetime: :class:`~obspy.UTCDateTime`

    :return data: instrument response
    :rtype  data: :class:`~obspy.core.response.Response`
    """
    assert obs_type in chan_maps
    try:
        basepath = os.path.dirname(os.path.abspath(inspect.getfile(
                                       inspect.currentframe())))
        inv_file = os.path.join(basepath, f'data/{obs_type}.station.xml')
        inv = read_inventory(inv_file)
    except:
        print(f'Could not read inventory file {inv_file}')
        sys.exit()
    return inv.get_response(f'XX.SPOB2.00.{channel}',datetime)


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
                  obs_type=args.obs_type)
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


# def _read_data_old(starttime, endtime, fp):
#     """
#     Return data
#
#     Reads block by block
#
#     :param starttime: start time
#     :type  starttime: :class: `~obspy.UTCDateTime`
#     :param endtime: end time
#     :type  endtime: :class: `~obspy.UTCDateTime`
#     :param fp: file pointer
#     :type  fp: class `file`
#     """
#     starttime, endtime = _get_time_limits(starttime, endtime, fp)
#     lcHeader = LCDiskHeader()
#     block = LCDataBlock()
#     lcHeader.readHeader(fp)
#     samp_rate = lcHeader.realSampleRate
#     n_start_block = _get_block_number(starttime, fp)
#     n_end_block = _get_block_number(endtime, fp)
#     n_chans = lcHeader.numberOfChannels
#     stream = Stream()
#
#     block.seekBlock(fp, n_start_block)
#     arrays = []
#     stats = []
#     # Stuff first record for each channel
#     for i in range(0, n_chans):
#         block.readBlock(fp)
#         arrays.append(np.array(block.convertDataTo24BitValues()))
#         stats.append({'sampling_rate': samp_rate,
#                       'starttime': UTCDateTime(block.getDateTime())})
#     # Read the rest
#     for n_block in range(n_start_block+n_chans, n_end_block + n_chans):
#         block.readBlock(fp)
#         i_channel = block.muxChannel
#         arrays[i_channel] = np.append(arrays[i_channel],
#             np.array(block.convertDataTo24BitValues()))
#     i = 0
#     for array in arrays:
#         stream.append(Trace(data=array.astype('int32'), header=stats[i]))
#         i += 1
#     return stream.slice(starttime, endtime)

# ---------------------------------------------------------------------------
# Run 'main' if the script is not imported as a module
# ---------------------------------------------------------------------------
# if __name__ == '__main__':
#     main()
