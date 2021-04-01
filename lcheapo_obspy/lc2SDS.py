#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read LCHEAPO data into an obspy stream
"""
# from __future__ import (absolute_import, division, print_function,
#                         unicode_literals)
# from future.builtins import *  # NOQA @UnusedWildImport

import argparse
import warnings
# import os
import sys
import datetime
import inspect
from pathlib import Path

from obspy.core import UTCDateTime
import lcheapo.sdpchain as sdpchain
from progress.bar import IncrementalBar

from .chan_maps import chan_maps
from .lcread import read as lcread, get_data_timelimits
from .version import __version__

# chan_maps = {'SPOBS1': ['SH3:00', 'BDH:00'],
#              'SPOBS2': ['BDH:00', 'SH2:00', 'SH1:00', 'SH3:00'],
#              'BBOBS1': ['BH2:00', 'BH1:00', 'BHZ:00', 'BDH:00'],
#              'HYDROCT': ['BDH:00', 'BDH:01', 'BDH:02', 'BDH:03']}


def lc2SDS():
    """
    Convert fixed LCHEAPO data to SeisComp Data Structure

    SIMPLE drift correction: only at beginning of each daily file, no
    offset information put in header, no modification of data quality field.
    Writes to directory SDS in the output directory.
    """
    print(lc2SDS.__doc__)
    parser = argparse.ArgumentParser(
        description=inspect.cleandoc(lc2SDS.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("infiles", nargs='+',
                        help="Input filename(s).  If there are captured "
                             "wildcards (put in '' so that they aren't "
                             "interpreted by the shell), will expand them "
                             "in the input directory")
    parser.add_argument("-t", "--obs_type", default='SPOBS2',
                        help="obs type.  Controls channel and location codes",
                        choices=[s for s in chan_maps])
    parser.add_argument("--station", default='SSSSS',
                        help="station code for this instrument")
    parser.add_argument("--network", default='XX',
                        help="network code for this instrument")
    parser.add_argument("-s", "--start_times", nargs='+', type=UTCDateTime,
                        metavar=("REF_START", "INST_START"),
                        help="Start datetimes for the reference (usually GPS) "
                             "and instrument.  If only one value is provided, "
                             "it will be used for both")
    parser.add_argument("-e", "--end_times", nargs=2, type=UTCDateTime,
                        metavar=("REF_END", "INST_END"),
                        help="End datetimes for the reference and instrument")
    parser.add_argument("-d", dest="base_dir", metavar="BASE_DIR",
                        default='.', help="base directory for files")
    parser.add_argument("-i", dest="in_dir", metavar="IN_DIR", default='.',
                        help="input file directory (absolute, " +
                             "or relative to base_dir)")
    parser.add_argument("-o", dest="out_dir", metavar="OUT_DIR", default='.',
                        help="output file directory (absolute, " +
                             "or relative to base_dir)")
    parser.add_argument("-v", "--verbose", action='store_true',
                        help="verbose output")
    args = parser.parse_args()
    args.in_dir, args.out_dir = sdpchain.setup_paths(args.base_dir,
                                                     args.in_dir,
                                                     args.out_dir)
    # Expand captured wildcards
    print(f'{args.infiles=}')
    args.infiles = [x.name for f in args.infiles
                    for x in Path(args.in_dir).glob(f)]
    print(f'expanded {args.infiles=}')

    startTimeStr = datetime.datetime.strftime(datetime.datetime.utcnow(),
                                              '%Y-%m-%dT%H:%M:%S')
    for infile in args.infiles:
        lc_start, lc_end = get_data_timelimits(Path(args.in_dir) / infile)

        if args.start_times and args.end_times:
            ref_start = args.start_times[0]
            if len(args.start_times) > 1:
                inst_start = args.start_times[1]
            else:
                inst_start = ref_start
            ref_end, inst_end = args.end_times
            if inst_start == 0:
                inst_start = ref_start
            inst_start_offset = inst_start - ref_start
            inst_drift = ((inst_end - ref_end) - inst_start_offset)\
                / (ref_end - inst_start)
            print('instrument start offset = {:g}s, drift rate = {:.4g}'
                  .format(inst_start_offset, inst_drift))
            # quality_flag = 'Q'  # Don't know how to put this in miniSEED
        else:
            ref_start, inst_start = lc_start, lc_start
            inst_start_offset = 0
            inst_drift = 0
            warnings.warn('Could not calculate clock drift, assuming zero!')
            # quality_flag = 'D'  # Don't know how to put this in miniSEED

        lc_start_day = lc_start.replace(hour=0, minute=0, second=0,
                                        microsecond=0)
        lc_end_day = lc_end.replace(hour=0, minute=0, second=0, microsecond=0)
        stime = lc_start_day
        bar = IncrementalBar(f'Processing {infile}', max=(lc_end_day-lc_start_day)/86400)
        while stime <= lc_end_day:
            inst_offset = inst_start_offset + inst_drift * (stime - ref_start)
            starttime = stime + inst_offset
            endtime = starttime + 86400
            if args.verbose:
                print('{}, inst_offset = {:.3f}s: reading {}-{}'.format(
                    stime.strftime('%Y-%m-%d'), inst_offset,
                    starttime.isoformat(), endtime.isoformat()))
            stream = lcread(Path(args.in_dir) / infile,
                            starttime=starttime,
                            endtime=endtime,
                            network=args.network,
                            station=args.station,
                            obs_type=args.obs_type)

            for tr in stream:
                s = tr.stats
                # Correct drift
                s.starttime -= inst_offset
                # s.mseed['dataquality'] = quality_flag

                # Write file
                dirname = Path(args.out_dir) / 'SDS' / str(stime.year) /\
                    s.network / s.station / f'{s.channel}.D'
                fname = '{}.{}.{}.{}.D.{}.{:03d}'.format(
                    s.network, s.station, s.location, s.channel,
                    stime.year, stime.julday)
                dirname.mkdir(parents=True, exist_ok=True)
                tr.write(str(dirname / fname), format='MSEED',
                         encoding='STEIM1', reclen=4096)
            bar.next()
            stime += 86400
        bar.finish()

    return_code = 0
    parameters = vars(args)
    parameters['start_times'] = [str(x) for x in parameters['start_times']]
    parameters['end_times'] = [str(x) for x in parameters['end_times']]
    sdpchain.make_process_steps_file(
        args.in_dir,
        args.out_dir,
        'lc2SDS_weak',
        'Create or add to an SDS archive from an lcheapo file',
        __version__,
        " ".join(sys.argv),
        startTimeStr,
        return_code,
        exec_parameters=vars(args))
    sys.exit(return_code)


# def _to_mseed_command():
#     """
#     Command-line conversion to miniSEED
#     """
#     parser = argparse.ArgumentParser(
#         description=__doc__)
#     parser.add_argument("infile", help="Input filename(s)")
#     parser.add_argument(
#         "starttime", default=0,
#         help="start time (ISO8601, or seconds from file start)")
#     parser.add_argument("endtime", default=None,
#                         help="end time (ISO8601, or secs from starttime)")
#     parser.add_argument("-g", "--granularity", type=int, default=86400,
#                         help="granularity for reading (seconds)")
#     parser.add_argument("-t", "--obs_type", default='SPOBS2',
#                         help="obs type",
#                         choices=[s for s in chan_maps])
#     parser.add_argument("-n", "--network", default='XX', help="network code")
#     parser.add_argument("--station", default='SSSSS',
#                         help="station code")
#     args = parser.parse_args()
#
#     stream = read(args.infile,
#                   _normalize_time_arg(args.starttime),
#                   _normalize_time_arg(args.endtime),
#                   network=args.network,
#                   station=args.station,
#                   obs_type=args.obs_type)
#
#     for tr in stream:
#         fname = tr.stats.starttime.strftime('%Y-%m-%dT%H%M%S') +\
#                 "{}.{}.{}.{}.mseed".format(tr.stats.network,
#                                            tr.stats.station,
#                                            tr.stats.channel,
#                                            tr.stats.location)
#         tr.write(fname, format='MSEED', encoding=3, reclen=4096)

# ---------------------------------------------------------------------------
# Run 'main' if the script is not imported as a module
# ---------------------------------------------------------------------------
# if __name__ == '__main__':
#     main()
