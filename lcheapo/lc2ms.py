#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read LCHEAPO data into an obspy stream
"""
# from __future__ import (absolute_import, division, print_function,
#                         unicode_literals)
# from future.builtins import *  # NOQA @UnusedWildImport

import argparse
# import os
import sys
import datetime
import inspect
from pathlib import Path

import lcheapo.sdpchain as sdpchain

from .chan_maps import chan_maps
from .lcread import read as lcread
from .version import __version__


def lc2ms():
    """
    Convert fixed LCHEAPO data to basic miniSEED files

    NO drift or leapsecond correction:
    """
    print(lc2ms.__doc__)
    parser = argparse.ArgumentParser(
        description=inspect.cleandoc(lc2ms.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("infile", help="Input filename")
    parser.add_argument("-t", "--obs_type", default='SPOBS2',
                        help="obs type.  Controls channel and location codes",
                        choices=[s for s in chan_maps])
    parser.add_argument("--station", default='SSSSS',
                        help="station code for this instrument")
    parser.add_argument("--network", default='XX',
                        help="network code for this instrument")
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
    parser.add_argument("--version", action='store_true',
                        help="Print version number and quit")
    args = parser.parse_args()
    parameters = vars(args).copy()
    if args.version is True:
        print(f"Version {__version__}")
        sys.exit(0)

    # ADJUST INPUT PARAMETERS
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
    stream = lcread(Path(args.in_dir) / args.infile, network=args.network,
                    station=args.station, obs_type=args.obs_type)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for tr in stream:
        s = tr.stats
        fname = str(out_dir / '{}.{}.{}.{}.mseed'.format(
                    s.network, s.station, s.location, s.channel))
        tr.write(fname, format='MSEED', encoding='STEIM1', reclen=4096)
    return_code = 0
    sdpchain.make_process_steps_file(
        args.in_dir,
        args.out_dir,
        'lc2ms_weak',
        'create miniSEED file(s) from LCHEAPO file(s)',
        __version__,
        " ".join(sys.argv),
        startTimeStr,
        return_code,
        exec_parameters=parameters)
    sys.exit(return_code)


# ---------------------------------------------------------------------------
# Run 'main' if the script is not imported as a module
# ---------------------------------------------------------------------------
# if __name__ == '__main__':
#     main()
