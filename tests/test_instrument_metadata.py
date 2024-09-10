#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
# from __future__ import (absolute_import, division, print_function,
#                         unicode_literals)
# from future.builtins import *  # NOQA @UnusedWildImport

import os
from pathlib import Path
import unittest
import filecmp
import inspect
import difflib
# import json
import glob
import subprocess

from lcheapo.instrument_metadata import load_station
from obspy.core import UTCDateTime
from obspy.core.inventory import Station


class TestAllMethods(unittest.TestCase):
    """
    Test suite for lcheapo_obspy.
    """
    def setUp(self):
        self.path = Path(__file__).parent.resolve()

    def test_load_station(self):
        """
        test load_station

        Read in an lcheapo file, write it to miniSEED and compare to existing
        miniSEED file
        """
        starttime, endtime = UTCDateTime(2024, 1, 1), UTCDateTime(2024, 1, 1)
        for sample_rate in [62.5, 125., 250., 500.]:
            # for obs_type in ['BBOBS1', 'SPOBS1', 'SPOBS2', 'HYDROCT']:
            for obs_type in ['BBOBS1', 'SPOBS1', 'SPOBS2', 'HYDROCT1']:
                sta = load_station(obs_type, sample_rate, starttime=starttime,
                                   endtime=endtime)
                assert isinstance(sta, Station)
                if obs_type=='SPOBS1':
                    assert len(sta.channels) == 2
                else:
                    assert len(sta.channels) == 4


def suite():
    return unittest.makeSuite(TestAllMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
