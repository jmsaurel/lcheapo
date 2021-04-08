#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
# from __future__ import (absolute_import, division, print_function,
#                         unicode_literals)
# from future.builtins import *  # NOQA @UnusedWildImport

import os
import unittest
import filecmp
import inspect
import difflib
# import json
import glob

from lcheapo_obspy.lcread import read as lcread, band_code_sps
from lcheapo_obspy.yaml_json import validate
from lcheapo_obspy.lc2SDS import _adjust_leapseconds, _leap_correct
from obspy.core import UTCDateTime


class TestAllMethods(unittest.TestCase):
    """
    Test suite for lcheapo_obspy.
    """
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(inspect.getfile(
            inspect.currentframe())))
        self.testing_path = os.path.join(self.path, "data")
        self.exec_path = os.path.split(self.path)[0]
        self.examples_path = os.path.join(self.exec_path, '_examples')

    def assertTextFilesEqual(self, first, second, msg=None):
        first_f = open(first)
        first_str = first_f.read()
        second_f = open(second)
        second_str = second_f.read()
        first_f.close()
        second_f.close()

        if first_str != second_str:
            first_lines = first_str.splitlines(True)
            second_lines = second_str.splitlines(True)
            delta = difflib.unified_diff(
                first_lines, second_lines,
                fromfile=first, tofile=second)
            message = ''.join(delta)

            if msg:
                message += " : " + msg

            self.fail("Multi-line strings are unequal:\n" + message)

    def assertBinFilesEqual(self, first, second, msg=None):
        """ Compares two binary files """
        self.assertTrue(filecmp.cmp(first, second))

    def test_read(self):
        """
        test lcread

        Read in an lcheapo file, write it to miniSEED and compare to existing
        miniSEED file
        """
        test_fname = 'XX.TEST.2019-11-07.mseed'
        infile = os.path.join(self.examples_path,
                              '20191107T14_SPOBS09_F02.raw.lch')
        compare_file = os.path.join(self.testing_path, test_fname)
        stream = lcread(infile, station='TEST', network='XX',
                        obs_type='SPOBS2')
        stream.write(test_fname, 'MSEED', encoding='STEIM1', byteorder='<')
        self.assertBinFilesEqual(test_fname, compare_file)
        os.remove(test_fname)

    def test_lctest_validate(self):
        """validate lctest YAML files in _examples directory"""
        for f in glob.glob(os.path.join(self.examples_path, '*.yaml')):
            validate(f, quiet=True)

    def test_band_code_sps(self):
        """validate lcread band_code_sps() function"""
        self.assertEqual("F", band_code_sps("B", 1000))
        self.assertEqual("G", band_code_sps("S", 1000))
        self.assertEqual("C", band_code_sps("B", 500))
        self.assertEqual("D", band_code_sps("S", 500))
        self.assertEqual("C", band_code_sps("B", 250))
        self.assertEqual("D", band_code_sps("S", 250))
        self.assertEqual("H", band_code_sps("B", 200))
        self.assertEqual("E", band_code_sps("S", 200))
        self.assertEqual("H", band_code_sps("B", 100))
        self.assertEqual("E", band_code_sps("S", 100))
        self.assertEqual("B", band_code_sps("B", 50))
        self.assertEqual("S", band_code_sps("S", 50))
        self.assertEqual("B", band_code_sps("B", 10))
        self.assertEqual("S", band_code_sps("S", 10))
        self.assertEqual("L", band_code_sps("B", 1))
        # Try alternative input band code
        self.assertEqual("B", band_code_sps("C", 50))
        self.assertEqual("S", band_code_sps("D", 50))
        # Test invalid input band code
        self.assertRaises(NameError, band_code_sps, *["N", 100])
        # Test short period band code for < 10 sps
        self.assertRaises(NameError, band_code_sps, *["S", 9])

    def test_leap_seconds(self):
        """ test leapsecond routines in lc2SDS """
        lstimes = ["2020-12-31:23:59:60", "2021-06-30T23:59:58"]
        self.assertRaises(IndexError, _adjust_leapseconds, *[lstimes, '+++'])
        self.assertRaises(ValueError, _adjust_leapseconds, *[lstimes, '+p'])
        lstm, lstp = _adjust_leapseconds(lstimes, '+')
        self.assertEqual(lstp, '++')
        self.assertIsInstance(lstm[0], UTCDateTime)
        lstm, lstp = _adjust_leapseconds(lstimes, '+-')
        self.assertEqual(_leap_correct(UTCDateTime('2020-12-31'), lstm, lstp), 0)
        self.assertEqual(_leap_correct(UTCDateTime('2021-01-01'), lstm, lstp), -1)
        self.assertEqual(_leap_correct(UTCDateTime('2021-04-01'), lstm, lstp), -1)
        self.assertEqual(_leap_correct(UTCDateTime('2021-07-01'), lstm, lstp), 0)
        
    # TESTS TO ADD:
    # def test_lcplot(self):
    #   """
    #   Verify that the plot is what we expect
    #   """

    # def test_lc2ms(self):
    #   """
    #   Verify that the output miniSEED file is as expected
    #   """

    # def test_lc2SDS(self):
    #   """
    #   Verify correct drift correction and SDS output
    #   """

    # def test_lctest(self):
    #   """
    #   Verify that the output of an lctest run is what we expect
    #   """


def suite():
    return unittest.makeSuite(TestAllMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
