#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to test the lcheapo functions
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

import os
import unittest
import filecmp
import inspect
import difflib
import json
import glob

from lcheapo_obspy.lcread import read as lcread
from lcheapo_obspy.yaml_json import validate


class TestLCHEAPOMethods(unittest.TestCase):
    """
    Test suite for nordic io operations.
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
        
        Currently tests only for failure, should verify that what is read
        corresponds to an expected waveform stream
        """
        lcread(os.path.join(self.examples_path,
                            '20191107T14_SPOBS09_F02.raw.lch'))

    def test_lctest_validate(self):
        """validate lctest YAML files in _examples directory"""
        for f in glob.glob(os.path.join(self.examples_path, '*.yaml')):
            validate(f, quiet=True)

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
    return unittest.makeSuite(TestLCHEAPOMethods, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
