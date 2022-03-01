v0.1
------

The Original

v0.4 (0.74?)
------
First distributed

- 0.4.2: Add "--network" argument to lc2SDS_weak
- 0.4.3: Minor code tightening, better test case for read/write
- 0.4.4: Reorganize lc2SDS, add sdpchain compatibility
- 0.4.5: Fixed an error in SDS channel directory names
- 0.4.6: Puts correct band codes for given sampling rate
- 0.4.7: Added leapsecond handling to lc2SDS
- 0.4.8: Fixed lcread to only return data within the requested time bounds
         (inclusive start, exclusive end).  Also stops lc2SDS daily files
         from overlapping
- 0.4.9: lcread/plot no longer quits if there is a bad input file or
         date range is outside the range of some of the files

v1.0
------
b0: Combined lcheapo and lcheapo_obspy
.0: Fixed some bugs associated with creation of the ProcessSteps class
1.0.1: Fixed bugs where lcheapo_obspy was still called, harmonized start_time
       in input file
1.0.2: Fixed a bug in accessing sdpchain library (affected lcinfo lcfix...)
1.0.3: Patch because 1.0.2 on PyPI was bad
1.0.4: Further cleaned up references to sdpchain, changed setup.py to
       automatically find modules (manual system used before didn't put
       sdpchain/ on PyPI)
1.0.5: 
- Changed required python version to 3.8 (some files use '=' specifier
  in f-strings).
- Integrated input_filename wildcard expansion into ProcessSteps.setup_path().
- lcinfo now tries to work even if there is no header
- Add example plots to lctest description in README.md
- Fix bug in lc_test particle_motion part