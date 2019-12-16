TO DO
======================

- lcplot
    - Allow multiple files to be specified, each with a station code
    - Specify which channels to plot for each file (select parameters)
- lctest
    - Schema and validator for YAML file
    - Add integration of data/LCHEAPO.station.xml transfer functions
- lcfix
    - Allow header file
    - Make main code more modular (more functions/subroutines) for clarity
- spectral
    - Allow overlap plot
    - Add codes for noise cleaning
- Add:
    - lc2ms, code to convert to miniSEED
        - Can apply clock corrections (on daily basis) if network file is provided
            * Set data quality to "M" (so that it's not "Q")? or some non-sensical
              value?
        - Always sets network code = "XX"
    - lccut: dd-like command to cut lcheapo files
Use `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ to modify this file.
