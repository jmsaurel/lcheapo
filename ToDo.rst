TO DO
======================

- lcplot
    - Allow multiple files to be specified, each with a station code
    - Specify which channels to plot for each file (select parameters)
- lctest
    - Schema and validator for YAML file
    - Add SPOBS1 and HYDROCT StationXML files
- lc2ms_w
    - Add apply clock corrections (on daily basis) if network file is provided
        * Set data quality to "M" (so that it's not "Q")? or some non-sensical
          value?
    - Always sets network code = "XX"
- spectral
    - Allow overlap plot
    - Add codes for noise cleaning
Use `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ to modify this file.
