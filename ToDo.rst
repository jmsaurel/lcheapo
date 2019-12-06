TO DO
======================

- Allow multiple files to be specified, each with a station code
- Specify which channels to plot for each file
- Add lctest, code to plot laboratory test results
    - With YAML file as input (validates the file)
        - List of events, and whether each event is for one or all instruments
        - List of instruments: filename, station_name, inst_type
    - Plots
      * overall time series for each station
      * "quiet" time series for each station
      * "event" (taps, jumps, etc) time series for each station
      * Overlapped spectra for each station
      * Z spectra for all stations
    - Options
      * Maximum length of data to read/plot [1 day] 
      * simply validate YAML file
      * Channels to use for multi-station spectra ['Z,3']
- Add lc2ms, code to convert to miniSEED
    - Can apply clock corrections (on daily basis) if network file is provided
        * Set data quality to "M" (so that it's not "Q")? or some non-sensical
          value?
    - Always sets network code = "XX"
Use `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ to modify this file.
