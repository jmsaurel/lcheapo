# lcheapo

Viewing and modifying LCHEAPO OBS data

## Overview

### Command-line programs

Type ``{command} -h`` to get a list of parameters and options

#### Non-modifying programs
View the files and/or their data

| Program     | description                                           |
| ----------- | ----------------------------------------------------- |
| lcdump      | dump raw information from LCHEAPO files               |
| lcinfo      | return basic information about an LCHEAPO file        |
| lcplot      | plot an LCHEAPO file                                  |
| lctest      | plot LCHEAPO tests                                    |
| lc_examples | create a directory with examples of lcplot and lctest |

#### Modifying programs
Change the files and/or their data.  Use the *SDPCHAIN* protocols for FAIR-
compliant data:

- Create/append to a process-steps.json file
- Read from input directory (-i) and output to (-o)

| Program     | description                                                                   |
| ----------- | ----------------------------------------------------------------------------- |
| lccut       | extract section of an LCHEAPO file                                            |
| lcfix       | fix common bugs in an LCHEAPO file                                            |
| lcheader    | create an LCHEAPO header + directory                                          |
| sdpcat      | concatenate data files                                                        |
| sdpstep     | run a command line tool and save info to process-steps file                   |
| lc2ms_weak  | converts LCHEAPO file to basic miniSEED files                                 |
| lc2SDS_weak | converts LCHEAPO file to SeisComp Data Structure, with basic drift correction |

## lctest control files

lctest uses YAML-format control files to indicate what kind of plots to
output.  There are 4 main sections in each file:

  - **``input``**: input data parameters
      - **``start_time``**: start time to read (0 means from the beginning
        of each file)

      - **``end_time``**: end time to read (0 means to the end of each file)

      - **``datafiles``**: a list of the LCHEAPO files to read

          - **``name``**: the filename
 
          - **``obs_type``**: the obs type (used for spectral instrument 
            responses,  possible values are given in the help for ``lcplot``
 
          - **``station``**: station name to use for this file in the plots

  - **``output``**: output plot file parameters
  
      - **``show``**: show the plots?  If False, just save them to files
      - **``filebase``**: each output file will start with this

  - **``plots``**: the plots to make
      - **``time_series``**: a list of standard waveform plots
      
          - **``description``**: text to put in plot title
 
          - **``select``**: parameters to use to select a subset of all the
            waveforms (see obspy.core.stream.Stream.select())

      - **``start_time``** and **``end_time``**: allow you to plot a 
        subwindow of the read data.  Set = `null` to plot all
        read data
 
      - **``particle_motion``**: list of particle motion plot types

      - **``spectra``**: list of spectra plots
      
          - **``description``**: as in ``time_series``
          - **``select``**: as in ``time_series``
          - **``start_time``** and **``end_time``**: as in ``time_series``

      - **``stack``**: plot waveforms from different times on the same
        channel together.  This useful for single stations where you did
        the same thing (tap, lift, etc) several times, to be sure that the
        input and response were consistent.

          - **``description``**: as in ``time_series``
 
          - **``orientation_codes``**: list of orientation codes to plot
            (one plot will be made for each orientation code (the last
            letter of the channel name))

          - **``times``**: list of times to plot at (each one
            "yyyy-mm-ddTHH:MM:SS")

          - **``offset_before.s``**: start_time will be this many seconds
            before each ``time``

          - **``offset_after.s``**: end_time will be this many seconds after
            each ``time``

      - **``particle_motion``**: list of particle motion plots to make.  Used
        to confirm the polarity/orientation of the channels.  Each plot will
        contain time series subplots for each channel and a particle motion
        plot combining the two.  Generally, the time span of the time-series
        plot should envelope that of the particle motion plot.  A bit silly
        that it uses ``times``, like ``stacks`` does, since it doesn't plot
        all of the times together.  I think I did it because we generally
        look at particle motions for the same taps that we do stacks on.

          - **``description``**: as in ``time_series``
 
          - **``orientation_code_x``**: orientation code to plot on the x axis

          - **``orientation_code_y``**: orientation code to plot on the y axis

          - **``times``**: list of times to plot at (each one
            "yyyy-mm-ddTHH:MM:SS")

          - **``offset_before.s``**: start_time for particle motion plot will
            be this many seconds before each ``time``

          - **``offset_after.s``**: end_time will for particle motion plot be
            this many seconds after each ``time``

          - **``offset_before_ts.s``**: as above, but for the time series plot

          - **``offset_after_ts.s``**: ditto

  - **``plot_globals``**: Default values for each type of plot.  Use the
    same names and values as for **``plots``**

### Example plots

#### Time series
```
        -   description: "Entire time series"
            select: {station: "*"}
            start_time: "2022-02-22T10:00:01"
            end_time: "2022-02-25T15:25:25"
```
![](README_images/BB02-V1_3-tests_Entire_time_series_ts.png)

```
        -   description: "Quiet time"
            select: {station: "*"}
            start_time: "2022-02-23T21:00:00"
            end_time: "2022-02-24T03:00:00"
```
#### spectra
```
plot_globals:
    spectra:
        window_length.s: 1024
plots:
    spectra:
        -   description: "Quiet time"
            select: {station: "*"}
            start_time: "2022-02-23T21:00:00"
            end_time: "2022-02-24T03:00:00"
```
#### stack
```
        -   description: "Stack, Jump South"
            orientation_codes: ["Z"]
            offset_before.s: 0.3
            offset_after.s: 1
            times:
            -    "2022-02-25T13:57:00.66"
            -    "2022-02-25T13:58:00.53"
            -    "2022-02-25T13:59:00.2"
```

#### particle_motion
```
        -   description: "Stack, Jump South"
            orientation_code_x: "2"
            orientation_code_y: "1"
            offset_before.s: 0.00
            offset_after.s: 0.03
            offset_before_ts.s: 0.2
            offset_after_ts.s: 1
            times:
            -    "2022-02-25T13:57:00.66"
            -    "2022-02-25T13:58:00.53"
            -    "2022-02-25T13:59:00.2"
```

### Examples

#### 1: Analysing one station

``` yaml
---
input:
    start_time: 0
    end_time: 0
    datafiles:
        -   name: "Data_BB07_04_10_12.raw.lch"
            obs_type: "BBOBS"
            station: "TEST"
    description: "Tests on BBOBS"
output:
    show: True
    filebase: "BB07tests"
plot_globals:
    spectra:
        window_length.s: 1024
plots:
    time_series:
        -   description: "Entire time series"
            select: {station: "*"}
            start_time: null
            end_time: null
        -   description: "Quiet time"
            select: {station: "*"}
            start_time: "2012-10-05T02:00:00"
            end_time: "2012-10-05T03:05:00"
        -   description: "Fake stack time"
            select: {station: "*"}
            start_time: "2012-10-05T05:12:00"
            end_time: "2012-10-05T05:21:00"
    spectra:
        -   description: "Quiet time"
            select: {station: "*"}
            start_time: "2012-10-05T02:00:00"
            end_time: "2012-10-05T03:05:00"
    stack:
        -   description: "Fake stack, no tests run"
            orientation_codes: ["Z"]
            offset_before.s: 1
            offset_after.s: 5
            times: 
                - "2012-10-05T05:12:10"
                - "2012-10-05T05:13:50"
                - "2012-10-05T05:17:55"
                - "2012-10-05T05:20:25"
    particle_motion:
        -   description: "rubber hammer tap from S* to N*"
            orientation_code_y: "1"
            orientation_code_x: "2"
            times: 
                - "2019-11-07T14:09:16.65"
                - "2019-11-07T14:09:26.55"
                - "2019-11-07T14:09:36.5"
                - "2019-11-07T14:09:46.75"
                - "2019-11-07T14:09:56.75"
        -   description: "rubber hammer tap from W* to E*"
            orientation_code_y: "1"
            orientation_code_x: "2"
            times: 
                - "2019-11-07T14:10:49.35"
                - "2019-11-07T14:10:59.55"
                - "2019-11-07T14:10:39.53"
```

#### 2: Comparing several stations

```yaml
---
input:
    starttime: null
    endtime: null
    datafiles:
        - 
            name: "20191107T14_SPOBS09_F02.raw.lch"
            obs_type: "SPOBS2"
            station: "09F2"
        - 
            name: "20191107T14_SPOBS09_F02.raw.lch"
            obs_type: "SPOBS2"
            station: "09c1"
        - 
            name: "20191107T14_SPOBS09_F02.raw.lch"
            obs_type: "SPOBS2"
            station: "09c2"
    description: "Simulation of multi-instrument test"
output:
    show: True
    filebase: "MAYOBS6"
plot_globals:
    stack:
        offset_before.s: 0.5
        offset_after.s:  1.5
        plot_span: False
    particle_motion:
        offset_before.s: 0.00
        offset_after.s: 0.03
        offset_before_ts.s: 0.1
        offset_after_ts.s: 0.2
    spectra:
        window_length.s: 100
plots:
    time_series:
        -
            description: "Entire time series"
            select: {station: "*"}
            start_time: null
            end_time: null
        -
            description: "Quiet period"
            select: {channel: "*3"}
            start_time: null
            end_time: "2019-11-07T13:57"
        -
            description: "Rubber hammer taps"
            select: {station: "*"}
            start_time: "2019-11-07T14:08"
            end_time: "2019-11-07T14:11:10"
    spectra:
        -
            description: "Entire time series"
            select: {component: "3"}
            start_time: null
            end_time: null
        -
            description: "Quiet period"
            select: {channel: "*3"}
            start_time: null
            end_time: "2019-11-07T13:57"
```
