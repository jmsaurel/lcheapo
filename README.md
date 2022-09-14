# lcheapo

Viewing and modifying LCHEAPO OBS data

## Overview

### Command-line programs

Type ``{command} -h`` to get a list of parameters and options

#### Programs that don't modify files

| Program     | description                                           |
| ----------- | ----------------------------------------------------- |
| lcdump      | dump raw information from LCHEAPO files               |
| lcinfo      | return basic information about an LCHEAPO file        |
| lcplot      | plot an LCHEAPO file                                  |
| lctest      | plot LCHEAPO tests                                    |
| lc_examples | create a directory with examples of lcplot and lctest |

#### Programs that modify files

These programs use the *SDPCHAIN* protocols for FAIR-
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
output.  The datetime ranges in the `plots` sections must be within those given
in the `input` section, as data is only read using the `input` time bounds.

for details on the control file format, type:
```
  python
    > import lcheapo
    > help(lcheapo.lctest)
```

to put example lctest control files in the current directory, type:
```
    lctest --examples
```

### Example plots

### Examples

#### 1: Analysing one station

``` yaml
---
input: 
    start_time: "2022-02-22T10:00:01"
    end_time: "2022-02-25T15:25:25"
    datafiles:
        -   name: "TestAcq-BB02-ProgV1-3.raw.lch"
            obs_type: 'BBOBS1'
            station: 'TEST'
    description: "Tests on BBOBS"
output:
    show: True
    filebase: 'BB02-V1_3-tests'
plot_globals:
    spectra:
        window_length.s: 1024
plots:
    time_series:
        -   description: "Entire time series"
            select: {station: "*"}
            start_time: "2022-02-22T10:00:01"
            end_time: "2022-02-25T15:25:25"
        -   description: "Quiet time"
            select: {station: "*"}
            start_time: "2022-02-23T21:00:00"
            end_time: "2022-02-24T03:00:00"
        -   description: "Stack time"
            select: {station: "*"}
            start_time: "2022-02-25T13:54:00"
            end_time: "2022-02-25T14:03:00"
    spectra:
        -   description: "Quiet time"
            select: {station: "*"}
            start_time: "2022-02-23T21:00:00"
            end_time: "2022-02-24T03:00:00"
    stack:
        -   description: "Stack, Jump South"
            orientation_codes: ["Z"]
            offset_before.s: 0.3
            offset_after.s: 1
            times:
            -    "2022-02-25T13:57:00.66"
            -    "2022-02-25T13:58:00.53"
            -    "2022-02-25T13:59:00.2"
        -   description: "Stack, Jump Est"
            orientation_codes: ["Z"]
            offset_before.s: 0.3
            offset_after.s: 1
            times:
            -    "2022-02-25T14:00:00.4"
            -    "2022-02-25T14:01:00.15"
            -    "2022-02-25T14:02:00.18"
    particle_motion:
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
        -   description: "Stack, Jump Est"
            orientation_code_x: "2"
            orientation_code_y: "1"
            offset_before.s: 0.1
            offset_after.s: 0.2
            offset_before_ts.s: 0.3
            offset_after_ts.s: 1
            times:
            -    "2022-02-25T14:00:00.4"
            -    "2022-02-25T14:01:00.15"
            -    "2022-02-25T14:02:00.18"
```
##### Output plots
###### time_series
![](https://github.com/WayneCrawford/lcheapo/raw/main/README_images/BB02-V1_3-tests_Entire_time_series_ts.png)
![](https://github.com/WayneCrawford/lcheapo/raw/main/README_images/BB02-V1_3-tests_Quiet_time_ts.png)

###### spectra
![](https://github.com/WayneCrawford/lcheapo/raw/main/README_images/BB02-V1_3-tests_Quiet_time_spect.png)

###### stack
![](https://github.com/WayneCrawford/lcheapo/raw/main/README_images/BB02-V1_3-tests_Stack_Jump_South_stack.png)

###### particle_motion
![](https://github.com/WayneCrawford/lcheapo/raw/main/README_images/BB02-V1_3-tests_Stack_Jump_South_pm.png)


#### 2: Comparing several stations

```yaml
---
input:
    start_time: null
    end_time: null
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
