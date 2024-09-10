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
| lctest      | plot LCHEAPO tests  (**use obstest instead**)         |
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
