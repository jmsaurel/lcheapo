"""
Plot LCHEAPO test results

Reads parameters from a lctest.yaml file
Plots:
    - time_series of one or more stations/channels
    - stacks of time series from one station/channel
    - spectra from  multiple stations/channels
    - particle_motions between two channels
"""
import numpy as np
import yaml
import sys
import re
import argparse
import pkg_resources
import os
import json

import jsonref
import jsonschema
from matplotlib import pyplot as plt
from obspy.core import UTCDateTime, Stream
from obspy.signal import PPSD

from .lcread import read as lcread
from .spectral import PSDs


def main():
    """
    Read in the yaml file and plot the specified tests
    """
    args = get_arguments()
    root = read_lctest_yaml(args.yaml_file)
    plot_globals = get_plot_globals(root)
    stream = read_files(root['input'])
    show = root['output']['show']
    filebase = root['output']['filebase']
    if 'time_series' in root['plots']:
        for plot in root['plots']['time_series']:
            plot_time_series(stream,
                             _UTCDateTimeorNone(plot['start_time']),
                             _UTCDateTimeorNone(plot['end_time']),
                             plot['select'],
                             description=plot['description'],
                             filebase=filebase,
                             show=show)
    if 'stack' in root['plots']:
        for plot in root['plots']['stack']:
            globs = plot_globals.get('stack', None)
            print(f'globs={globs}')
            print(f'plot={plot}')
            for o_code in plot['orientation_codes']:
                trace = _trace_component(stream, o_code)
                plot_stack(
                    trace,
                    [UTCDateTime(t) for t in plot['times']],
                    plot['description'],
                    offset_before=_get_val('offset_before.s', plot, globs),
                    offset_after=_get_val('offset_after.s', plot, globs),
                    plot_span=_get_val('plot_span', plot, globs),
                    filebase=filebase,
                    show=show)
    if 'particle_motion' in root['plots']:
        for plot in root['plots']['particle_motion']:
            tracex = _trace_component(stream, plot['orientation_code_x'])
            tracey = _trace_component(stream, plot['orientation_code_y'])
            # tracex = stream.select(component=plot['orientation_code_x'])[0]
            # tracey = stream.select(component=plot['orientation_code_y'])[0]
            globs = plot_globals['particle_motion']
            plot_particle_motion(
                tracex, tracey, [UTCDateTime(t) for t in plot['times']],
                plot['description'],
                offset_before=_get_val('offset_before.s', plot, globs),
                offset_after=_get_val('offset_after.s', plot, globs),
                offset_before_ts=_get_val('offset_before_ts.s', plot, globs),
                offset_after_ts=_get_val('offset_after_ts.s', plot, globs),
                plot_span=_get_val('plot_span', plot, globs),
                filebase=filebase,
                show=show)
    if 'spectra' in root['plots']:
        for plot in root['plots']['spectra']:
            spect = calc_spect(stream, plot,
                               plot_globals.get('spectra', None))
            if filebase:
                spect.plot(outfile=filebase + '_'
                           + get_valid_filename(plot['description'])
                           + '_spect.png')
            if show:
                spect.plot()
            # inv = read_inventory(arguments.sta_file)


def _trace_component(stream, component):
    """ Return the first trace corresponding to the given component """
    try:
        return stream.select(component=component)[0]
    except IndexError:
        print(f'Did not find component {component}')
        sys.exit()


def calc_spect(stream, plot_parms, glob_parms=None):
    """
    Calculate spectra for in_stream
    """
    # Select appropriate channels
    if plot_parms['select']:
        stream = stream.select(**plot_parms['select'])
    starttime = _UTCDateTimeorNone(plot_parms['start_time'])
    endtime = _UTCDateTimeorNone(plot_parms['end_time'])
    if starttime or endtime:
        stream = stream.slice(starttime=starttime, endtime=endtime)
    # Calculate spectra
    wl = _get_val('window_length.s', plot_parms, glob_parms)
    if wl:
        spect = PSDs.calc(stream, window_length=wl)
    else:
        spect = PSDs.calc(stream)
    return spect


def plot_time_series(in_stream, starttime=None, endtime=None,
                     selectargs=None, description=None, filebase=None,
                     show=True):
    """
    Plot a time series

    in_stream: obspy stream object
    starttime: starttime (UTCDateTime or None)
    endtime: endtime (UTCDateTime or None)
    selectargs: dictionary of kwargs to pass to stream.select()
    filebase: base filename
    """
    print(f'Plotting time series "{description}"')
    stream = in_stream.slice(starttime=starttime, endtime=endtime)
    if selectargs:
        stream = stream.select(**selectargs)
    outfile = None
    if filebase:
        outfile = filebase

    fig = plt.figure()
    fig.suptitle(description)

    outfile = False
    if filebase:
        outfile = filebase + '_' + get_valid_filename(description) + '_ts.png'
        stream.plot(size=(800, 600), fig=fig, outfile=outfile)
    if show:
        stream.plot(size=(800, 600), fig=fig)
        plt.show()


def plot_stack(trace, times, description, offset_before=0.5, offset_after=1.5,
               plot_span=False, filebase=None, show=True):
    """
    Plot a stack of time series from one trace

    trace is an obspy Trace object
    title is the title to put on the graph
    outfile is the filename to write the plot out to (None = print to screen)
    times is a list of UTCDateTimes
    plot_span: whether to plot a time series spanning the first to last time
    """
    title = description + f', orientation_code="{trace.stats.channel[-1]}"'
    assert offset_before >= 0,\
        f'plot_stack "{description}": offset_before < 0 ({offset_before:g})'
    assert offset_after > 0,\
        f'plot_stack "{description}": offset_after <= 0 ({offset_after:g})'
    if plot_span:
        _plot_span(times, Stream(traces=[trace]))
    print(f'Plotting stack "{title}"')
    colors = plt.cm.rainbow(np.linspace(0, 1, len(times)))
    offset_vertical = 0
    # time_zero = UTCDateTime(times[0])
    ax = plt.axes()
    max_val = 0
    # Set up y axis range
    for time in times:
        temp = trace.slice(time - offset_before, time + offset_after)
        if abs(temp.max()) > max_val:
            max_val = abs(temp.max())
    # Plot the subtraces
    for time, c in zip(times, colors):
        offset_vertical += max_val
        t = trace.slice(time - offset_before, time + offset_after)
        ax.plot(t.times("utcdatetime") - time,
                t.data - offset_vertical, color=c,
                label=time.strftime('%H:%M:%S') +
                '.{:02d}'.format(int(time.microsecond / 1e4)))
        offset_vertical += max_val
    ax.set_title(title)
    ax.grid()
    ax.legend()
    if filebase:
        plt.savefig(filebase + '_' + get_valid_filename(description) +
                    '_stack.png')
    if show:
        plt.show()


def plot_particle_motion(tracex, tracey, times, description,
                         offset_before=0.0, offset_after=0.5,
                         offset_before_ts=0.5, offset_after_ts=1.5,
                         plot_span=False, filebase=None, show=True):
    """
    Plot particle motions

    tracex is an obspy Stream object to plot on x axis
    tracey is an obspy Stream object to plot on y axis
    title is the title to put on the graph
    outfile is the filename to write the plot out to (None = print to screen)
    times is a list of UTCDateTimes
    offset_before, offset_after: seconds before and after "time" to plot
    offset_before_ts, offset_after_st: seconds before and after "time"
                                       to plot as time series
    plot_span: whether to plot a time series spanning the first to last time
    """
    tx_comp = tracex.stats.channel[-1]
    ty_comp = tracey.stats.channel[-1]
    title = description + f', compare "{tx_comp}" to "{ty_comp}"'
    if plot_span:
        _plot_span(times, Stream(traces=[tracex, tracey]))
    print(f'Plotting particle motion "{title}"')
    # Setup axis grid
    fig = plt.figure()
    gs = fig.add_gridspec(2, 3, hspace=0, wspace=0)
    # two columns, one row:
    axx = fig.add_subplot(gs[0, :2])
    axy = fig.add_subplot(gs[1, :2], sharex=axx, sharey=axx)
    # one column, one row
    axxy = fig.add_subplot(gs[1, 2], sharey=axy)
    print(type(offset_before_ts),type(offset_after_ts))
    for time in times:
        # time series plots
        _plot_one_ts(axx, tracex, tx_comp, time, offset_before_ts, offset_after_ts)
        _plot_one_ts(axy, tracey, ty_comp, time, offset_before_ts, offset_after_ts)

        # partical motion plot
        tx = tracex.slice(time - offset_before, time + offset_after)
        ty = tracey.slice(time - offset_before, time + offset_after)
        axxy.plot(tx.data, ty.data)
        axxy.axvline(0)
        axxy.axhline(0)
        axxy.set_aspect('equal', 'datalim')
        axxy.set_xlabel(tx_comp)
        # axxy.set_yticklabels([])
    fig.suptitle(title)
    if filebase:
        plt.savefig(filebase + '_' + get_valid_filename(description) +
                    '_pm.png')
    if show:
        plt.show()


def _plot_one_ts(ax, trace, comp, time, offset_before, offset_after):
    t = trace.slice(time - offset_before, time + offset_after)
    ax.plot(t.times("utcdatetime") - time, t.data)
    ax.axvline(-offset_before, color='k', linestyle='--')
    ax.axvline(offset_after, color='k', linestyle='--')
    ax.set_ylabel(comp)


def plot_PPSD(trace, sta, start_time, interval=7200, filebase=None,
              show=True):
    """
    Plot a Probabilistic Power Spectral Desnsity for the trace

    trace = obspy Trace objet
    sta = obspy Inventory/Station object corresponding to the trace
    start_time = time at which to start spectra
    interval=offset between PSDs (seconds, minimum=3600)
    """
    now_time = trace.start_time
    first_read = True
    while now_time < trace.end_time-interval:
        if first_read:
            if trace.stats.component[1] == 'D':
                ppsd = PPSD(trace.stats, metadata=sta,
                            special_handling='hydrophone')
            else:
                ppsd = PPSD(trace.stats, metadata=sta)
            first_read = False
        ppsd.add(trace)
        now_time += interval

    ppsd.save_npz(f'{filebase}_PPSD.npz')
    if filebase:
        description = '{}.{}.{}.{}'.format(trace.stats.network,
                                           trace.stats.station,
                                           trace.stats.location,
                                           trace.stats.channel)
        ppsd.plot(filebase + '_' + description
                  + '_PPSD.png')
    if show:
        plt.plot()
    # ppsd.plot_temporal([0.1,1,10])
    # ppsd.plot_spectrogram()
    return 0


def _UTCDateTimeorNone(input):
    """
    Converts input string to UTCDateTime or None if it doesn't work
    """
    try:
        return UTCDateTime(input)
    except TypeError:
        return None


def _get_val(key, dict1, dict2):
    """
    Returns key value from dict1 or dict2

    If it's in dict1, returns from dict1,
    If it's not in dict1 but in dict2, returns from dict2
    If it's in neither, returns False
    """
    if key in dict1:
        return dict1[key]
    elif dict2:
        if key in dict2:
            return dict2[key]
    return False


def get_plot_globals(root):
    """
    Return global plot values
    """
#     globals = {'stack': {'offset_before.s': False,
#                          'offset_after.s': False,
#                          'plot_span': False},
#                'particle_motion': {'offset_before.s': False,
#                                    'offset_after.s': False,
#                                    'offset_before_ts.s': False,
#                                    'offset_after_ts.s': False,
#                                    'plot_span': False}}
#     if 'plot_globals' in root:
#         plot_globals = root['plot_globals']
#         for type in globals:
#             if type in plot_globals:
#                 for parm in globals[type]:
#                     if parm in plot_globals[type]:
#                         globals[type][parm] = plot_globals[type][parm]
    globals = root.get('plot_globals',None)
    return globals


def read_files(inputs):
    """
    Read in data from one or more datafiles

    datafile_list: list of dict(name=, obs_type=, station=)
    """
    stream = Stream()
    for df in inputs['datafiles']:
        s = lcread(df['name'],
                   starttime=inputs.get('starttime', False),
                   endtime=inputs.get('endtime', False),
                   obs_type=df['obs_type'])
        for t in s:
            t.stats.station = df['station']
        stream += s
    return stream


def get_arguments():
    """
    Get command line arguments
    """
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('yaml_file', help='YAML parameter file')
    p.add_argument("-v", "--verbose", default=False, action='store_true',
                   help="verbose")
    p.add_argument("--inv_file", help="inventory file")
    return p.parse_args()


def get_valid_filename(s):
    s = str(s).strip().replace(' ', '_')
    return re.sub(r'(?u)[^-\w.]', '', s)


def _plot_span(times, stream):
    first = min(times)
    last = max(times)
    span = last-first
    stream.slice(first - span / 2, last + span / 2).plot()


def read_lctest_yaml(filename):
    """
    Verify and read in an lctest yaml file
    """
    with open(filename) as f:
        root = yaml.safe_load(f)
        if not validate_schema04(root, 'lctest'):
            sys.exit()
    return root


def validate_schema04(instance: dict, type: str = 'lctest'):
    """
    Validates a data structure against a draft 04 jsonschema

    :param instance: data structure read from a yaml or json filebase
    :param type: type of the data structure
    """
    schema_file = pkg_resources.resource_filename(
        "lcheapo_obspy", f"data/{type}.schema.json")
    base_path = os.path.dirname(schema_file)
    base_uri = f"file:{base_path}/"
    with open(schema_file, "r") as f:
        try:
            schema = jsonref.loads(f.read(), base_uri=base_uri,
                                   jsonschema=True)
        except json.decoder.JSONDecodeError as e:
            print(f"JSONDecodeError: Error loading schema file: {schema_file}")
            print(str(e))
            return False
        # except:
        #     print(f"Error loading JSON schema file: {schema_file}")
        #     print(sys.exc_info()[1])
        #     return False

    # Lazily report all errors in the instance
    # ASSUMES SCHEMA IS DRAFT-04
    try:
        v = jsonschema.Draft4Validator(schema)

        if not v.is_valid(instance):
            print("")
            for error in sorted(v.iter_errors(instance), key=str):
                print("\t\t", end="")
                for elem in error.path:
                    print(f"['{elem}']", end="")
                print(f": {error.message}")
            print("\tFAILED")
            return False
        else:
            print("OK")
            return True
    except jsonschema.ValidationError as e:
        print("")
        print("\t" + e.message)
        return False


if __name__ == '__main__':
    sys.exit(main())
