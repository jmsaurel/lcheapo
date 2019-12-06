"""
Plot LCHEAPO test results

Reads parameters from a lctest.yaml file
Plots:
    - time_series of one or more stations/channels
    - stacks of time series from one station/channel
    - spectra from  multiple stations/channels
    - particle_motions between two channels
"""
from lc_read import read
from obspy.core import UTCDateTime, Stream
from obspy.signal import PPSD
import yaml
import sys
import argparse
from matplotlib import pyplot as plt


def main():
    """
    Read in the yaml file and plot the specified tests
    """
    args = get_arguments()
    root = read_lctest_yaml(args.yaml_file)
    stack_offset_before, stack_offset_after, stack_plot_span = get_globals(root)
    stream = read_files(root['datafiles'])
    for test in root['tests']:
        plot_type = test['plot_type']
        if plot_type == 'time_series':
            try :
                starttime = UTCDateTime(test['start_time'])
            except:
                starttime = None
            try :
                endtime = UTCDateTime(test['end_time'])
            except:
                endtime = None
            plot_time_series(stream, starttime, endtime, test['select'],
                             description=test['description'])
        elif plot_type == 'stack':
            for o_code in test['orientation_codes']:
                trace = stream.select(component=o_code)[0]
                plot_stack(trace,
                           [UTCDateTime(t) for t in test['times']],
                           test['description'],
                           offset_before=test.get('offset_before',
                                                  stack_offset_before),
                           offset_after=test.get('offset_after',
                                                 stack_offset_after),
                           plot_span=test.get('plot_span',stack_plot_span))
        elif plot_type == 'spectra':
            print('NOT YET IMPLEMENTED')
            # inv = read_inventory(arguments.sta_file)
        elif plot_type == 'particle_motion':
            print('NOT YET IMPLEMENTED')
        else:
            print(f'UNKNOWN PLOT TYPE "{plot_type}", skipping...')


def get_globals(root):
    """
    Return global values
    
    For now, just offset_before and offset_after
    """
    offset_before, offset_after = None, None
    plot_span = False
    if 'globals' in root:
        if 'stack' in root['globals']:
            offset_before = root['globals']['stack'].get('offset_before', None)
            offset_after = root['globals']['stack'].get('offset_after', None)
            plot_span = root['globals']['stack'].get('plot_span', False)
    return offset_before, offset_after, plot_span


def read_files(datafile_list):
    """
    Read in data from one or more datafiles

    datafile_list: list of dict(name=, obs_type=, station=)
    """
    stream = Stream()
    for df in datafile_list:
        s = read(df['name'], chan_map=df['obs_type'])
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


def plot_time_series(in_stream, starttime=None, endtime=None,
                     selectargs=None, description=None, filebase=None):
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
        stream.select(**selectargs)
    outfile = None
    if filebase:
        outfile = filebase + '.png'

    fig=plt.figure()
    fig.suptitle(description)
    stream.plot(size=(800, 600), outfile=outfile, fig=fig)
    plt.show()


def plot_stack(trace, times, description, outfile=None,
               offset_before=0.5, offset_after=1.5, plot_span=False):
    """
    Plot a stack of time series from one trace

    trace is an obspy Trace object
    title is the title to put on the graph
    outfile is the filename to write the plot out to (None = print to screen)
    times is a list of UTCDateTimes
    plot_span: whether to plot a time series spanning the first to last time
    """
    title = description + f', orientation_code="{trace.stats.channel[-1]}"'
    if plot_span:
        first = min(times)
        last = max(times)
        span = last-first
        trace.slice(first-span/2,last+span/2).plot()
    print(f'Plotting stack "{title}"')
    colors = 'rbgmcyk'  # should be more intelligent about color cycling
    i = 0               # then won't need "i"
    offset_vertical = 0
    # time_zero = UTCDateTime(times[0])
    ax = plt.axes()
    max_val = 0
    for time in times:
        temp=trace.slice(time - offset_before, time + offset_after)
        if abs(temp.max()) > max_val:
            max_val = abs(temp.max())
    for time in times:
        offset_vertical += max_val
        t = trace.slice(time - offset_before, time + offset_after)
        ax.plot(t.times("utcdatetime") - time,
                t.data - offset_vertical,
                colors[i % len(colors)],
                label=time.strftime('%H:%M:%S') +
                '.{:02d}'.format(int(time.microsecond / 1e4)))
        offset_vertical += max_val
        i += 1
    ax.set_title(title)
    ax.grid()
    ax.legend()
    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()


def read_lctest_yaml(filename='_examples/events_secondtest.yaml'):
    """
    Verify and read in an lctest yaml file
    """
    print('VERIFY LCTEST.YAML NOT YET IMPLEMENTED!')
    with open(filename) as f:
        root = yaml.safe_load(f)
    return root


def plot_PPSD(trace, sta, start_time, interval=7200):
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

    filebase = '{}.{}.{}.{}'.format(trace.stats.network_code,
                                    trace.stats.station,
                                    trace.stats.location_code,
                                    trace.stats.component)
    ppsd.save_npz(f'{filebase}_PPSD.npz')
    ppsd.plot(f'{filebase}_PPSD.png')
    # ppsd.plot_temporal([0.1,1,10])
    # ppsd.plot_spectrogram()
    return 0


if __name__ == '__main__':
    sys.exit(main())
