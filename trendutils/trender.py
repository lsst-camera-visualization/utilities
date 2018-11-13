#!/usr/bin/env python
""" trending data app: gets specified channels for requested time period
"""
import re
import os
import argparse
import textwrap
import logging
import sys
import time
import datetime as dt
import requests
from lxml import etree
from dateutil.tz import gettz
import numpy as np
from astropy import stats
from astropy.time import Time
import matplotlib.pyplot as plt
import trendutils as tu

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Fetch trending data from camera database
                                    '''),
        epilog=textwrap.dedent('''
        This application expects to access the CCS trending database,
        either directly, at lsst-mcm:8080, or via an ssh tunnel which
        must be set up in advance to use localhost:8080.

        Alternately trending data can come from a local file saved
        in an earlier run.  See the "input_file" options.  This allows
        multiple output runs (stats, plots, etc.) using the same trending
        data w/out re-querying the server.

        The "channel_source"(s) are either files or patterns (re's).  Files
        are best constructed using the sibling application
        "trendingChannels.py" and editing the resuling files to choose
        which channels to activate.

        The interval and start/stop date specifiers an be nearly any
        variant of a full date/time spec as output by the unix
        "date <options>" command.  One particular choice is the format
        from "date --iso-8601=seconds".
                               '''))
    # Input args
    parser.add_argument("channel_source", nargs='+',
                        help="ListFile|Pattern specifying channels to report")
    parser.add_argument("--input_file", nargs='+',
                        help="XML file with trending data, =>no db query")
    # Output options
    ogroup = parser.add_mutually_exclusive_group()
    ogroup.add_argument("--xml", action='store_true',
                        help="Print formatted xml from trending to stdout")
    ogroup.add_argument("--text", action='store_true',
                        help="Print (timestamp, value, path) colum text")
    ogroup.add_argument("--stats", action='store_true',
                        help="Print statistics for each channel")
    ogroup.add_argument("--matchingchannels", action='store_true',
                        help="print list of matching channels and exit")
    # Plotting specifications
    parser.add_argument("--plot", action='store_true',
                        help="produce plotting")
    parser.add_argument("--overlaytime", action='store_true',
                        help="all time intervals on same plot")
    parser.add_argument("--overlayunits", action='store_true',
                        help="all units on same plot")
    # Time interval specifications
    dgroup = parser.add_mutually_exclusive_group()
    dgroup.add_argument("--interval", metavar=('dstart', 'dstop'),
                        nargs=2, action='append',
                        help="Pair of date/time specifiers")
    dgroup.add_argument("--start", metavar='tstart',
                        help="Date/time specifier(s)", nargs='+')
    dgroup.add_argument("--stop", metavar='tstop',
                        help="Date/time specifier(s)", nargs='+')
    parser.add_argument("--duration", metavar='seconds',
                        default=None, type=int,
                        help="duration in seconds for start/stop spec")
    parser.add_argument("--debug", action='store_true',
                        help="Print additional debugging info")
    return parser.parse_args()


def main():
    """main logic"""
    # get command args and options
    optlist = parse_args()

    # logging level
    if optlist.debug:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.DEBUG)
        logging.debug("printing extra debugging...")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.INFO)

    logging.debug('optlist: %s', optlist)

    # get list of time intervals to process
    intervals = get_unique_time_intervals(optlist)
    if not intervals:
        logging.error('time interval spec failed')
        exit(1)
    # deal with interval accounting
    intcnt = len(intervals)
    inttot = int(sum([t[1] - t[0] for t in intervals])/1000)
    tmin = intervals[0][0]
    tmax = intervals[-1][1]

    trending_server = tu.get_trending_server()
    if trending_server:
        data_url = "http://{}:8080/rest/data/dataserver".format(
            trending_server)
    else:
        logging.error('failed to contact trending server')
        exit(1)

    # from https://stackoverflow.com/questions/16710076
    # regex to split a string preserving quoted fields
    #
    rpat = re.compile(r"""             #
                      (?:[^\s"']+)|    # match non-delimiter
                      "(?:\\.|[^"]*)"| # match double quoted
                      '(?:\\.|[^']*)'  # match single quoted
                      """, re.X)

    # loop over input sources to define the channels for query/output
    # 2 choices:
    #     1.) a file with 3 fields per line (after comment removal) where
    #         the line format is '^[01]\s+<CCS path>\s+<channel_id>'
    #     2.) a pattern that can be matched against the cached full list of
    #         channels in ~/.trender/listchannels.xml
    #
    # After evaluating input, channels needed are stored in "oflds"
    # oflds[channel_id_num] = trending_full_path
    # eg: oflds[2372] = aliveness-raft/R00.Reb2.RGL
    #
    oflds = dict()  # dict to hold channel information
    #
    channel_source_type = None
    chid_dict = None
    for csource in optlist.channel_source:
        logging.debug('channel_source= %s', csource)
        #  test to determine type of channel_source
        #
        if channel_source_type == 'file' or not channel_source_type:
            logging.debug('test for formatted input file...')
            try:
                cf = open(csource, mode='r')
            except OSError as e:
                logging.debug('open(%s) failed: %s', csource, e)
            else:
                # populate oflds[id] with the corresponding channel path
                for line in cf:
                    if re.match(r"^\s*#", line):  # skip block comment
                        continue
                    if re.match(r"^\s*$", line):  # skip white space line
                        continue
                    # strip inline cmnt
                    sline = re.sub(r"""(#[^\'^"]*$)""", '', line)
                    # tokenize what remains
                    flds = [''.join(t) for t in rpat.findall(sline)]
                    if len(flds) != 3:
                        logging.warning('bad input line: %s', line)
                        continue
                    if int(flds[0]) == 1:
                        oflds[flds[2]] = flds[1]  # oflds[id] = path
                cf.close()
            if len(oflds):  # only set if some lines were good
                channel_source_type = 'file'

        if channel_source_type == 'pattern' or not channel_source_type:
            logging.debug('eval pattern for matching channels...')
            # 0: datachannels [-]
            #   1: datachannel [-]
            #     2: path [-]
            #       3: pathelement [-]
            #   2: id [-]
            #   2: metadata [name, value]
            channel_file = tu.update_trending_channels_xml()
            if not chid_dict:
                # build a channel dictionary for all channels
                # just do this once
                # this should be a function of course
                logging.debug('building full channel dictionary')
                chid_dict = dict()
                path_dict = dict()
                tree = etree.parse(channel_file)
                root = tree.getroot()
                for dchan in root.iterfind('datachannel'):
                    chid = dchan.find('id').text
                    # build path
                    parr = []  # list to hold path elements
                    pp = dchan.find('path')
                    for pe in pp.iterfind('pathelement'):
                        parr.append(pe.text)
                    path = '/'.join(parr)
                    if path and chid:  # this is supposed to be 1-1 map
                        chid_dict[chid] = path
                        path_dict[path] = chid

            # search the entire catalog for each pattern, this may be slow
            cpat = re.compile(csource)
            for chid in chid_dict:
                if cpat.search(chid_dict[chid]):
                    oflds[chid] = chid_dict[chid]

            if len(oflds):
                channel_source_type = 'pattern'

    # end of loop over channel sources

    if optlist.debug:
        logging.debug('Found matching channels:')
        for chid in oflds:
            logging.debug('id= %5d  path= %s', int(chid), oflds[chid])
        if optlist.matchingchannels:
            exit(0)

    if optlist.matchingchannels:
        print('Found matching channels:')
        for chid in oflds:
            print("   id: {}  path: {}".format(chid, oflds[chid]))
        exit(0)

    if optlist.debug:
        logging.debug('Found matching channels:')
        for chid in oflds:
            logging.debug('id= %5d  path= %s', int(chid), oflds[chid])

    if optlist.input_file:
        # get input from these files rather than trending service
        # an issue is that the input file need not have the
        # same set of time intervals as given on the command line.
        # We expect to apply the intervals to the output only by taking
        # only times in the intersection
        responses = []
        for ifile in optlist.input_file:
            logging.debug('using %s for input', ifile)
            logging.debug('test for well-formed xml...')
            try:
                tree = etree.parse(ifile)
            except etree.ParseError as e:
                logging.debug('parsing %s failed: %s', ifile, e)
                exit(1)
            except etree.XMLSyntaxError as e:
                logging.debug('parsing %s failed: %s', ifile, e)
                exit(1)
            else:
                logging.debug('successfully parsed %s', ifile)

            responses.append(etree.tostring(tree.getroot(), encoding="UTF-8",
                xml_declaration=True, pretty_print=False))
    else:
        # put the rest server query responses into a list
        # join the ids requested as "id0&id=id1&id=id2..." for query
        idstr = '&id='.join(id for id in oflds)
        responses = []
        for interval in intervals:
            res = query_rest_server(interval[0], interval[1], data_url, idstr)
            responses.append(res)

    # Output a well formed xml tree aggregating all the xml received
    # Main use is to save to local file, and re-use for subsequent queries
    # for statistics, plots etc. with subset of channels and time periods
    # Also useful fo debugging and verification of data
    if optlist.xml:
        xml_dec = b'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        os.write(1, xml_dec)
        os.write(1, b'<datas>\n')
        for res in responses:
            root = etree.fromstring(res)
            for data in root.iter('data'):
                os.write(1, etree.tostring(
                    data, encoding="UTF-8",
                    xml_declaration=False, pretty_print=True))
        os.write(1, b'</datas>')
        exit(0)

    # Translate the xml responses into internal arrays etc.
    # XML Tree structure looks like this:
    # 1: data [id, path]
    # 2: trendingresult [-]
    #     3: channelmetadata [-]
    #         4: channelmetadatavalue [tstart, tstop, name, value]
    #     3: trendingdata [-]
    #         4: datavalue [name, value]
    #         4: axisvalue [name, value, loweredge, upperedge]
    # where [id, path] could appear multiple times and input time intervals are
    # allowed to overlap
    #
    chanspec = dict()  # where keys are chids, element is also a dict
    chanmd = dict()    # key is chid, elements will be dicts holding arrays
    chandata = dict()  # key is chid, element is list of (time, value) tuples
    datacnt = 0
    for res in responses:
        root = etree.fromstring(res)
        for data in root.iter('data'):
            datacnt += 1
            chid = data.attrib.get('id')
            path = data.attrib.get('path')
            # verify this element's (chid, path) matches the input list
            # logging.debug('id=%s  path=%s', chid, path)
            if chid not in oflds:
                continue
            if oflds[chid] != path:
                logging.error('inputpath(id=%s): %s != %s = xmlpath',
                              chid, oflds[chid], path)
            # first check for existing....
            if chid in chanspec and chanspec[chid]['path'] != path:
                logging.warning('path mismatch for channel_id= %d', chid)
                logging.warning('  %s != %s, skipping....',
                                chanspec[chid]['path'], path)
            else:
                chanspec[chid] = dict()
                chanspec[chid]['path'] = path
                chanspec[chid]['units'] = "none"

            # channelmetadata:
            # each element is a name, value and time interval
            # a name can appear multiple times with distinct time intervals
            # convert to a list, per name, of ordered pairs (value,time)
            # that could be plotted using those points
            #
            if chid not in chanmd:  # check if already exists
                chanmd[chid] = dict()
            for mdval in data.iter('channelmetadatavalue'):
                if mdval.keys():  # empty sequence is false
                    mdname = mdval.attrib.get('name')  # key
                    mdvalue = mdval.attrib.get('value')  # value
                    mdstart = mdval.attrib.get('tstart')
                    mdstop = mdval.attrib.get('tstop')
                if mdname in chanmd[chid]:
                    chanmd[chid][mdname].append((mdstart, mdvalue))
                    chanmd[chid][mdname].append((mdstop, mdvalue))
                else:  # first assignment
                    chanmd[chid][mdname] = [(mdstart, mdvalue),
                                            (mdstop, mdvalue)]
            # trendingdata:
            # extract timestamp, value pairs in axisvalue, datavalue tags
            if chid not in chandata:  # first time
                chandata[chid] = []  # empty list
            for tdval in data.iter('trendingdata'):
                dataval = tdval.find('datavalue')
                if dataval is not None:
                    tvalue = dataval.attrib.get('value')
                else:
                    continue
                axisval = tdval.find('axisvalue')
                if axisval is not None:
                    tstamp = axisval.attrib.get('value')
                else:
                    continue
                # if tstamp is in intervals then append
                for ival in intervals: # slow, but no other way?
                    if ival[0] < int(tstamp) < ival[1]:
                        chandata[chid].append((tstamp, tvalue))
                        break

    # done with responses, delete all the raw xml responses
    logging.debug('processed %d xml channel responses', len(responses))
    logging.debug('processed %d uniq channel requests', len(chanspec))
    logging.debug('processed %d total channel queries', datacnt)
    del responses

    # chanspec = dict()  # where keys are chids, values are ccs paths
    # chanmd = dict()  # key is chid, elements will be dicts holding arrays
    # chandata = dict() # key is chid, elements are (time, value) pair arrays
    # so all responses processed, now have data organized by a set of dicts
    # with the the index on channel id.  Multiple queries for a given channel
    # id are grouped together and there could be duplicate values.
    #
    # To facilitate operating on the data, transform from list[] based data
    # (which was easy to append to) to np.array based data.
    chandt = np.dtype({'names': ['tstamp', 'value'],
                       'formats': ['int', 'float']})
    for chid in chanspec:
        path = chanspec[chid]['path']
        logging.debug('id=%s  path=%s', chid, path)
        for mdname in chanmd[chid]:
            # pick out and process the md's we want
            if mdname == 'units' and chanspec[chid]['units'] == 'none':
                unitstr = chanmd[chid][mdname][-1][1]  # use last one
                chanspec[chid]['units'] = unitstr

        logging.debug('    units=%s', chanspec[chid]['units'])
        # sort and remove duplicates from chandata[chid] where:
        # chandata[chid] = [(t0, v0), (t1, v1), ...., (tn, vn)]
        # and convert to np array
        tmparr = np.array(chandata[chid], dtype=chandt)
        chandata[chid] = np.unique(tmparr)
        logging.debug('    chandata: %d uniq/sorted values from %d entries',
                      np.size(chandata[chid]), np.size(tmparr))

    # print to stdout a text dump of the data, in time order per channel
    #
    if optlist.text:
        # print a header for the text
        #
        print("# CCS trending dump at {}".format(
            dt.datetime.now(gettz()).isoformat(timespec='seconds')))
        print("# Data for {} total seconds from {} intervals".format(
            inttot, intcnt), end="")
        print(" over {} (h:m:s) from:".format(
            dt.timedelta(seconds=(tmax/1000 - tmin/1000))))
        print("#     tmin={}: \"{}\"".format(tmin, dt.datetime.fromtimestamp(
            tmin/1000, gettz()).isoformat(timespec='seconds')))
        print("#     tmax={}: \"{}\"".format(tmax, dt.datetime.fromtimestamp(
            tmax/1000, gettz()).isoformat(timespec='seconds')))
        print("#{:<{wt}s}  {:>{wv}s}  {:<{wu}s}  {:<s}"
              .format(" \'timestamp (ms)\'", "\'value\'",
                      "\'unit\'", "\'channel CCS path\'",
                      wt=17, wv=9, wu=6))
        # loop over all channels
        # for chid in chanspec
        for chid in chanspec:
            path = chanspec[chid]['path']
            unitstr = chanspec[chid]['units']
            for (tstamp, value) in chandata[chid]:
                try:
                    print(
                        "{:<{wt}d} {:>{wv}g} {:>{wu}s} {:<s}".format(
                        int(tstamp), float(value), unitstr, path,
                        wt=18, wv='9.3', wu=6))
                except IOError:
                    # 'Broken pipe' IOError when stdout is closed
                    pass

    # print some statistics for each channel
    #
    if optlist.stats:
        # print a header for the stats
        #
        print("# CCS trending stats at {}".format(
            dt.datetime.now(gettz()).isoformat(timespec='seconds')))
        print("# Data for {} total seconds from {} intervals".format(
            inttot, intcnt), end="")
        print(" over {} (h:m:s) from:".format(
            dt.timedelta(seconds=(tmax/1000 - tmin/1000))))
        print("#     tmin=\"{}\"".format(
            dt.datetime.fromtimestamp(tmin/1000,
                                      gettz()).isoformat(timespec='seconds')))
        print("#     tmax=\"{}\"".format(
            dt.datetime.fromtimestamp(tmax/1000,
                                      gettz()).isoformat(timespec='seconds')))
        print("# {:>5s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s}  ".format(
            'nelem', 'mean', 'median', 'stddev',
            'rmean', 'rmedian', 'rstddev'), end="")
        print(" {:<{wt}s} {:>{wu}s}".format('path', 'units', wt=40, wu=6))

        for chid in chanspec:
            path = chanspec[chid]['path']
            unitstr = chanspec[chid]['units']
            tstamp = chandata[chid]['tstamp']
            # nelem = tstamp[(tstamp > tmin) & (tstamp < tmax)].size
            #         if ival[0] < tstamp < ival[1]:
            nelem = tstamp.size
            if nelem > 0:
                # y = chandata[chid]['value'][(tstamp > tmin) & (tstamp < tmax)]
                y = chandata[chid]['value']
                avg = np.mean(y)
                med = np.median(y)
                std = np.std(y)
                rmean, rmedian, rstd = stats.sigma_clipped_stats(y)
            else:
                avg = med = std = rmean = rmedian = rstd = 0
            try:
                print("{:>7g} {:>9.3g} {:>9.3g} {:>9.2g} ".format(
                    nelem, avg, med, std), end="")
                print("{:>9.3g} {:>9.3g} {:>9.2g}   ".format(
                    rmean, rmedian, rstd), end="")
                print("{:<{wt}s} {:>{wu}s}".format(
                    path, unitstr, wt=40, wu=6))
            except IOError:
                # 'Broken pipe' IOError when stdout is closed
                pass

    # Plotting: initial rudimentary ability...
    #
    if optlist.plot:
        # make one or more plots of the time series data
        # if optlist.overlaytime:
        # if optlist.overlayunits:
        total_axes = len(chanspec)
        # if optlist.overlaytime and optlist.overlayunits:
        fig, axes = plt.subplots(1, 1)

        ax = axes
        for chid in chandata:
            tstamp = chandata[chid]['tstamp'] / 1000.0  # to unix time
            # x = chandata[chid]['tstamp'][(tstamp > tmin) & (tstamp < tmax)]
            # y = chandata[chid]['value'][(tstamp > tmin) & (tstamp < tmax)]
            x = Time(tstamp, format='unix')
            logging.debug('x[0]= %f x[20]= %f',x[0].plot_date,x[20].plot_date)
            # x = chandata[chid]['tstamp']
            y = chandata[chid]['value']
            ax.plot_date(x.plot_date, y, '.-', label="{}".format(chanspec[chid]['path']))

    plt.show()

    # end of main()
    exit(0)


def query_rest_server(ts1, ts2, data_url, idstr):
    """get xml from restful interface for the requested channels
       with a single request
       inputs:
           t1, t2 are start/stop time in ms
           data_url is url w/out args
           idstr is string with all channel ids
       output: raw xml response from the service is returned
    """
    s = requests.Session()
    options = {'t1': int(ts1), 't2': int(ts2), 'flavor': 'raw', 'n': 1}
    uri = "{}/data/?id={}".format(data_url, idstr)
    t_start = time.time()
    try:
        resp = s.get(uri, params=options)
    except requests.ConnectionError as e:
        logging.error('ConnectionError: %s', e)
        logging.error('check status of ssh tunnel to trending server')
    if resp.status_code != 200:
        logging.error('invalid response %s from Trending Server',
                      resp.status_code)
        return None
    # logging.debug('URL=%s', resp.url)
    logging.debug('channels: %s', re.sub(r"(id=)?([0-9]+)&*", r"\2 ", idstr))
    logging.debug('dt=%.3f seconds', (time.time() - t_start))
    s.close()
    return resp.content


def get_unique_time_intervals(optlist):
    """
    Input: Command line options defining a set of intervals
           using pairs or start/stop with duration,
           empty input will return 1 defaut interval

    The set of intervals as ordered pairs in seconds are processed
    to merge overlapping periods yielding distinct intervals.

    Output: A ordered list of non-overlapping periods are returned
            as [[t00,t01], [t10,t11],...,[tn0,tn1]]
    """
    intervals = []
    duration = optlist.duration
    if optlist.interval:
        for interval in optlist.interval:
            (t1, t2) = tu.get_time_interval(interval[0], interval[1])
            intervals.append([t1, t2])
    elif optlist.start:
        for start in optlist.start:
            (t1, t2) = tu.get_time_interval(start, None, duration)
            intervals.append([t1, t2])
    elif optlist.stop:
        for stop in optlist.stop:
            (t1, t2) = tu.get_time_interval(None, stop, duration)
            intervals.append([t1, t2])
    else:
        (t1, t2) = tu.get_time_interval(None, None, duration)
        intervals.append([t1, t2])

    for interval in intervals:
        interval[0] = interval[0]
        interval[1] = interval[1]
        if interval[0] is None or interval[1] is None:
            logging.error("Date assignment failed")
            return None

    if optlist.debug:
        i = 0
        for interval in intervals:
            logging.debug(
                'time interval[%d] (before merge): %d -- %d (%d sec)',
                i, interval[0], interval[1],
                (interval[1]-interval[0])/1000)
            i += 1

    # merge overlaps to generate list of distinct intervals
    intervals.sort()  # sorts so that intervals[i][0] <= intervals[i+1][0]
    i = 1
    while i < len(intervals):  # loop over pairs of intervals
        if intervals[i-1][1] >= intervals[i][0]:
            intervals[i][0] = intervals[i-1][0]  # move left edge down
            if intervals[i-1][1] > intervals[i][1]:
                intervals[i][1] = intervals[i-1][1]  # move right edge up
            del intervals[i-1]  # delete the 1st of the pair
        else:
            i += 1  # no overlap so move to next pair

    if optlist.debug:
        i = 0
        for interval in intervals:
            logging.debug(
                'time interval[%d] (after merge): %d -- %d (%d sec)',
                i, interval[0], interval[1],
                (interval[1]-interval[0])/1000)
            i += 1

    return intervals


if __name__ == '__main__':
    main()
