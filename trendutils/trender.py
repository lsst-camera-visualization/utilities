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
import math
import datetime as dt
import requests
from lxml import etree
from dateutil.tz import gettz
import numpy as np
from astropy import stats
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
import trendutils as tu

# plt.rcParams['text.usetex'] = True
# plt.rcParams['figure.autolayout'] = True

if sys.version_info[0] < 3 or sys.version_info[1] < 7:
    raise Exception("Must be using Python >=3.7")


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Fetch and present trending data from camera database
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
                        help="ListFile|Regex specifying channels to report")
    parser.add_argument("--input_file", nargs='+',
                        help="XML file with trending data, =>no db query")
    #
    # Output options for text based outputs
    #
    ogroup = parser.add_mutually_exclusive_group()  # all to stdout
    ogroup.add_argument("--xml", action='store_true',
                        help="Print formatted xml from trending to stdout")
    ogroup.add_argument("--text", action='store_true',
                        help="Print (timestamp, value, path) colum text")
    ogroup.add_argument("--stats", action='store_true',
                        help="Print statistics for each channel")
    ogroup.add_argument("--match", action='store_true',
                        help="print list of matching channels and exit")
    parser.add_argument("--rstats", action='store_true',
                        help="Print additional robust stats per channel")
    #
    # Plotting specifications
    #
    parser.add_argument("--plot", action='store_true',
                        help="produce plots for each channel vs. time")
    parser.add_argument("--overlay", action='store_true',
                        help="all channels on a single plot, eg. no subplots")
    parser.add_argument("--overlayunits", action='store_true',
                        help="channels grouped by units on same (sub)plot")
    #
    pgroup = parser.add_mutually_exclusive_group()  # all to stdout
    pgroup.add_argument("--overlaytime", action='store_true',
                        help="time axis is union of all time intervals")
    pgroup.add_argument("--overlaystart", action='store_true',
                        help="overlay channel data vs time from tstart")
    pgroup.add_argument("--overlaystop", action='store_true',
                        help="overlay channel data vs time until tstop")
    #
    parser.add_argument("--nolegends", action='store_true',
                        help="Don't place any legends")
    parser.add_argument("--fmt", nargs=1, metavar='format_str',
                        help="Matplotlib format string (eg. \'o-\')")
    parser.add_argument("--title", nargs='?', metavar='Title (or blank)',
                        const='auto', help="specify Title String")
    #
    # Time interval specifications
    #
    dgroup = parser.add_mutually_exclusive_group()
    dgroup.add_argument("--interval", metavar=('tstart', 'tstop'),
                        nargs=2, action='append',
                        help="Pair of date/time specifiers")
    dgroup.add_argument("--start", metavar='tstart',
                        help="Date/time specifier(s)", nargs='+')
    dgroup.add_argument("--stop", metavar='tstop',
                        help="Date/time specifier(s)", nargs='+')
    parser.add_argument("--duration", metavar='seconds', default=None,
                        help="duration [s]|(*s,*m,*h,*d) start/stop spec")
    parser.add_argument("--timebins", nargs='?', const=0, type=int,
                        metavar='nBins (blank for autosize)', default=None,
                        help="plot time avg'd bins (esp for long durations)")
    parser.add_argument("--debug", action='store_true',
                        help="Print additional debugging info")
    #
    parser.add_argument("--forceupdate", action='store_true',
                        help="Force update on cached channel file")

    #
    #
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
    intmax = int(max([t[1] - t[0] for t in intervals])/1000)
    tmin = intervals[0][0]
    tmax = intervals[-1][1]
    tz_trending = tu.tz_trending

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
    # 2 exclusive choices:
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
    channel_file = None
    regexes = []
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
            if oflds:  # only set if some lines were good
                channel_source_type = 'file'

        if channel_source_type == 'pattern' or not channel_source_type:
            logging.debug('eval pattern for matching channels...')
            if not channel_file:
                if optlist.forceupdate:
                    channel_file = tu.update_trending_channels_xml()
                else:
                    channel_file = tu.update_trending_channels_xml(
                        tmin/1000, tmax/1000)

            if not chid_dict:
                chid_dict = get_channel_dict(channel_file)

            # add csource as a new parameter in oflds[csource][chid]?
            # search the entire catalog for each pattern, this may be slow
            cpat = re.compile(csource)
            for chid in chid_dict:
                if cpat.search(chid_dict[chid]):
                    oflds[chid] = chid_dict[chid]
                    #  oflds[csource][chid] = chid_dict[chid]

            if oflds:
                channel_source_type = 'pattern'
                regexes.append(csource)

    del chid_dict
    # end of loop over input sources to define the channels for query/output
    # now have info needed to query the CCS trending db

    if optlist.match:
        print('Found matching channels:')
        for chid in oflds:
            print("   id: {}  path: {}".format(chid, oflds[chid]))
        exit(0)

    if optlist.debug:
        logging.debug('Found matching channels:')
        for chid in oflds:
            logging.debug('id= %5d  path= %s', int(chid), oflds[chid])

    #  Get the trending data either from local saved files or via
    #  trending db queries to the rest service
    if optlist.input_file:
        # get input from these files rather than trending service
        # an issue is that the input file need not have the
        # same set of channels or time intervals as requested on command line.
        # We expect to apply the intervals to the output only by taking
        # only times in the intersection.  Up to user to use valid file.
        responses = []
        parser = etree.XMLParser(remove_blank_text=True)
        for ifile in optlist.input_file:
            logging.debug('using %s for input', ifile)
            logging.debug('test for well-formed xml...')
            try:
                tree = etree.parse(ifile, parser)
            except etree.ParseError as e:
                logging.debug('parsing %s failed: %s', ifile, e)
                exit(1)
            except etree.XMLSyntaxError as e:
                logging.debug('parsing %s failed: %s', ifile, e)
                exit(1)
            else:
                logging.debug('successfully parsed %s', ifile)

            logging.debug('appending to responses...')
            responses.append(etree.tostring(
                tree.getroot(), encoding="UTF-8",
                xml_declaration=True, pretty_print=False))

            logging.debug('deleting the etree')
            del tree

    else:
        # query the rest server query and place responses into a list
        # join the ids requested as "id0&id=id1&id=id2..." for query
        idstr = '&id='.join(id for id in oflds)
        responses = []
        for ival in intervals:  # only one interval per query allowed
            if optlist.timebins == 0:  # autosize it
                logging.debug('timebins=0')
                if int((ival[1] - ival[0])/1000) < 300:  # just get raw data
                    logging.debug('ival[1]= %d, ival[0]= %d', ival[1], ival[0])
                    optlist.timebins = None
                else:
                    logging.debug('ival[1]= %d, ival[0]= %d', ival[1], ival[0])
                    optlist.timebins = int(((ival[1] - ival[0])/1000.0) / 60.0)
            logging.debug('timebins= %d', optlist.timebins)
            res = query_rest_server(ival[0], ival[1],
                                    data_url, idstr, optlist.timebins)
            responses.append(res)
    # Now have the data from trending service

    # Output to stdout a well formed xml tree aggregating the xml received
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
        try:
            os.write(1, b'</datas>')
        except OSError:
            # 'Broken pipe' OSError when stdout is closed
            pass

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
            if path is None or oflds[chid] != path:
                logging.warning(
                    'inputpath(id=%s): %s != %s (xmlpath), using %s',
                    chid, oflds[chid], path, oflds[chid])
                path = oflds[chid]
                # continue
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
            # metadata:
            # parse all but only using units for now
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
                for ival in intervals:  # slow, but no other way?
                    if ival[0] < int(tstamp) < ival[1]:
                        chandata[chid].append((tstamp, tvalue))
                        break

    # Done ranslating the xml responses into internal arrays etc.
    # Delete all the raw xml responses
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
    trimids = []
    for chid in chanspec:
        path = chanspec[chid]['path']
        logging.debug('id=%s  path=%s', chid, path)
        for mdname in chanmd[chid]:
            # pick out and process the md's we want
            if mdname == 'units' and chanspec[chid]['units'] == 'none':
                chanspec[chid]['units'] = chanmd[chid][mdname][-1][1]

        logging.debug('    units=%s', chanspec[chid]['units'])
        # sort and remove duplicates from chandata[chid] where:
        # chandata[chid] = [(t0, v0), (t1, v1), ...., (tn, vn)]
        # and convert to np array
        tmparr = np.array(chandata[chid], dtype=chandt)
        chandata[chid] = np.unique(tmparr)
        logging.debug('    chandata: %d uniq/sorted values from %d entries',
                      np.size(chandata[chid]), np.size(tmparr))
        del tmparr
        if np.size(chandata[chid]) == 0:  # append chid to trimid list
            logging.debug('%s has no data', chanspec[chid]['path'])
            # arrange to trim empty data
            trimids.append(chid)

    for chid in trimids:
        del chandata[chid]
        del chanmd[chid]
        del chanspec[chid]

    # print to stdout a text dump of the data, in time order per channel
    #
    if optlist.text:
        # print a header for the text
        #
        print("#")
        print("# {}".format(optlist.title))
        print("#")
        print("# CCS trending dump at {}".format(
            dt.datetime.now(gettz()).isoformat(timespec='seconds')))
        print("# Data for {} total seconds from {} intervals".format(
            inttot, intcnt), end="")
        print(" over {} (h:m:s) from:".format(
            dt.timedelta(seconds=(tmax/1000 - tmin/1000))))
        print("#     tmin={}: \"{}\"".format(tmin, dt.datetime.fromtimestamp(
            tmin/1000, gettz(tz_trending)).isoformat(timespec='seconds')))
        print("#     tmax={}: \"{}\"".format(tmax, dt.datetime.fromtimestamp(
            tmax/1000, gettz(tz_trending)).isoformat(timespec='seconds')))
        print("#{:<{wt}s}  {:>{wv}s}  {:<{wu}s}  {:<s}"
              .format(" \'timestamp (ms)\'", "\'value\'",
                      "\'unit\'", "\'channel CCS path\'",
                      wt=17, wv=9, wu=6))
        # loop over all channels
        # for chid in chanspec
        for chid in chanspec:
            path = chanspec[chid]['path']
            unitstr = chanspec[chid]['units']
            if np.size(chandata[chid]) == 0:
                continue
            for (tstamp, value) in chandata[chid]:
                try:
                    print(
                        "{:<{wt}d} {:>{wv}g} {:>{wu}s} {:<s}".format(
                            int(tstamp), float(value), unitstr, path,
                            wt=18, wv='9.4', wu=6))
                except IOError:
                    # 'Broken pipe' IOError when stdout is closed
                    pass

    # print some statistics for each channel
    #
    if optlist.stats:
        # print a header for the stats
        #
        print("#")
        print("# {}".format(optlist.title))
        print("#")
        print("# CCS trending stats at {}".format(
            dt.datetime.now(gettz()).isoformat(timespec='seconds')))
        print("# Data for {} total seconds from {} intervals".format(
            inttot, intcnt), end="")
        print(" over {} (h:m:s) from:".format(
            dt.timedelta(seconds=(tmax/1000 - tmin/1000))))
        print("#     tmin=\"{}\"".format(
            dt.datetime.fromtimestamp(
                tmin/1000, gettz(tz_trending)).isoformat(timespec='seconds')))
        print("#     tmax=\"{}\"".format(
            dt.datetime.fromtimestamp(
                tmax/1000, gettz(tz_trending)).isoformat(timespec='seconds')))
        print("# {:>4s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} ".format(
            'cnt', 'mean', 'median', 'stddev', 'min', 'max'), end="")
        if optlist.rstats:
            print("{:>8s} {:>8s} {:>8s}  ".format(
                'rmean', 'rmedian', 'rstddev'), end="")
        print(" {:<{wt}s} {:>{wu}s}".format('path', 'units', wt=40, wu=6))

        for chid in chanspec:
            path = chanspec[chid]['path']
            unitstr = chanspec[chid]['units']
            tstamp = chandata[chid]['tstamp']
            nelem = tstamp.size
            if nelem > 0:
                y = chandata[chid]['value']
                avg = np.mean(y)
                med = np.median(y)
                std = np.std(y)
                npmin = np.min(y)
                npmax = np.max(y)
                if optlist.rstats:
                    rmean, rmedian, rstd = stats.sigma_clipped_stats(y)
            else:
                avg = med = std = npmin = npmax = 0
                grad = rmean = rmedian = rstd = 0
            try:
                print("{:>6g} {:>8.3g} {:>8.3g} {:>8.2g} ".format(
                    nelem, avg, med, std,), end="")
                print("{:>8.3g} {:>8.3g} ".format(npmin, npmax), end="")
                if optlist.rstats:
                    print("{:>8.3g} {:>8.3g} {:>8.2g}   ".format(
                        rmean, rmedian, rstd), end="")
                print("{:<{wt}s} {:>{wu}s}".format(
                    path, unitstr, wt=40, wu=6))
            except IOError:
                # 'Broken pipe' IOError when stdout is closed
                pass

    # Plotting:
    #
    if optlist.plot:
        # make one or more plots of the time series data
        # default is plot per channel per interval
        # option to combine intervals and channels by units
        # and to overlay all

        # figure out how many distinct plots and windows to make
        # subplots layout and shape are determined
        # nax will store the number of actual plots
        # the nxm array of plots may be larger
        #
        if optlist.overlayunits:  # axis set per unit
            unit_order = dict()
            unit_idx = 0  # counts types of units
            for chid in chanspec:
                unit = chanspec[chid]['units']
                if unit not in unit_order:
                    unit_order[unit] = unit_idx
                    unit_idx += 1
            nax = len(unit_order)
            logging.debug('unit_order=%s', unit_order)
        else:                     # nominal, one per channel
            nax = len(chanspec)

        if nax == 0:
            logging.error('no data to plot, check inputs?')
            logging.error('try running with --debug')
            exit(1)

        if not optlist.overlaytime and\
            not optlist.overlaystart and\
                not optlist.overlaystop:
            nax = nax * len(intervals)  # per interval, per channel

        if optlist.overlay:  # just one plot box (axis)
            nax = nrows = ncols = 1
        else:
            # arrange for plot array shape ~landscape
            ncols = int(math.ceil(math.sqrt(nax)))
            nrows = int(math.ceil(float(nax)/ncols))
            while nrows > 1 and float(ncols)/float(nrows) < 1.3:
                ncols += 1
                nrows = int(math.ceil(float(nax)/ncols))

        logging.debug('nax=%d', nax)
        logging.debug('nrows= %d  ncols=%d', nrows, ncols)

        if optlist.overlay:  # just one plot box
            fig, axes = plt.subplots(1, 1)
        else:
            (fig_width, fig_height) = plt.rcParams["figure.figsize"]
            fig, axes = plt.subplots(nrows, ncols, sharex=False,
                                     figsize=(fig_width*1.6, fig_height))

        # loop over data channels and plot them on correct axis
        chids = list(chanspec.keys())
        logging.debug('chids=%s', chids)
        unit = None
        for chidx in range(0, len(chids)):  # channels
            chid = chids[chidx]
            unit = chanspec[chid]['units']
            tstamp = chandata[chid]['tstamp']
            mcolor = None
            for idx in range(0, len(intervals)):
                # choose on which axis to plot
                if optlist.overlayunits:
                    axcnt = unit_order[unit]  # map unit to correct axis
                else:
                    axcnt = chidx
                if not optlist.overlaytime and\
                        not optlist.overlaystart and\
                        not optlist.overlaystop:
                    # stride is number of intervals
                    axcnt = axcnt*len(intervals) + idx

                if optlist.overlay:
                    axcnt = 0
                ax = np.ravel(axes)[axcnt]
                # ax.grid(True)
                logging.debug('axcnt= %d  idx= %d', axcnt, idx)
                # mask the tstamps outside of the interval
                mask = ((intervals[idx][0] < tstamp) &
                        (tstamp < intervals[idx][1]))
                mask_start = intervals[idx][0] / 1000.0
                mask_stop = intervals[idx][1] / 1000.0
                x = chandata[chid]['tstamp'][mask] / 1000.0  # time in seconds
                y = chandata[chid]['value'][mask]
                # deal with point/line format
                if optlist.fmt:
                    fmt = optlist.fmt[0]
                else:
                    if optlist.timebins:
                        fmt = '|-'
                    else:
                        fmt = 'o-'
                # do the actual plotting
                if not optlist.overlaystart and not optlist.overlaystop:
                    #
                    # convert time axis to matplotlib dates sequence
                    #
                    dates = [dt.datetime.fromtimestamp(ts, gettz(tz_trending))
                             for ts in x]
                    mds = mdate.date2num(dates)
                    if idx == 0:  # label on first interval, store color
                        mlabel = "{}".format(chanspec[chid]['path'])
                        line = ax.plot_date(mds, y, fmt, label=mlabel,
                                            markersize=4.0,
                                            tz=gettz(tz_trending))
                        mcolor = line[0].get_color()
                        logging.debug("mcolor= %s", mcolor)
                    else:         # no label on later intervals, use same color
                        line = ax.plot_date(mds, y, fmt, label=None,
                                            color=mcolor, markersize=4.0,
                                            tz=gettz(tz_trending))
                    if not ax.get_ylabel():
                        ax.set_ylabel("{}".format(unit))
                    # xlabel for this axis
                    if not ax.get_xlabel():
                        xlabel_str = "{} (tstart)".format(
                            dt.datetime.fromtimestamp(
                                intervals[idx][0]/1000,
                                gettz(tz_trending)).isoformat(
                                    timespec='seconds'))
                        if optlist.timebins:
                            xlabel_str = "{} [{} timebins]".format(
                                xlabel_str, optlist.timebins)
                        logging.debug('ax.set_xlabel(%s)', xlabel_str)
                        ax.set_xlabel("{}".format(xlabel_str),
                                      position=(0., 1e6),
                                      horizontalalignment='left')
                        fig.autofmt_xdate()
                else:
                    # convert x to duration axis units (s+/-)
                    #
                    if optlist.overlaystart:
                        x = x - mask_start
                    elif optlist.overlaystop:
                        x = x - mask_stop
                    else:
                        logging.error('overlaystart/stop problem')
                        exit(1)

                    mlabel = "{}[{}]".format(chanspec[chid]['path'], idx)
                    line = ax.plot(x, y, fmt, label=mlabel, markersize=4.0)
                    mcolor = line[0].get_color()
                    logging.debug("mcolor= %s", mcolor)
                    if not ax.get_ylabel():
                        ax.set_ylabel("{}".format(unit))
                    # xlabel for this axis
                    if not ax.get_xlabel():
                        if optlist.overlaystart:
                            xstr = "tstart"
                            xid = 0
                        if optlist.overlaystop:
                            xstr = "tstop"
                            xid = 1
                        xlabel_str = "{} ({}[0])".format(
                            dt.datetime.fromtimestamp(
                                intervals[0][xid]/1000,
                                gettz(tz_trending)).isoformat(
                                    timespec='seconds'), xstr)
                        if len(intervals) > 1:
                            xlabel_str1 = "{} ({}[{}])".format(
                                dt.datetime.fromtimestamp(
                                    intervals[-1][xid]/1000,
                                    gettz(tz_trending)).isoformat(
                                        timespec='seconds'), xstr,
                                len(intervals) - 1)
                            if len(intervals) > 2:
                                xlabel_str1 = "...{}".format(xlabel_str1)
                            xlabel_str = "{}\n{}".format(xlabel_str,
                                                         xlabel_str1)
                        ax.set_xlabel("{}".format(xlabel_str,
                                                  position=(0., 1e6),
                                                  horizontalalignment='left'))

        # plot array padded with invisible boxes
        for pcnt in range(nax, nrows*ncols):
            ax = np.ravel(axes)[pcnt]
            ax.grid(False)
            ax.set_frame_on(False)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

        # make the legends for each plot in the array
        # with legend placement/size adjustments
        if not optlist.nolegends:
            for pcnt in range(0, nax):
                ax = np.ravel(axes)[pcnt]
                handles, labels = ax.get_legend_handles_labels()
                if not handles or not labels:
                    continue
                # sort the labels for easier reading
                # using https://stackoverflow.com/questions/9764298
                logging.debug('len(handles)=%d  len(labels)=%d',
                              len(handles), len(labels))
                labels, handles = (
                    list(t) for t in zip(*sorted(zip(labels, handles))))
                if labels:
                    if len(handles) < 4:  # place in the box
                        ax.legend(handles, labels, loc='best', fancybox=True,
                                  framealpha=0.5, fontsize='x-small')
                    elif len(handles) < 8:  # place to right, small
                        ax.legend(handles, labels, loc='upper left',
                                  bbox_to_anchor=(1, 1), fontsize='small')
                    elif len(handles) < 25:  # place to right, xx-small
                        ax.legend(handles, labels, loc='upper left',
                                  bbox_to_anchor=(1, 1), fontsize='x-small')
                    else:                    # to right, truncate legend
                        labels[24] = "truncated labels[{}:{}]".format(
                            24, len(handles))
                        ax.legend(handles[:25], labels[:25], loc='upper left',
                                  bbox_to_anchor=(1, 1), fontsize='x-small')

        if optlist.title:
            logging.debug("using title=%s", optlist.title)
            fig.suptitle("{}".format(optlist.title))
            fig.set_tight_layout({"pad": 2.4})
            plt.show()
        else:
            fig.set_tight_layout(True)
            plt.show()

    # end of main()
    exit(0)


def query_rest_server(ts1, ts2, data_url, idstr, nbins):
    """get xml from restful interface for the requested channels
       with a single request
       inputs:
           t1, t2 are start/stop time in ms
           data_url is url w/out args
           idstr is all channel ids as '&id='.join(id for id in oflds)
           nbins: request raw or binned data from DB
       output: raw xml response from the service is returned
    """
    s = requests.Session()
    if nbins is None:  # raw data
        options = {'t1': int(ts1), 't2': int(ts2), 'flavor': 'raw', 'n': 1}
    else:  # CCS stat data
        options = {'t1': int(ts1), 't2': int(ts2),
                   'flavor': 'stat', 'n': int(nbins)}
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
            as [[t00,t01], [t10,t11], ...,[tn0,tn1]]
    """
    intervals = []
    if optlist.interval:
        for interval in optlist.interval:
            (t1, t2) = tu.get_time_interval(interval[0], interval[1])
            intervals.append([t1, t2])
    elif optlist.start:
        for start in optlist.start:
            (t1, t2) = tu.get_time_interval(start, None, optlist.duration)
            intervals.append([t1, t2])
    elif optlist.stop:
        for stop in optlist.stop:
            (t1, t2) = tu.get_time_interval(None, stop, optlist.duration)
            intervals.append([t1, t2])
    else:
        (t1, t2) = tu.get_time_interval(None, None, optlist.duration)
        intervals.append([t1, t2])

    for interval in intervals:
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


def get_channel_dict(channel_file):
    """
    build a channel dictionary for all channels
    struction of channel_file is:
    0: datachannels [-]
      1: datachannel [-]
        2: path [-]
          3: pathelement [-]
      2: id [-]
      2: metadata [name, value]
    """
    logging.debug('building full channel dictionary from %s', channel_file)
    cdict = dict()
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
        if path and chid:  # create entry in dict
            cdict[chid] = path
    logging.debug('channel dict contains %d active channels', len(cdict))
    del tree
    return cdict


if __name__ == '__main__':
    main()
