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

        The "channel_file"(s) are best constructed using the sibling
        application "trendingChannels.py" and editing the resuling files
        to choose which channels to activate.

        The start/stop date specifiers could be nearly any variant of a
        full date/time spec as output by the unix "date <options>" command.
        One particular choice is the format from "date --iso-8601=seconds".
                               '''))
    parser.add_argument("channel_file", nargs='+',
                        help="File specifying channels to obtain")
    parser.add_argument("--debug", action='store_true',
                        help="Print additional debugging info")
    parser.add_argument("--rawxml", action='store_true',
                        help="Print raw xml from trending to stdout")

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
    return parser.parse_args()


def main():
    """main logic"""
    # get command args and options
    optlist = parse_args()

    # print("optlist={}".format(optlist))
    # exit(0)

    # logging level
    if optlist.debug:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.DEBUG)
        logging.debug("printing extra debugging...")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.INFO)

    # get list of time intervals to process
    intervals = get_intervals(optlist)
    # exit(0)

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

    # loop over input files containing 1 channel/line
    # oflds[channel_id_num] = trending_full_path
    # eg: oflds[2372] = aliveness-raft/R00.Reb2.RGL
    #
    oflds = dict()  # dict to hold channel information
    #
    # for the case where we are taking the requested channels from
    # one or more input "channel" files
    #
    for cfile in optlist.channel_file:
        logging.debug('channel_file= %s', cfile)
        try:
            cf = open(cfile, mode='r')
        except OSError as e:
            logging.error('open(%s) failed: %s', cfile, e)
            exit(1)

        # populate oflds[id] with the corresponding channel path
        for line in cf:
            if re.match(r"^\s*#", line):  # skip block comment
                continue
            if re.match(r"^\s*$", line):  # skip white space line
                continue
            sline = re.sub(r"""(#[^\'^"]*$)""", '', line)  # strip inline cmnt
            flds = [''.join(t) for t in rpat.findall(sline)]
            if len(flds) != 3:
                logging.warning('bad input line: %s', line)
                continue
            if int(flds[0]) == 1:
                oflds[flds[2]] = flds[1]  # oflds[id] = path

    # join the ids requested as "id0&id=id1&id=id2..." for query
    idstr = '&id='.join(id for id in oflds)

    # put the rest server query responses into a list
    responses = []
    for interval in intervals:
        res = query_rest_server(interval[0], interval[1], data_url, idstr)
        responses.append(res)

    # Output a well formed xml tree aggregating all the xml received
    # Mainly for debugging but useful for external processing too
    if optlist.rawxml:
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
    #     4: channelmetadatavalue [tstart, tstop, name, value]
    #     3: trendingdata [-]
    #     4: datavalue [name, value]
    #     4: axisvalue [name, value, loweredge, upperedge]
    # where [id, path] could appear multiple times
    #

    chanspec = dict()  # where keys are channel ids, values are ccs paths
    chanmd = []  # 1st index channel id, elements will be array of dicts
    for res in responses:
        root = etree.fromstring(res)
        for data in root.iter('data'):
            chid = data.attrib.get('id')
            path = data.attrib.get('path')
            # verify this element's (chid, path) matches the input list
            logging.debug('id=%s  path=%s', chid, path)
            if oflds[chid] != path:
                logging.error('inputpath(id=%s): %s != %s = xmlpath',
                              chid, oflds[chid], path)
            chanspec[chid] = (path)
            chanmd[chid] = []

    exit(0)

        pcols.clear()
        #  channel path goes in first
        pcols['channel'] = ch[0]
        root = etree.fromstring(ch[2])
        # load up meta data name:value pairs
        for mdval in root.iter('channelmetadatavalue'):
            if mdval.keys():  # empty sequence is false
                # check if pcols[name] exists, if so, check tstart and only
                # update if this instance is newer
                mdname = mdval.attrib.get('name')  # key
                mdvalue = mdval.attrib.get('value')  # value
                mdstart = mdval.attrib.get('tstart')
                mdstop = mdval.attrib.get('tstop')
                if mdname not in pcols:             # new entry
                    pcols[mdname] = [mdvalue, mdstart]
                else:                               # update only if newer
                    if mdstart > pcols[mdname][1]:
                        pcols[mdname] = [mdvalue, mdstart]
    # ----------------------------------------------------------------
#   print("# CCS trending dump at {}".format(
#       dt.datetime.now(gettz()).isoformat(timespec='seconds')))
#   print("# Trending data from:")
#   print("#          t1=\"{}\"".format(
#       time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t1/1000))))
#   print("#          t2=\"{}\"".format(
#       time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t2/1000))))
#   print("#{:<{wt}s}  {:>{wv}s}  {:<{wu}s}  {:<s}"
#         .format(" \'timestamp (ms)\'", "\'value\'",
#                 "\'unit\'", "\'channel description\'",
#                 wt=17, wv=9, wu=6))
    # resp
    # loop over channels and get xml from restful interface
    pcols = dict()  # local dict for processing
    for ch in oflds:
        logging.debug('processing xml for %s', ch[0])
        # print "ch[0]={}".format(ch[0]) # channel string (eg. CCS path)
        # print "ch[1]={}".format(ch[1]) # channel index in TDB
        # print "ch[2]={}".format(ch[2]) # xml as a string
        #  loop over channelmetadatavalue elements
        #  build columns to print as ordered dictionary
        #  columns are: (numbers formated by format attribute)
        # ---- channel: ch[1]
        # ---- chanval: (last instance of trendingdata(value))
        # ---- alarmLow:warnLow
        # ---- warnHigh:alarmHigh
        # ---- state: NOM|WAR|ALA
        # ---- tstamp: convert to seconds UTC, later TAI
        # ---- description
        pcols.clear()
        #  channel path goes in first
        pcols['channel'] = ch[0]
        root = etree.fromstring(ch[2])
        # load up meta data name:value pairs
        for mdval in root.iter('channelmetadatavalue'):
            if mdval.keys():  # empty sequence is false
                # check if pcols[name] exists, if so, check tstart and only
                # update if this instance is newer
                mdname = mdval.attrib.get('name')  # key
                mdvalue = mdval.attrib.get('value')  # value
                mdstart = mdval.attrib.get('tstart')
                mdstop = mdval.attrib.get('tstop')
                if mdname not in pcols:             # new entry
                    pcols[mdname] = [mdvalue, mdstart]
                else:                               # update only if newer
                    if mdstart > pcols[mdname][1]:
                        pcols[mdname] = [mdvalue, mdstart]
        # iterate over trending data and print results
        for tdval in root.iter('trendingdata'):
            dataval = tdval.find('datavalue')
            if dataval is not None:
                value = dataval.attrib.get('value')
            else:
                continue
            axisval = tdval.find('axisvalue')
            if axisval is not None:
                tstamp = axisval.attrib.get('value')
            else:
                continue
            if 'units' in pcols:
                unitstr = "{:<s}".format(pcols['units'][0])
            else:
                unitstr = "none"
            # chanw = 44
            # if len(pcols['channel']) > chanw:
            #     chanstr = "...{:<s}".format(pcols['channel'][-(chanw-3):])
            # else:
            #  chanstr = "{:<{width}s}".format(pcols['channel'], width=chanw)
            #
            try:
                print("{:<{wt}d} {:>{wv}g} {:>{wu}s} {:<s}".format(
                    int(tstamp), float(value), unitstr, pcols['channel'],
                    wt=18, wv='9.3', wu=6))
            except IOError:
                # 'Broken pipe' IOError when stdout is closed
                pass


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
    logging.debug('URL=%s', resp.url)
    logging.debug('#---->%s: dt=%.3f', "channels", (time.time() - t_start))
    s.close()
    return resp.content


def get_intervals(optlist):
    #  set (t1,t2) time intervals for trending query
    #  intervals will be a list of (t1, t2) tuples
    intervals = []
    duration = optlist.duration
    if optlist.interval:
        for interval in optlist.interval:
            (t1, t2) = tu.get_time_interval(interval[0], interval[1])
            intervals.append((t1, t2))
    elif optlist.start:
        for start in optlist.start:
            (t1, t2) = tu.get_time_interval(start, None, optlist.duration)
            intervals.append((t1, t2))
    elif optlist.stop:
        for stop in optlist.stop:
            (t1, t2) = tu.get_time_interval(None, stop, optlist.duration)
            intervals.append((t1, t2))
    else:
        (t1, t2) = tu.get_time_interval(None, None, optlist.duration)
        intervals.append((t1, t2))

    for interval in intervals:
        (t1, t2) = interval
        if t1 is None or t2 is None:
            logging.error("Date assignment failed")
            return None
        logging.debug('time interval: %d -- %d (%d sec)', t1, t2, (t2-t1)/1000)
    return intervals


if __name__ == '__main__':
    main()
