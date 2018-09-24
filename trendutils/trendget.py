#!/usr/bin/env python
""" trending data app: gets specified channels for requested time period
"""
import re
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
    parser.add_argument("--startdate", metavar='datestr',
                        help="Almost any unambiguous format work")
    parser.add_argument("--stopdate", metavar='datestr',
                        help="Almost any unambiguous format work")
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

    #  set (t1,t2) time interval for trending query
    (t1, t2) = tu.get_time_interval(optlist.startdate, optlist.stopdate)
    if t1 is None or t2 is None:
        logging.error("Date assignment failed...exiting")
        exit(1)
    logging.debug('time interval: %d -- %d (%d sec)', t1, t2, (t2-t1)/1000)

    trending_server = tu.get_trending_server()
    if trending_server:
        data_url = "http://{}:8080/rest/data/dataserver".format(
            trending_server)
    else:
        logging.error('failed to contact trending server')
        exit(1)

    # loop over input files containing 1 channel/line
    #
    oflds = []  # list to hold channel information
    for cfile in optlist.channel_file:
        logging.debug('channel_file= %s', cfile)
        try:
            cf = open(cfile, mode='r')
        except OSError as e:
            logging.error('open(%s) failed: %s', cfile, e)
            exit(1)

        # populate oflds[] with the channels we want
        for line in cf:
            if re.match(r"^\s*#", line):
                continue
            if re.match(r"^$", line):
                continue
            flds = line.strip().split()
            if len(flds) != 3:
                logging.warn("bad input file line: %s", line)
                continue
            if int(flds[0]) == 1:
                oflds.append([flds[1], flds[2]])

    # loop over channels and get xml from restful interface
    # to be parallelized later, sequential for now
    #
    s = requests.Session()
    payload = {'t1': int(t1), 't2': int(t2), 'flavor': 'raw', 'n': 1}
    for ch in oflds:
        # query = "data/{}?t1={}&t2={}&flavor={}&n={}".format(
        #    ch[1], int(t1), int(t2), "raw", 1)
        t_start = time.time()
        uri = "{}/data/?id={}".format(data_url, ch[1])
        try:
            resp = s.get(uri, params=payload)
        except requests.ConnectionError as e:
            logging.error('ConnectionError: %s', e)
            logging.error('channel= %s:%s', ch[0], ch[1])
            logging.error('check status of ssh tunnel to trending server')
        if resp.status_code != 200:
            logging.error('invalid response %s from Trending Server',
                          resp.status_code)
            exit(1)
        ch.append(resp.content)
        logging.debug('#---->%s: dt=%.3f', ch[0], (time.time() - t_start))
    s.close()

    print("# CCS trending dump at {}".format(
        dt.datetime.now(gettz()).isoformat(timespec='seconds')))
    print("# Trending data from:")
    print("#          t1=\"{}\"".format(
        time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t1/1000))))
    print("#          t2=\"{}\"".format(
        time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t2/1000))))
    print("#{:<{w1}s} {:>{w2}s} {:<{w3}s}".format(" \'timestamp (ms)\'",
                                                  "\'value\'",
                                                  "\'channel description\'",
                                                  w1=15, w2=9, w3=40))

    # loop over channels and get xml from restful interface
    pcols = dict()  # local dict for processing
    for ch in oflds:
        logging.debug('processing xml for %s', ch[0])
        # print "ch[0]={}".format(ch[0]) # channel string
        # print "ch[1]={}".format(ch[1]) # channel index in TDB
        # print "ch[2]={}".format(ch[2]) # xml in string
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
                unitstr = ""
            chanw = 40
            if len(pcols['channel']) > chanw:
                chanstr = "\' ...{:<s} ({})\'".format(
                    pcols['channel'][-(chanw-4):], unitstr)
            else:
                chanstr = "\'{:<s} ({})\'".format(
                    pcols['channel'], unitstr)
            print("{:<d} {:>9.3g} {:<s}".format(
                int(tstamp), float(value), chanstr))


if __name__ == '__main__':
    main()
