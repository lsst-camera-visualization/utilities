#!/usr/bin/env python
# code to illustrate parsing of XML files
# importing the required modules
import csv
import requests
import socket
import re
import argparse
import logging
import collections
import os.path, time
from stat import *
from lxml import etree
import dateutil.parser as dp
import sys
if sys.version_info[0] < 3:
        raise Exception("Must be using Python 3")

def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="Display Trending data CCS localdb")
    parser.add_argument("channelFile", nargs='+',
                        help="file specifying channels to obtain")
    parser.add_argument("--debug", action='store_true',
                        help="print additional debugging info")
    parser.add_argument("--date", metavar='datestr',
                        help="where almost any unambiguous format is okay")
    #- note iso-8601 is from unix cmd "date --iso-8601=sconds"
    return parser.parse_args()

def main():
    """main logic"""

    optlist = parse_args()

    if optlist.debug:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.DEBUG)
        logging.debug("printing extra debugging...")
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.INFO)

    #- get time spec for the trending data pull
    if optlist.date:  #- parse the date string almost any format
        logging.debug("datestr={}".format(optlist.date))
        pt = dp.parse(optlist.date)
        if not pt.tzinfo:               #- no time zone in string
            t2 = int(time.mktime(pt.timetuple()))
        elif pt.tzinfo and pt.tzname():  #- tz name defined
            t2 = int(time.mktime(pt.timetuple()))
        elif pt.tzinfo and not pt.tzname():  #- tz name not defined
            t2 = int(time.mktime(pt.timetuple())) - 3600*time.daylight
        else:
            logging.error("could not parse the datestring...exiting")
            exit(1)
        logging.debug("t2 as localtime: {}".format(
            time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t2))))
    else:
        t2 = int(time.time())  #- right now
        logging.debug("t2 as localtime: {}".format(
            time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime(t2))))

    #- seconds to milliseconds
    t2 = t2 * 1000
    t1 = t2 - 20*1000

    data_url = "http://{}:8080/rest/data/dataserver".format(getTrendingServer())
    try:
        resp = requests.head("{}/listchannels".format(data_url))
    except requests.ConnectionError as e:
        logging.error("ConnectionError: {}".format(e))
        logging.error("check status of ssh tunnel to trending server")
        exit(1)

    #- loop over input files containing 1 channel/line
    #
    oflds = [] #- list to hold channel information
    records = dict() #- dict to hold output data
    for cfile in optlist.channelFile:
        logging.debug("channelFile= {}".format(cfile))
        try:
            cf = open(cfile, mode='r')
        except OSError as e:
            logging.error("open() failed", e)
            exit(1)

        #- populate oflds[] with the channels we want
        for line in cf:
            if re.match(r"^#", line):
                continue
            flds = line.strip().split()
            if len(flds) != 3:
                continue
            if int(flds[0]) == 1:
                oflds.append([flds[1], flds[2]])

        #- loop over channels and get xml from restful interface
        #- to be parallelized later, sequential for now
        #
        s = requests.Session()
        payload = {'t1': t1, 't2': t2, 'flavor' : 'raw', 'n' : 1}
        for ch in oflds:
            #query = "data/{}?t1={}&t2={}&flavor={}&n={}".format(
            #    ch[1], t1, t2, "raw", 1)
            t_start = time.time()
            uri = "{}/data/{}".format(data_url, ch[1])
            try:
                resp = s.get(uri, params=payload)
            except requests.ConnectionError as e:
                logging.error("ConnectionError: {}".format(e))
                logging.error("channel= {}:{}".format(ch[0],ch[1]))
                logging.error("check status of ssh tunnel to trending server")
            if resp.status_code != 200:
                logging.error("error: invalid response {} from Trending Server".format(resp.status_code))
                exit(1)
            #print resp.content
            ch.append(resp.content)
            logging.debug("#---->{}: dt={:>.3f}".format(ch[0],time.time() - t_start))
        logging.debug("t2={}".format(t2))
        s.close()

        print("# CCS trending report at {}".format(
            time.strftime("%a, %d %b %Y %H:%M:%S",
                          time.localtime(t2/1000))))
        #---- header line
        (chanw ,valuew ,alarmw ,unitw ,statew ,toffw) = (24, 9, 5, 8, 3, 4)
        hdstr1 = "{:<{w1}s} {:>{w2}s}  {:>{w3}s}".format(
                "# Channel", "Value", "Alo", w1=chanw, w2=valuew, w3=alarmw)
        hdstr2 = " {:<{w1}s} {:{w2}s} {:{w3}s}  {:>{w4}s}".format(
                "Ahi", "Units", "St", "Toff", w1=alarmw,
                w2=unitw, w3=statew, w4=toffw)
        print("{} {}".format(hdstr1, hdstr2))

        #- loop over channels and get xml from restful interface
        pcols = dict() #- columns to print out for channel
        for ch in oflds:
            logging.debug("processing xml for {}".format(ch[0]))
            #print "ch[0]={}".format(ch[0]) #- channel string
            #print "ch[1]={}".format(ch[1]) #- channel index in TDB
            #print "ch[2]={}".format(ch[2]) #- xml in string
            #- loop over channelmetadatavalue elements
            #- build columns to print as ordered dictionary
            #- columns are: (numbers formated by format attribute)
            #---- channel: ch[1]
            #---- chanval: (last instance of trendingdata(value))
            #---- alarmLow:warnLow
            #---- warnHigh:alarmHigh
            #---- state: NOM|WAR|ALA
            #---- tstamp: convert to seconds UTC, later TAI
            #---- description
            pcols.clear()  # empty out the dictionary and start over
            #- channel path goes in first
            pcols['channel'] = ch[0]
            root = etree.fromstring(ch[2])
            #- load up meta data name:value pairs
            for mdval in root.iter('channelmetadatavalue'):
                if len(mdval.keys()) > 0:
                    #- check if pcols[name] exists, if so, check tstart and only
                    #- update if this instance is newer
                    mdname = mdval.attrib.get('name') #- key
                    mdvalue = mdval.attrib.get('value') #- value
                    mdstart = mdval.attrib.get('tstart')
                    mdstop = mdval.attrib.get('tstop')
                    if mdname not in pcols:             #- new entry
                        pcols[mdname] = [mdvalue, mdstart]
                    else:                               #- update only if newer
                        if mdstart > pcols[mdname][1]:
                            pcols[mdname] = [mdvalue, mdstart]
            #- iterate over trending data and take the latest value
            for tdval in root.iter('trendingdata'):
                dataval = tdval.find('datavalue')
                if dataval is not None: value = dataval.attrib.get('value')
                else: continue
                axisval = tdval.find('axisvalue')
                if axisval is not None: tstamp = axisval.attrib.get('value')
                else: continue
                if 'chanval' not in pcols:              #- new entry
                    pcols['chanval'] = [value, tstamp]
                else:
                    if tstamp > pcols['chanval'][1]:
                        pcols['chanval'] = [value, tstamp]

            #- print stuff out now
            #print "{}".format(pcols)
            if 'chanval' not in pcols: #- no data, skip
                logging.debug("no data found for {}".format(pcols['channel']))
                continue
            else:
                chvalstr = "{:>{width}.3g} ".format(float(pcols['chanval'][0]),width=valuew)
            if len(pcols['channel']) > chanw:
                chanstr = " ...{:<{width}s}".format(pcols['channel'][-(chanw-4):],width=(chanw-4))
            else:
                chanstr = "{:<{width}s}".format(pcols['channel'],width=chanw)
            if 'alarmLow' in pcols:
                alarmLowstr = "{:>{width}g}".format(float(pcols['alarmLow'][0]),width=alarmw)
            else:
                alarmLowstr = "{:>{width}s}".format("",width=alarmw)
            if 'alarmHigh' in pcols:
                alarmHighstr = "{:<{width}g}".format(float(pcols['alarmHigh'][0]),width=alarmw)
            else:
                alarmHighstr = "{:<{width}s}".format("",width=alarmw)
            if 'units' in pcols:
                unitstr = "{:<{width}s}".format(pcols['units'][0],width=unitw)
            else:
                unitstr = "{:<{width}s}".format("",width=unitw)
            if 'state' in pcols:
                statestr = "{:<{width}.{width}s}".format(pcols['state'][0],width=statew)
            else:
                statestr = "{:<{width}.{width}s}:".format("",width=statew)
            tstampstr = "{:>{width}d}s".format(int((int(pcols['chanval'][1])-t2)/1000),width=toffw)
            print("{0} {1} {2}  {3} {4} {5} {6}".format(
                        chanstr, chvalstr, alarmLowstr,
                        alarmHighstr, unitstr, statestr, tstampstr))

#<?xml version="1.0" encoding="UTF-8" standalone="yes"?><data><trendingresult><channelmetadata/></trendingresult></data>
#- need a data_url and a query and then pull data_url/query
#- our structure for channels list
#0: datachannels [-]
#  1: datachannel [-]
#    2: path [-]
#      3: pathelement [-]
#    2: id [-]
#    2: metadata [name, value]
# we want to specify the channels as:
#   subsystem:/pathelem[0]/pathelem[1]/.../pathelem[ne] id#
# then will come up with a way to pull values of those items
#------------------------------------------------------------

def printChannelStructure(xmlcontent):
    """
    """
    #
    xml_root = etree.fromstring(xmlcontent)
    raw_tree = etree.ElementTree(xml_root)
    nice_tree = collections.OrderedDict()
    for tag in xml_root.iter():
        path = re.sub('\[[0-9]+\]', '', raw_tree.getpath(tag))
        if path not in nice_tree:
            nice_tree[path] = []
            if len(tag.keys()) > 0:
                nice_tree[path].extend(
                    attrib for attrib in tag.keys() if attrib not in nice_tree[path])
                for path, attribs in nice_tree.items():
                    indent = int(path.count('/') - 1)
                    print('{0}{1}: {2} [{3}]'.format('    ' * indent,
                                                     indent, path.split('/')[-1],
                                                     ', '.join(attribs) if len(attribs) > 0 else '-'))

def printChannelContent(xmlcontent, ssnames):
    """Walk the tree, find items with subsystem and trending matches
    print out the path and trending-ID
    """
    root = etree.fromstring(xmlcontent)
    tree = etree.ElementTree(root)
    for ssname in ssnames:
        for dchan in root.iterfind('datachannel'):
            channel_match = 0
            for md in dchan.iterfind('metadata'):
                if md.attrib.get('name') == 'subsystem':
                    if md.attrib.get('value') == ssname:
                        channel_match += 1
                if md.attrib.get('name') == 'type':
                    if md.attrib.get('value') == 'trending':
                        channel_match += 1

            if  channel_match != 2:
                continue
            parr = []
            for pp in dchan.iterfind('path'):
                for pelem in pp.iterfind('pathelement'):
                    if not len(parr):
                        parr.append("{}".format(pelem.text))
                    else:
                        parr.append("/{}".format(pelem.text))

            for eid in dchan.iterfind('id'):
                parr.append("  {}".format(eid.text))
            print("1 {}".format(''.join(parr)))

#0: datachannels [-]
#  1: datachannel [-]
#    2: path [-]
#      3: pathelement [-]
#    2: id [-]
#    2: metadata [name, value]

def getAllChannels():
    """pulls all the channels
    """
    trending_server = getTrendingServer()
    listpath="8080/rest/data/dataserver/listchannels"
    url = "http://{}:{}".format(trending_server, listpath)

    # creating HTTP response object from given url
    try:
        resp = requests.head(url)
    except requests.ConnectionError as e:
        logging.error("ConnectionError: {}".format(e))
        logging.error("check status of ssh tunnel to trending server")
        return
    if resp.status_code != 200:
        logging.error("error: invalid response {} from Trending Server".format(resp.status_code))
        return
    resp = requests.get(url)
    return resp.content

def getTrendingServer():
    """trending server is lsst-mcm or localhost via
       ssh tunnel if not on slac.stanford.edu domain
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    res = re.search(r"134\.79\.[0-9]*\.[0-9]*", s.getsockname()[0])
    if res:
        trending_server = "lsst-mcm.slac.stanford.edu"
        logging.debug("on SLAC network: trending server is lsst-mcm")
    else:
        trending_server = "localhost"
        logging.debug("off SLAC network: trending server is localhost")
    return trending_server

def updateTrendingChannelsXML():
    """maintain local cache of trending channels in xml file with
       with age limit triggering a refresh (1 day)
    """

    logging.debug("updateTrendingChannelsXML()...")
    cachedir = "{}/.trender".format(os.environ.get('HOME'))
    channelFile = "{}/.trender/listchannels.xml".format(os.environ.get('HOME'))
    #- check channelFile exists, get mtime, update if need be
    if not os.path.exists(cachedir): #- make cachdir if not exist
        os.mkdirs(cachedir)
    if not os.path.isdir(cachedir):    #- is not a directory
        logging.error("{} is not a directory, exiting...".format(cachedir))
        exit(1)
    if os.path.exists(channelFile): #- file exists or not
        statinfo = os.stat(channelFile)
        mode = statinfo.st_mode
        if not S_IWUSR&mode: #- not writeable
            os.chmod(channelFile, mode|S_IWUSR)
        delta =  time.time() - statinfo.st_mtime
        logging.debug("found existing cache age: {}s".format(delta))
        if delta < 8:
        #if delta < 86400:
            return channelFile
        logging.info("updating cached channelFile...")
        logging.debug("initial file: age={} (s)".format(deltamtime(channelFile)))
    xmlstring = getAllChannels()
    if xmlstring:
        chfile = open(channelFile, mode='w')
        chfile.write(xmlstring)
        chfile.close()
    else:
        print("failed to update cached channelFile, continuing")
    logging.debug("final file: age={} (s)".format(deltamtime(channelFile)))
    return channelFile

def deltamtime(pathname):
        return (time.time() - os.stat(pathname)[ST_MTIME])

if __name__ ==  '__main__':
    main()


