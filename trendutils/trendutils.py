import socket
import re
import logging
import collections
import os.path
import time
import stat
from datetime import datetime
import requests
from lxml import etree
import dateutil.parser as dp
from dateutil.tz import gettz
from dateutil.tz import tzutc

# constants? for lack of a better term (default to slac)
trendnetre = r"134\.79\.[0-9]*\.[0-9]*"
default_trending_server = "lsst-mcm.slac.stanford.edu"
tz_trending = 'America/Los_Angeles'


def deltamtime(pathname):
    """ time since last mod
    """
    return time.time() - os.stat(pathname)[stat.ST_MTIME]


def get_all_channels(maxidle=-1):
    """pulls all the channels from rest server
    """
    trending_server = get_trending_server()
    if not trending_server:
        return None
    listpathroot = "8080/rest/data/dataserver/listchannels?maxIdleSeconds="
    if maxidle:
        listpath = "{}{:d}".format(listpathroot, maxidle)
    else:
        listpath = "{}{:s}".format(listpathroot, "-1")
    url = "http://{}:{}".format(trending_server, listpath)
    logging.debug('url=%s', url)

    # creating HTTP response object from given url
    resp = requests.get(url)
    return resp.content


def parse_datestr(datestr):
    """parse date string to return time in seconds
       or now if the string is None
       The tz_trending timezone specification is used
    """
    if datestr:  # parse the date string almost any format
        logging.debug('datestr=%s', datestr)
        try:
            dt = dp.parse(datestr)
        except ValueError as e:
            logging.error('ValueError: %s', e)
            logging.error('could not parse the datestring')
            return None

        if not dt.tzinfo:  # no time zone, attach tz_trending
            dt = dt.replace(tzinfo=gettz(tz_trending))  # default

        # convert timezone to tz_trending
        dt = dt.astimezone(gettz(tz_trending))
    else:
        dt = datetime.now(tz=gettz(tz_trending))

    return (dt - datetime(1970, 1, 1, tzinfo=tzutc())).total_seconds()


def get_time_interval(startstr, stopstr, duration=None):
    """return the timeinterval boundaries (ms) from datestrings
    """
    if startstr:
        t1 = parse_datestr(startstr)
        logging.debug('t1 as trending db local time: %s',
                      datetime.fromtimestamp(t1, gettz(tz_trending)))
        t1 *= 1000

    if stopstr:
        t2 = parse_datestr(stopstr)
        logging.debug('t2 as trending db local time: %s',
                      datetime.fromtimestamp(t2, gettz(tz_trending)))
        t2 *= 1000

    if not duration:
        duration = 600

    # cases for t1, t2
    if startstr and stopstr:
        if t1 > t2:
            logging.error('starttime must be earlier than stop')
            t1 = None
            t2 = None
    elif startstr and not stopstr:  # get durations interval
        t2 = t1 + duration*1000
    elif not startstr and stopstr:  # get durations interval
        t1 = t2 - duration*1000
    elif not startstr and not stopstr:
        t2 = int(time.time()) * 1000
        t1 = t2 - duration*1000
    else:
        t1 = t2 = None

    return (t1, t2)


def get_trending_server():
    """trending server is lsst-mcm or localhost via
       ssh tunnel if not on slac.stanford.edu domain
       returns valid service or None on failure
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    res = re.search(trendnetre, s.getsockname()[0])
    s.close()
    if res:
        trending_server = default_trending_server
        logging.debug('on SLAC network: trending server is lsst-mcm')
    else:
        trending_server = "localhost"
        logging.debug('outside trending network: trending server is localhost')
    # test connection to trending server, just get head to verify
    #
    try:
        requests.head("http://{}:8080".format(trending_server))
    except requests.ConnectionError as e:
        logging.error('ConnectionError: %s', e)
        logging.error('check status connection to trending server: %s',
                      trending_server)
        return None
    return trending_server


def geturi(uri):
    """return the xml file from internet or locally
    """
    res = re.match(r"http.*", uri)
    if res:
        #  we have a url, use request to return it
        try:
            resp = requests.head(uri)
        except requests.ConnectionError as e:
            logging.error('ConnectionError: %s', e)
            return None
        if resp.status_code != 200:
            logging.error(
                'error: invalid response %s from Server', resp.status_code)
            return None
        resp = requests.get(uri)
        return resp.content

    #  it may be a file, try to open as a path
    if not os.path.exists(uri):
        logging.error('%s is not a path', uri)
        return None
    try:
        xmlfile = open(uri, mode='rb')
    except IOError as e:
        logging.error('I/O error(%s): %s', e.errno, e.strerror)
        return None
    xmlstring = xmlfile.read()
    xmlfile.close()
    return xmlstring


def print_channel_content(xmlcontent, ssnames):
    """Walk the tree, find items with subsystem and trending matches
       print out the path and trending-ID
    """
    root = etree.fromstring(xmlcontent)
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

            if channel_match != 2:
                continue
            parr = []
            for pp in dchan.iterfind('path'):
                for pelem in pp.iterfind('pathelement'):
                    if not parr:
                        parr.append("{}".format(pelem.text))
                    else:
                        parr.append("/{}".format(pelem.text))

            for eid in dchan.iterfind('id'):
                parr.append("  {}".format(eid.text))
            print("1 {}".format(''.join(parr)))


def print_channel_structure(xmlcontent):
    """ print out xml tree structure using algorithm from stackoverflow.com
    """
    #
    xml_root = etree.fromstring(xmlcontent)
    raw_tree = etree.ElementTree(xml_root)
    nice_tree = collections.OrderedDict()
    for tag in xml_root.iter():
        path = re.sub(r'\[[0-9]+\]', '', raw_tree.getpath(tag))
        if path not in nice_tree:
            nice_tree[path] = []
            if tag.keys():
                nice_tree[path].extend(
                    attrib for attrib in tag.keys()
                    if attrib not in nice_tree[path])
                for path, attribs in nice_tree.items():
                    indent = int(path.count('/') - 1)
                    print('{0}{1}: {2} [{3}]'.format(
                        '    ' * indent, indent, path.split('/')[-1],
                        ', '.join(attribs) if attribs else '-'))


def update_trending_channels_xml(maxidle=None):
    """maintain local cache of trending channels in xml file with
       with age limit triggering a refresh (1 day)
    """
    logging.debug("update_trending_channels_xml()...")
    cachedir = "{}/.trender".format(os.environ.get('HOME'))
    channel_file = "{}/.trender/listchannels.xml".format(
        os.environ.get('HOME'))
    # check channel_file exists, get mtime, update if need be
    if not os.path.exists(cachedir):   # make cachdir if not exist
        os.mkdir(cachedir)
    if not os.path.isdir(cachedir):    # is not a directory
        logging.error('%s is not a directory, exiting...', cachedir)
        exit(1)
    if os.path.exists(channel_file) or maxidle:
        statinfo = os.stat(channel_file)
        mode = statinfo.st_mode
        if not stat.S_IWUSR & mode:  # not writeable
            os.chmod(channel_file, mode | stat.S_IWUSR)
        delta = time.time() - statinfo.st_mtime
        logging.debug('found existing cache age: %s (s)', delta)
        update = True if delta > 86400 else False
        if not update:
            logging.debug('returning existing channel_file=%s',
                          channel_file)
    else:  # file does not exist
        update = True
    if update:
        logging.info('updating cached channel_file...')
        logging.debug('initial file: age=%s (s)', deltamtime(channel_file))
        xmlstring = get_all_channels(maxidle)
        if xmlstring:
            chfile = open(channel_file, mode='wb')
            chfile.write(xmlstring)
            chfile.close()
        else:
            logging.info('failed update cached channel_file, continuing')
        logging.debug('final file: age=%s (s)', deltamtime(channel_file))
    return channel_file
