import socket
import re
import logging
import collections
import os.path
import time
import stat
import requests
from lxml import etree
import dateutil.parser as dp

# constants? for lack of a better term
slacnetworkre = r"134\.79\.[0-9]*\.[0-9]*"
slac_trending_server = "lsst-mcm.slac.stanford.edu"


def deltamtime(pathname):
    """ time since last mod
    """
    return time.time() - os.stat(pathname)[stat.ST_MTIME]


def get_all_channels():
    """pulls all the channels
    """
    trending_server = get_trending_server()
    listpath = "8080/rest/data/dataserver/listchannels?maxIdleSeconds=-1"
    url = "http://{}:{}".format(trending_server, listpath)

    # creating HTTP response object from given url
    try:
        resp = requests.head(url)
    except requests.ConnectionError as e:
        logging.error('ConnectionError: %s', e)
        logging.error('check status of ssh tunnel to trending server')
        return None
    if resp.status_code != 200:
        logging.error('error: invalid response %s from Trending Server',
                      resp.status_code)
        return resp.status_code
    resp = requests.get(url)
    return resp.content


def parse_datestr(datestr):
    """parse date string to return time in seconds
       or now if the string is None
    """
    if datestr:  # parse the date string almost any format
        logging.debug('datestr=%s', datestr)
        pt = dp.parse(datestr)
        if not pt.tzinfo:               # no time zone in string
            ts = int(time.mktime(pt.timetuple()))
        elif pt.tzinfo and pt.tzname():  # tz name defined
            ts = int(time.mktime(pt.timetuple()))
        elif pt.tzinfo and not pt.tzname():  # tz name not defined
            ts = int(time.mktime(pt.timetuple())) - 3600*time.daylight
        else:
            logging.error('could not parse the datestring')
            ts = None
    else:
        ts = int(time.time())  # right now
    return ts


def get_time_interval(startstr, stopstr):
    """return the timeinterval boundaries (ms) from datestrings
    """
    t1 = parse_datestr(startstr)
    if t1:
        logging.debug('t1 as localtime: %s', time.strftime(
            "%a, %d %b %Y %H:%M:%S", time.localtime(t1)))
        t1 *= 1000

    t2 = parse_datestr(stopstr)
    if t2:
        logging.debug('t2 as localtime: %s', time.strftime(
            "%a, %d %b %Y %H:%M:%S", time.localtime(t2)))
        t2 *= 1000

    # cases for t1, t2
    if startstr and stopstr:
        if t1 > t2:
            logging.error('starttime must be earlier than stop')
            t1 = None
            t2 = None
    elif startstr and not stopstr:  # get 300s interval
        t2 = t1 + 300*1000
    elif not startstr and stopstr:  # get 300s interval
        t1 = t2 - 300*1000
    elif not startstr and not stopstr:
        t1 = t2 - 300*1000
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
    res = re.search(slacnetworkre, s.getsockname()[0])
    if res:
        trending_server = slac_trending_server
        logging.debug('on SLAC network: trending server is lsst-mcm')
    else:
        trending_server = "localhost"
        logging.debug('off SLAC network: trending server is localhost')
    # test connection to trending server, just get head to verify
    #
    try:
        requests.head("http://{}:8080".format(trending_server))
    except requests.ConnectionError as e:
        logging.error('ConnectionError: %s', e)
        logging.error('check status of ssh tunnel to trending server')
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
    else:
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



def update_trending_channels_xml():
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
    if os.path.exists(channel_file):  # file exists or not
        statinfo = os.stat(channel_file)
        mode = statinfo.st_mode
        if not stat.S_IWUSR & mode:  # not writeable
            os.chmod(channel_file, mode | stat.S_IWUSR)
        delta = time.time() - statinfo.st_mtime
        logging.debug('found existing cache age: %s (s)', delta)
        if delta < 86400:
            return channel_file
        logging.info('updating cached channel_file...')
        logging.debug('initial file: age=%s (s)', deltamtime(channel_file))
    xmlstring = get_all_channels()
    if xmlstring:
        chfile = open(channel_file, mode='wb')
        chfile.write(xmlstring)
        chfile.close()
    else:
        logging.info('failed to update cached channel_file, continuing')
    logging.debug('final file: age=%s (s)', deltamtime(channel_file))
    return channel_file
