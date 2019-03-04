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
    """get list of channels from the server
       maxidle recovers channels active within maxidle seconds
    """
    trending_server = get_trending_server()
    if not trending_server:
        return None
    listpathroot = "8080/rest/data/dataserver/listchannels?maxIdleSeconds="
    if maxidle:
        listpath = "{}{:d}".format(listpathroot, int(maxidle))
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
        logging.debug('parse_date_str(datestr=%s)', datestr)
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
    else:
        duration = convert_to_seconds(duration)

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


def update_trending_channels_xml(tstart=None, tstop=None):
    """maintain local cache of trending channels in xml file
       arg: tstart -- channels active since tstart (seconds since the epoch)
    """
    logging.debug('update_trending_channels_xml(%s, %s)', tstart, tstop)
    cachedir = "{}/.trender".format(os.environ.get('HOME'))
    channel_file = "{}/.trender/listchannels.xml".format(
        os.environ.get('HOME'))
    update = True
    # check channel_file exists, get mtime, update if need be
    if not os.path.exists(cachedir):   # make cachdir if not exist
        os.mkdir(cachedir)
    if not os.path.isdir(cachedir):    # is not a directory
        logging.error('%s is not a directory, exiting...', cachedir)
        exit(1)
    #
    # Trigger an update based on whether the interval (tstart,tstop)
    # is covered in the existing cached file based on time of last update
    # and maxIdleSeconds given when the file was fetched.  The attribute
    # 'activeSinceDate' is then maxIdleSeconds prior to last update time
    #
    if os.path.exists(channel_file):  # trigger update based on mtime
        statinfo = os.stat(channel_file)
        mode = statinfo.st_mode
        if not stat.S_IWUSR & mode:  # not writeable
            os.chmod(channel_file, mode | stat.S_IWUSR)
        delta = int(time.time() - statinfo.st_mtime)
        logging.debug('existing cache age: %d (s)', delta)

        chfile = open(channel_file, mode='rb')
        xmlstring = chfile.read()
        chfile.close()
        root = etree.fromstring(xmlstring)
        active_since = root.attrib.get('activeSinceDate')
        if active_since:  # parse, convert and compare to tstart
            xml_start_epoch = parse_datestr(active_since)
            logging.debug('%s channels active_since: %s', channel_file,
                          datetime.fromtimestamp(xml_start_epoch,
                                                 gettz(tz_trending)))
        # If tstart is inside the interval: [xml_start, last_update]
        # then can guarantee desired channels were being published
        # and hence are already in the cached file.
        if tstart and xml_start_epoch < tstart < statinfo.st_mtime:
            if not tstop or tstop < (statinfo.st_mtime + 86400):
                update = False

    if update:
        logging.info('updating cached channel_file...')
        if tstart:
            xstart = tstart - 86400  # adjust to 24h earlier
            maxidle = int(time.time() - xstart)
        else:
            maxidle = 3600 * 24 * 7  # give it a week
            xstart = int(time.time() - maxidle)
        xmlstring = get_all_channels(maxidle)
        root = etree.fromstring(xmlstring)
        active_since = root.attrib.get('activeSinceDate')
        if not active_since:  # needed until service adds this attrib
            active_since_str = datetime.fromtimestamp(xstart).isoformat(
                timespec='seconds')
            root.set("activeSinceDate", active_since_str)
        if xmlstring:
            tree = etree.ElementTree(root)
            tree.write(channel_file, xml_declaration=True,
                       encoding='UTF-8', pretty_print=False,
                       standalone="yes")
    else:
        logging.debug('returning existing channel_file=%s', channel_file)
    return channel_file


def convert_to_seconds(duration_str):
    """
    return duration in seconds
    """

    seconds = 0
    if re.match(r"[0-9]+$", duration_str):
        seconds = int(duration_str[:-1])
    elif re.match(r"[0-9]+s$", duration_str):
        seconds = int(duration_str[:-1])
    elif re.match(r"[0-9]+m$", duration_str):
        seconds = 60 * int(duration_str[:-1])
    elif re.match(r"[0-9]+h$", duration_str):
        seconds = 3600 * int(duration_str[:-1])
    elif re.match(r"[0-9]+d$", duration_str):
        seconds = 84600 * int(duration_str[:-1])

    return seconds
