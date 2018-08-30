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

def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="list Trending channels from CCS localdb")
    parser.add_argument("ssname", nargs='+',
                        help="CCS subsystem name(s)")
    parser.add_argument("--debug", action='store_true',
                        help="print additional debugging info")
    parser.add_argument("--structure", action='store_true',
                        help="print structure of xml")
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

    channelFile = updateTrendingChannelsXML()
    chfile = open(channelFile, mode='rb')
    channels_xml = chfile.read()
    chfile.close()
    if optlist.structure:
        printChannelStructure(channels_xml)
        return
    printChannelContent(channels_xml, optlist.ssname)

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
        #if delta < 86400:
        if delta < 3600:
            return channelFile
        logging.info("updating cached channelFile...")
        logging.debug("initial file: age={} (s)".format(deltamtime(channelFile)))
    xmlstring = getAllChannels()
    if xmlstring:
        chfile = open(channelFile, mode='wb')
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


