#!/usr/bin/env python
# code to illustrate parsing of XML files
# importing the required modules
import csv
import os
import argparse
import logging
import collections
import requests
import re
from lxml import etree

def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="Parse and display structure of an xml file")
    parser.add_argument("target", metavar="URL|file",
                        help="filepath or URL")
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

    xmlstr = geturi(optlist.target)
    if xmlstr is None:
        logging.error("Failed to open uri: {}".format(uri))
        exit(1)
    printStructure(xmlstr)

def printStructure(xmlcontent):
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

def geturi(uri):
    """return the file either from internet or locally
    """
    res = re.match(r"http.*", uri)
    if res:
        #- we have a url, use request to return it
        try:
            resp = requests.head(uri)
        except requests.ConnectionError as e:
            logging.error("ConnectionError: {}".format(e))
            return
        if resp.status_code != 200:
            logging.error("error: invalid response {} from Server".format(resp.status_code))
            return
        resp = requests.get(uri)
        return resp.content
    else:
        #- it may be a file, try to open as a path
        if not os.path.exists(uri):
            logging.error("{} is not a path, exiting...".format(uri))
            return
        statinfo = os.stat(uri)
        mode = statinfo.st_mode
        try:
            xmlfile = open(uri, mode='r')
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            return
        xmlstring = xmlfile.read()
        xmlfile.close()
        return xmlstring

if __name__ ==  '__main__':
    main()


