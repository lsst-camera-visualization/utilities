#!/usr/bin/env python
"""Print fits file headers with a variety of options
"""

import argparse
from astropy.io import fits


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="Print headers (Primary+Image)")
    parser.add_argument("fitsfile", help="input fits file")
    parser.add_argument("--info", action='store_true',
                        help="print the info() table summarizing file")
    parser.add_argument("--hduname",
                        help="dump only the named extension header")
    parser.add_argument("--hduindex", type=int,
                        help="dump only id'd extension header")
    parser.add_argument("--unsorted", action='store_true',
                        help="print headers in index order")
    parser.add_argument("--all",
                        action='store_true', help="print all extensions")
    return parser.parse_args()
    # return opts


def header_print(opts):
    """print the headers according to the options
    """
    seglist = {}      # dict to hold segment name, index as k,v
    otherlist = {}    # dict to hold non-image HDU's (except for primary)
    # opts = parse_args()
    hdulist = fits.open(opts.fitsfile)
    # just print the image info with indexes, names, sizes
    if opts.info:
        hdulist.info()
    # print single hdu by name
    elif opts.hduname:
        index = hdulist.index_of(opts.hduname)
        hdu = hdulist[index]
        print "#--------{}---------".format(opts.hduname)
        print hdu.header.tostring(sep='\n', endcard=False, padding=False)
    # print single hdu by index
    elif opts.hduindex:
        hdu = hdulist[opts.hduindex]
        print "#--------extension {}---------".format(opts.hduindex)
        print hdu.header.tostring(sep='\n', endcard=False, padding=False)
    # print primary and Image headers, others optionally
    else:
        # build dicts etc. to facilitate processing
        for hdu in hdulist:
            index = hdulist.index_of(hdu.name)
            if isinstance(hdu, fits.ImageHDU):
                seglist[hdu.name] = index
            elif isinstance(hdu, fits.PrimaryHDU):
                pindex = index
            else:
                otherlist[hdu.name] = index
        # print the primary
        hdu = hdulist[pindex]
        print "#--------{}---------".format(hdu.name)
        print hdu.header.tostring(sep='\n', endcard=False, padding=False)
        # print the Image headers
        if opts.unsorted:
            # in original index order
            # for hdu in hdulist:
            ids = list(seglist.values())
            for index in sorted(ids):
                hdu = hdulist[index]
                print "#--------{}---------".format(hdu.name)
                print hdu.header.tostring(
                    sep='\n', endcard=False, padding=False)
        else:
            # sorted by name (default)
            for name, index in sorted(seglist.iteritems()):
                hdu = hdulist[index]
                print "#--------{}---------".format(name)
                print hdu.header.tostring(
                    sep='\n', endcard=False, padding=False)
        # print the other headers
        if opts.all:
            ids = list(otherlist.values())
            for index in sorted(ids):
                hdu = hdulist[index]
                print "#--------{}---------".format(hdu.name)
                print hdu.header.tostring(sep='\n',
                                          endcard=False, padding=False)
    hdulist.close()

if __name__ == '__main__':
    optlist = parse_args()
    header_print(optlist)
