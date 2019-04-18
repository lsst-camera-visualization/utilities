#!/usr/bin/env python
"""Print fits file headers with a variety of options
"""

import argparse
import logging
from astropy.io import fits


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="Print headers (Primary+Image)")
    parser.add_argument("fitsfile", nargs='+',
                        metavar="file", help="input fits file(s)")
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


def main():
    """print out headers"""
    optlist = parse_args()
    # processing -- loop over files
    for ffile in optlist.fitsfile:
        try:
            hdulist = fits.open(ffile)
        except IOError as ioerr:
            emsg = "IOError: {}".format(ioerr)
            logging.error(emsg)
            exit(1)
        if optlist.info:  # just print the image info and exit
            hdulist.info()
        else:
            header_print(optlist, hdulist)
        hdulist.close()
    exit(0)


def header_print(opts, hdulist):
    """print the headers according to the options
    """
    text = []         # list to hold all the text to print at end
    seglist = {}      # dict to hold segment name, index as k,v
    otherlist = {}    # dict to hold non-image HDU's (except for primary)
    if opts.hduname:
        index = hdulist.index_of(opts.hduname)
        hdr = hdulist[opts.hduindex].header
        text.append("#--------{}---------".format(opts.hduname))
        text.append(
            hdr.tostring(sep='\n', endcard=False, padding=False))
    # print single hdu by index
    elif opts.hduindex:
        hdr = hdulist[opts.hduindex].header
        text.append("#--------extension {}---------".format(opts.hduindex))
        text.append(
            hdr.tostring(sep='\n', endcard=False, padding=False))
    # print primary and Image headers, others optionally
    else:
        # build dicts etc. to facilitate processing
        for hdu in hdulist:
            index = hdulist.index_of(hdu.name)
            if isinstance(hdu, (fits.ImageHDU, fits.CompImageHDU)):
                seglist[hdu.name] = index
            elif isinstance(hdu, fits.PrimaryHDU):
                pindex = index
            else:
                otherlist[hdu.name] = index
        # print the primary
        hdr = hdulist[pindex].header
        text.append("#--------{}---------".format(hdu.name))
        text.append(
            hdr.tostring(sep='\n', endcard=False, padding=False))
        # print the Image headers
        if opts.unsorted:
            # in original index order
            # for hdu in hdulist:
            ids = list(seglist.values())
            for index in sorted(ids):
                hdr = hdulist[index].header
                text.append(hdr.tostring(
                    sep='\n', endcard=False, padding=False))
        else:
            # sorted by name (default)
            for name, index in sorted(seglist.items()):
                hdr = hdulist[index].header
                text.append("#--------{}---------".format(name))
                text.append(hdr.tostring(
                    sep='\n', endcard=False, padding=False))
        # print the other headers
        if opts.all:
            ids = list(otherlist.values())
            for index in sorted(ids):
                hdr = hdulist[index].header
                text.append("#--------{}---------".format(hdu.name))
                text.append(hdr.tostring(sep='\n',
                                         endcard=False, padding=False))
    hdulist.close()
    try:
        print('\n'.join(text))
    except IOError:
        # A 'Broken pipe' IOError may occur when stdout is closed prematurely
        pass


if __name__ == '__main__':
    main()
