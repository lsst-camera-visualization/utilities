#!/usr/bin/env python
"""Print fits file headers with a variety of options
"""

import argparse
import logging
from astropy.io import fits
import imutils as iu


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="Print headers (Primary+Image)")
    parser.add_argument("fitsfile", nargs='+',
                        metavar="file", help="input fits file(s)")
    parser.add_argument("--info", action='store_true',
                        help="print the info() table summarizing file")
    parser.add_argument("--hduname", nargs='+',
                        metavar='idn', help="process HDU list by names")
    parser.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="process HDU list by ids")
    parser.add_argument("--all",
                        action='store_true', help="print all extensions")
    return parser.parse_args()
    # return opts


def main():
    """print out headers"""
    optlist = parse_args()
    iu.init_logging
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
    hduids = []
    if opts.hduname:
        for hduname in opts.hduname:
            try:
                hduid = hdulist.index_of(hduname)
                if hduid not in hduids:
                    hduids.append(hduid)
            except KeyError as ke:
                logging.error('KeyError: %s', ke)
                logging.error('HDU[%s] not found, skipping', hduname)

    if opts.hduindex:
        for hduid in opts.hduindex:
            try:
                hdu = hdulist[hduid]
                if hduid not in hduids:
                    hduids.append(hduid)
            except IndexError:
                logging.error('HDU[%d] not found, skipping', hduid)

    if opts.hduname or opts.hduindex:
        # just dump the ones asked for
        for hduid in sorted(hduids):
            hdu = hdulist[hduid]
            hdr = hdu.header
            print_single_hdr(text, hduid, hdu.name, hdr)
    else:
        for hdu in hdulist:
            hdr = hdu.header
            hduid = hdulist.index(hdu)
            if isinstance(hdu, (fits.PrimaryHDU, fits.ImageHDU,
                                fits.CompImageHDU)) or opts.all:
                print_single_hdr(text, hduid, hdu.name, hdr)

    hdulist.close()
    try:
        print('\n'.join(text))
    except IOError:
        # A 'Broken pipe' IOError may occur when stdout is closed prematurely
       pass


def print_single_hdr(text, hduid, name, hdr):
    """
    """
    text.append("#--------{}:{}--------".format(hduid, name))
    text.append(
        hdr.tostring(sep='\n', endcard=False, padding=False))


if __name__ == '__main__':
    main()
