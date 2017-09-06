#!/usr/bin/env python
"""
Calculate statistical results for FITS images
"""

import os.path
import re
import argparse
from astropy.io import fits
import numpy as np

def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="Calculate statistical quantities for image")
    parser.add_argument("fitsfile", nargs='+',
                        metavar="file", help="input fits file")
    parser.add_argument("--info", action='store_true',
                        help="print the info() table summarizing file")
    parser.add_argument("--noheadings", action='store_true', default=False,
                        help="Don't print column heads for stats")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--hduname", nargs='+',
                       metavar='idn', help="process HDU list by names")
    group.add_argument("--hduindex", nargs='+', type=int,
                       metavar='idx', help="process HDU list by ids")
    parser.add_argument("--region", nargs='+', metavar='reg',
                        help="region fmt: \"x1:x2,y1:y2\"")
    parser.add_argument("--stats", nargs='+', metavar='stat',
                        help="select: mean median stddev min max")
    parser.add_argument("--quicklook", action='store_true',
                        help="estimate signal, noise, counts/sec in adus")
    return parser.parse_args()

def stat_print(hdulist):
    """print statistics for region according to options
    """
    hduids = []  #-- make a list of HDU ids to work on
    if optlist.hduname:
        for hduname in optlist.hduname:
            hduids.append(hdulist.index_of(hduname))
    elif optlist.hduindex:
        hduids = optlist.hduindex
    else: #- all segments
        for hdu in hdulist:
            if isinstance(hdu, fits.ImageHDU):
                hduids.append(hdulist.index_of(hdu.name))

    for hduid in hduids: # process each with optional region
        pix = hdulist[hduid].data
        name = hdulist[hduid].name
        reg = ""
        if not optlist.region:
            do_stats(hduid, name, pix, reg)
        else:
            for  reg in optlist.region:
                res = re.match(r"\[*([0-9]*):([0-9]+),([0-9]+):([0-9]+)\]*", reg)
                if res:
                    (x1, x2, y1, y2) = res.groups()
                    do_stats(hduid, name, pix[int(y1):int(y2),
                                            int(x1):int(x2)], reg)
                else: #- region = [*,y1:y2]
                    res = re.match(r"\[*(\*),([0-9]+):([0-9]+)\]*", reg)
                    if res:
                        (x, y1, y2) = res.groups()
                        do_stats(hduid, name, pix[int(y1):int(y2),:], reg)
                    else: # reg = [x1:x2,*]
                        res = re.match(r"\[*([0-9]+):([0-9]+),(\*)\]*", reg)
                        if res:
                            (x1, x2, y) = res.groups()
                            do_stats(hduid, name, pix[:,int(x1):int(x2)], reg)
                        else:
                            print "bad region spec {}, no match produced".\
                                format(res.string())
                            exit(1)

def do_stats(sid, name, buf, reg):
    """perform and print the given statistics quantities
    """
    if not optlist.stats:
        optlist.stats = ["mean", "median", "stddev", "min", "max"]
    if not optlist.noheadings and ncalls.counter == 0:
        print "#{:>3s} {:>9s}".format("id", "HDUname"),
        if "mean" in  optlist.stats:
            print "{:>9s}".format("mean"),
        if "median" in  optlist.stats:
            print "{:>9s}".format("median"),
        if "stddev" in  optlist.stats:
            print "{:>7s}".format("stddev"),
        if "min" in  optlist.stats:
            print "{:>8s}".format("min"),
        if "max" in  optlist.stats:
            print "{:>8s}".format("max"),
        if len(reg):
            print " {:21s}".format("region"),
        print #-- newline

    if not optlist.noheadings:
        print " {:3d} {:>9s}".format(sid, name),

    if "mean" in  optlist.stats:
        print "{:>9g}".format(np.mean(buf)),
    if "median" in  optlist.stats:
        print "{:>9g}".format(np.median(buf)),
    if "stddev" in  optlist.stats:
        print "{:>7.1f}".format(np.std(buf)),
    if "min" in  optlist.stats:
        print "{:>8g}".format(np.min(buf)),
    if "max" in  optlist.stats:
        print "{:>8g}".format(np.max(buf)),
    if len(reg):
        print " {:21s}".format(reg),
    print #-- newline
    ncalls() #-- track call count, acts like static variable

def quicklook_print(hdulist):
    """print quicklook for hdu according to options
    """
    hdr = hdulist[0].header
    try:
        expt = float(hdr['EXPTIME'])
    except KeyError as ke:
        print "Key error: {}".format(ke)
        print "adu/sec won't be available"
        expt = 0.0
    hduids = []  #-- ids of hdus to operate on
    if optlist.hduname:
        for hduname in optlist.hduname:
            hduids.append(hdulist.index_of(hduname))
    elif optlist.hduindex:
        hduids = optlist.hduindex
    else: #- all segments
        for hdu in hdulist:
            if isinstance(hdu, fits.ImageHDU):
                hduids.append(hdulist.index_of(hdu.name))

    for hduid in hduids:
        pix = hdulist[hduid].data
        name = hdulist[hduid].name
        hdr = hdulist[hduid].header
        try:
            dstr = hdr['DATASEC']
        except KeyError as ke:
            print "Key error: {}".format(ke)
            print "DATASEC required for quicklook, exiting..."
            exit(1)
                #print "datasec={}".format(dstr)
        res = re.match(r"\[*([0-9]*):([0-9]+),([0-9]+):([0-9]+)\]*", dstr)
        if res:
            datasec = res.groups()
        naxis1 = int(hdr['NAXIS1'])
        naxis2 = int(hdr['NAXIS2'])
        #sig_region = "[y1:y2,x1:x2]"
        y1 = max(int(int(datasec[3])/2) - 150, int(datasec[2]))
        y2 = min(y1+300, int(datasec[3]))
        #x1 = max(int(int(datasec[0])/2) - 150, int(datasec[0]))
        x1 = max(int(int(datasec[1]) - int(datasec[0]))/2 - 150,
                 int(datasec[0]))
        x2 = min(x1+300, int(datasec[1]))
        #print "sig_buf: x1={}, x2={}, y1={}, y2={}".format(x1,x2,y1,y2)
        sig_buf = pix[y1:y2,x1:x2]

        #bias_region = "[y1:y2,x1:x2]"
        y1 = max(int(int(datasec[3])/2) - 150, int(datasec[2]))
        y2 = min(y1+300, int(datasec[3]))
        x1 = int(datasec[1]) + 2
        x2 = naxis1 - 2
        #print "bias_buf x1={}, x2={}, y1={}, y2={}".format(x1,x2,y1,y2)
        if y1 > y2 or x1 > x2:
            print "No bias region available for datasec={} with naxis1={}, naxis2={}".\
                format(datasec,naxis1,naxis2)
        bias_buf = pix[y1:y2,x1:x2]
        do_quicklook(hduid, name, sig_buf, bias_buf, expt)

def do_quicklook(sid, name, sig_buf, bias_buf, expt):
    """perform and print the given statistics quantities
       fields are: mean, bias, signal, noise, adu/s
    """
    quick_fields = ["mean", "bias", "signal", "noise", "adu/sec"]
    if not optlist.noheadings and ncalls.counter == 0:
        print "#{:>3s} {:>9s}".format("id", "HDUname"),
        if "mean" in  quick_fields:
            print "{:>10s}".format("[adu] mean"),
        if "bias" in  quick_fields:
            print "{:>9s}".format("bias"),
        if "signal" in  quick_fields:
            print "{:>9s}".format("signal"),
        if "noise" in  quick_fields:
            print "{:>8s}".format("noise"),
        if "adu/sec" in  quick_fields and expt > 0:
            print "{:>8s}".format("adu/sec"),
        print #-- newline

    if not optlist.noheadings:
        print " {:3d} {:>9s}".format(sid, name),

    if "mean" in  quick_fields:
        sig_mean = np.median(sig_buf)
        print "{:>10g}".format(sig_mean),
    if "bias" in  quick_fields:
        bias_mean = np.median(bias_buf)
        print "{:>9g}".format(bias_mean),
    if "signal" in  quick_fields:
        signal = sig_mean - bias_mean
        print "{:>9g}".format(signal),
    if "noise" in  quick_fields:
        print "{:>8.2f}".format(np.std(bias_buf)),
    if "adu/sec" in  quick_fields and expt > 0:
        print "{:>8.2f}".format(signal/expt),
    print #-- newline
    ncalls() #-- track call count, acts like static variable

def ncalls():
    """maintain a counter
    """
    ncalls.counter += 1

if __name__ == '__main__':
    optlist = parse_args()
    ncalls.counter = 0
    # begin processing
    for ffile in optlist.fitsfile:
        try:
            hdulist = fits.open(ffile)
        except IOError as ioerr:
            print "I/O error: {}".format(ioerr.strerror)
            exit(1)
        if optlist.info: # just print the image info and exit
            hdulist.info()
            exit(0)
        if not optlist.noheadings:
            print "#"
            print "# {}".format(hdulist[0].header['filename'])
        if optlist.quicklook:
            quicklook_print(hdulist)
        else:
            stat_print(hdulist)
        ncalls.counter = 0




