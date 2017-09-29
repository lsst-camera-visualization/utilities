#!/usr/bin/env python
"""
Calculate statistical results for FITS images
"""

import re
import argparse
import logging
import os.path
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
    parser.add_argument("--debug", action='store_true',
                        help="print additional debugging messages")
    parser.add_argument("--noheadings", action='store_true',
                        default=False, help="Don't print column heads for stats")
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument("--hduname", nargs='+',
                        metavar='idn', help="process HDU list by names")
    hgroup.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="process HDU list by ids")
    sgroup = parser.add_argument_group("stats", "select statistics and regions"
                                       " (exclusive of quicklook)")
    sgroup.add_argument("--region", nargs='+', metavar='reg',
                        help="region fmt: \"x1:x2,y1:y2\",")
    sgroup.add_argument("--stats", nargs='+', metavar='stat',
                        help="select: mean median stddev min max")
    parser.add_argument("--quicklook", action='store_true',
                        help="estimate signal, noise, counts/sec in adus")
    return parser.parse_args()

def main():
    """main logic:"""
    optlist = parse_args()
    if optlist.debug:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.INFO)
    ncalls.counter = 0
    # begin processing -- loop over files
    for ffile in optlist.fitsfile:
        try:
            hdulist = fits.open(ffile)
        except IOError as ioerr:
            emsg = "IOError: {}".format(ioerr)
            logging.error(emsg)
            exit(1)
        if optlist.info: # just print the image info and exit
            hdulist.info()
            exit(0)
        if not optlist.noheadings: #- print filename
            print "#"
            #- save this for later review
            #print "# {}".format(hdulist[0].header['filename'])
            #print "fitsfile={}".format(ffile)
            #print os.path.basename(ffile)
            print "# {}".format(os.path.basename(ffile))
        if optlist.quicklook:
            quicklook(optlist, hdulist)
        else:
            stats_proc(optlist, hdulist)
        ncalls.counter = 0


def stats_proc(optlist, hdulist):
    """print statistics for region according to options
    """
    hduids = []  #-- make a list of HDU ids to work on
    if optlist.hduname:
        for hduname in optlist.hduname:
            try:
                hduids.append(hdulist.index_of(hduname))
            except KeyError as ke:
                emsg = "KeyError: {}".format(ke)
                logging.error(emsg)
                exit(1)

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
            stats_print(optlist, hduid, name, pix, reg)
        else:
            for reg in optlist.region:
                res = re.match(
                    r"\[*([0-9]*):([0-9]+),([0-9]+):([0-9]+)\]*",
                    reg)
                if res:
                    (x1, x2, y1, y2) = res.groups()
                    stats_print(optlist, hduid, name,
                                pix[int(y1)-1:int(y2),
                                    int(x1)-1:int(x2)], reg)
                    continue
                #- reg = [x0,y1:y2] -- single column
                res = re.match(r"\[*([0-9]+),([0-9]+):([0-9]+)\]*",
                               reg)
                if res:
                    (x0, y1, y2) = res.groups()
                    stats_print(optlist, hduid, name,
                                pix[int(y1)-1:int(y2), int(x0)-1], reg)
                    continue
                #- reg = [*,y1:y2]
                res = re.match(r"\[*(\*),([0-9]+):([0-9]+)\]*", reg)
                if res:
                    (x, y1, y2) = res.groups()
                    stats_print(optlist, hduid, name,
                                pix[int(y1)-1:int(y2), :], reg)
                    continue
                # reg = [x0,y0] -- single pixel
                res = re.match(r"\[*([0-9]+),([0-9]+)\]*", reg)
                if res:
                    (x0, y0) = res.groups()
                    stats_print(optlist, hduid, name,
                                pix[int(y0)-1, int(x0)-1], reg)
                    continue
                # reg = [x1:x2,y0] -- single row
                res = re.match(r"\[*([0-9]+):([0-9]+),([0-9]+)\]*", reg)
                if res:
                    (x1, x2, y0) = res.groups()
                    stats_print(optlist, hduid, name,
                                pix[int(y0)-1, int(x1)-1:int(x2)], reg)
                    myarr = pix[int(y0)-1, int(x1)-1:int(x2)]
                    print "npix={}".format(myarr.size)
                    print "pix={}".format(myarr)
                    continue
                # reg = [x1:x2,*]
                res = re.match(r"\[*([0-9]+):([0-9]+),(\*)\]*", reg)
                if res:
                    (x1, x2, y) = res.groups()
                    stats_print(optlist, hduid, name,
                                pix[:, int(x1)-1:int(x2)], reg)
                    continue
                # reg = [*,*] #- redundant, but for completeness
                res = re.match(r"\[*(\*),(\*)\]*", reg)
                if res:
                    (x, y) = res.groups()
                    stats_print(optlist, hduid, name, pix[:, :], reg)
                    continue
                else:
                    emsg = "bad region spec {}, no match produced".\
                        format(reg)
                    logging.error(emsg)
                    exit(1)

def stats_print(optlist, sid, name, buf, reg):
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
            print "{:>7s}".format("min"),
        if "max" in  optlist.stats:
            print "{:>7s}".format("max"),
        if reg:
            print " {:20s}".format("region"),
        print #-- newline

    if not optlist.noheadings:
        print " {:3d} {:>9s}".format(sid, name),

    if "mean" in  optlist.stats:
        print "{:>9g}".format(np.mean(buf)),
    if "median" in  optlist.stats:
        print "{:>9g}".format(np.median(buf)),
    if "stddev" in  optlist.stats:
        print "{:>7.1g}".format(np.std(buf)),
    if "min" in  optlist.stats:
        print "{:>7g}".format(np.min(buf)),
    if "max" in  optlist.stats:
        print "{:>7g}".format(np.max(buf)),
    if reg:
        print " {:20s}".format(reg),
    print #-- newline
    ncalls() #-- track call count, acts like static variable

def quicklook(optlist, hdulist):
    """print quicklook for hdu according to options
    """
    hdr = hdulist[0].header
    try:
        expt = float(hdr['EXPTIME'])
    except KeyError as ke:
        emsg = "KeyError: {}".format(ke)
        logging.warn(emsg)
        emsg = "adu/sec won't be available"
        logging.warn(emsg)
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
        #debugmsg = "shape(pixl)={}".format(np.shape(pix))
        #logging.debug(debugmsg)
        logging.debug("shape(pixl)={}".format(np.shape(pix)))
        name = hdulist[hduid].name
        hdr = hdulist[hduid].header
        try:
            dstr = hdr['DATASEC']
        except KeyError as ke:
            emsg = "KeyError: {}, required for quicklook mode".format(ke)
            logging.error(emsg)
            exit(1)
        debugmsg = "DATASEC={}".format(dstr)
        logging.debug(debugmsg)
        res = re.match(r"\[*([0-9]*):([0-9]+),([0-9]+):([0-9]+)\]*",
                       dstr)
        if res:
            datasec = res.groups()
        naxis1 = int(hdr['NAXIS1'])
        naxis2 = int(hdr['NAXIS2'])
        #sig_region = "[y1:y2,x1:x2]"
        x1 = int(datasec[0]) - 1
        x2 = int(datasec[1])
        y1 = int(datasec[2]) - 1
        y2 = int(datasec[3])
        sig_buf = pix[y1:y2, x1:x2]
        debugmsg = "sig_buf = pix[{}:{}, {}:{}]".format(y1, y2, x1, x2)
        logging.debug(debugmsg)
        debugmsg =  "shape(sig_buf)={}".format(np.shape(sig_buf))
        logging.debug(debugmsg)

        #bias_region = "[y1:y2,x1:x2]"
        x1 = int(datasec[1])
        x2 = naxis1
        y1 = int(datasec[2]) - 1
        y2 = int(datasec[3])
        if y1 > y2 or x1 > x2:
            emsg = "No bias region available for datasec={}"\
                " with naxis1={}, naxis2={}".\
                format(datasec, naxis1, naxis2)
            logging.error(emsg)
            exit(1)
        bias_buf = pix[y1:y2, x1:x2]
        debugmsg = "bias_buf = pix[{}:{}, {}:{}]".format(y1, y2, x1, x2)
        logging.debug(debugmsg)
        debugmsg = "shape(bias_buf)={}".format(np.shape(bias_buf))
        logging.debug(debugmsg)
        quicklook_print(optlist, hduid, name, sig_buf, bias_buf, expt)

def quicklook_print(optlist, sid, name, sig_buf, bias_buf, expt):
    """perform and print the given statistics quantities
       fields are: mean, bias, signal, noise, adu/s
    """
    #quick_fields = ["mean", "bias", "signal", "noise", "adu/sec"]
    quick_fields = ["mean", "bias", "signal", "noise", "adu/sec", "eper:cte"]
    if not optlist.noheadings and ncalls.counter == 0:
        print "#{:>3s} {:>9s}".format("id", "HDUname"),
        if "mean" in  quick_fields:
            print "{:>8s}".format("median"),
        if "bias" in  quick_fields:
            print "{:>8s}".format("bias"),
        if "signal" in  quick_fields:
            print "{:>9s}".format("signal"),
        if "noise" in  quick_fields:
            print "{:>8s}".format("noise"),
        if "adu/sec" in  quick_fields and expt > 0:
            print "{:>8s}".format("adu/sec"),
        if "eper:cte" in  quick_fields:
            print "{:>10s}".format("eper:cte"),
        print #-- newline

    if not optlist.noheadings:
        print " {:3d} {:>9s}".format(sid, name),

    if "mean" in  quick_fields:
        sig_mean = np.median(sig_buf)
        print "{:>8g}".format(sig_mean),
    if "bias" in  quick_fields:
        bias_mean = np.median(bias_buf)
        print "{:>8g}".format(bias_mean),
    if "signal" in  quick_fields:
        signal = sig_mean - bias_mean
        print "{:>9g}".format(signal),
    if "noise" in  quick_fields:
        (nrows, ncols) = np.shape(sig_buf)
        print "{:>8.2f}".format(np.std(
            bias_buf[int(nrows/2-100):int(nrows/2+100), 2:ncols-2])),
    if "adu/sec" in  quick_fields and expt > 0:
        print "{:>8.2f}".format(float(signal)/expt),
    if "eper:cte" in  quick_fields:
        (nrows, ncols) = np.shape(sig_buf)
        y1 = int(nrows/2-100)
        y2 = int(nrows/2+100)
        x0 = ncols-1
        debugmsg = "s_n=sig_buf[{}:{},{}]".format(y1,y2,x0)
        logging.debug(debugmsg)
        s_n = sig_buf[y1:y2, x0]
        l_n = np.mean(s_n) - bias_mean
        y1 = int(nrows/2-100)
        y2 = int(nrows/2+100)
        x0 = 0
        x1 = 2
        debugmsg = "b_n=bias_buf[{}:{},{}]".format(y1,y2,x0,x1)
        logging.debug(debugmsg)
        b_n = bias_buf[y1:y2, x0:x1]
        l_nn = np.mean(b_n) - bias_mean
        eper = 1 - (l_nn / (ncols * l_n))
        print "{:>10.6g}".format(eper),
    print #-- newline
    ncalls() #-- track call count, acts like static variable

def ncalls():
    """maintain a counter
    """
    ncalls.counter += 1

if __name__ == '__main__':
    main()
