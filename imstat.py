#!/usr/bin/env python
"""
Calculate statistical results for FITS images
"""

import re
import argparse
import logging
import os.path
from astropy.io import fits
from astropy import stats
import numpy as np

def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        description="Calculate statistical quantities for image")
    parser.add_argument("fitsfile", nargs='+',
                        metavar="file", help="input fits file(s)")
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
    parser.add_argument("--tearing", action='store_true',
                        help="add tearing metric to quicklook output")
    parser.add_argument("--dipoles", action='store_true',
                        help="add dipole metric to quicklook output")
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
        elif not optlist.noheadings: #- print filename
            print "#"
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
        print "{:>7.2g}".format(np.std(buf)),
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
        #- get header details to extract signal, bias regions
        pix = hdulist[hduid].data
        logging.debug("shape(pix)={}".format(np.shape(pix)))
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
        else:
            emsg = "DATASEC:{} parsing failed".format(dstr)
            logging.error(emsg)
            exit(1)
        naxis1 = int(hdr['NAXIS1'])
        naxis2 = int(hdr['NAXIS2'])
        #- define region to measure signal level
        x1 = int(datasec[0]) - 1
        x2 = int(datasec[1])
        y1 = int(datasec[2]) - 1
        y2 = int(datasec[3])
        sig_buf = pix[y1:y2, x1:x2]
        debugmsg = "sig_buf = pix[{}:{}, {}:{}]".format(y1, y2, x1, x2)
        logging.debug(debugmsg)
        debugmsg =  "shape(sig_buf)={}".format(np.shape(sig_buf))
        logging.debug(debugmsg)

        #Serial bias_region = "[y1:y2,x1:x2]"
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

        #parallel bias_region = "[y1:y2,x1:x2]"
        x1 = int(datasec[0]) - 1
        x2 = int(datasec[1])
        y1 = int(datasec[3])
        y2 = naxis2
        if y1 > y2 or x1 > x2:
            emsg = "No bias region available for datasec={}"\
                " with naxis1={}, naxis2={}".\
                format(datasec, naxis1, naxis2)
            logging.error(emsg)
            exit(1)
        p_bias_buf = pix[y1:y2, x1:x2]
        debugmsg = "p_bias_buf = pix[{}:{}, {}:{}]".format(y1, y2, x1, x2)
        logging.debug(debugmsg)
        debugmsg = "shape(p_bias_buf)={}".format(np.shape(p_bias_buf))
        logging.debug(debugmsg)
        quicklook_print(optlist, hduid, name, sig_buf,
                        bias_buf, p_bias_buf, expt)

def quicklook_print(optlist, sid, name, sig_buf,
                    bias_buf, p_bias_buf, expt):
    """perform and print the given statistics quantities
       fields are: mean, bias, signal, noise, adu/s
    """
    #quick_fields = ["mean", "bias", "signal", "noise", "adu/sec"]
    quick_fields = ["mean", "bias", "signal",
                    "noise", "adu/sec", "eper:s-cte", "eper:p-cte"]
    if optlist.tearing:
        quick_fields.append("tearing")
    if optlist.dipoles:
        quick_fields.append("dipoles")
    if not optlist.noheadings and ncalls.counter == 0:
        print "#{:>3s} {:>9s}".format("id", "HDUname"),
        if "mean" in  quick_fields:
            print " {:>6s}".format("median"),
        if "bias" in  quick_fields:
            print " {:>5s}".format("bias"),
        if "signal" in  quick_fields:
            print " {:>6s}".format("signal"),
        if "noise" in  quick_fields:
            print " {:>7s}".format("noise"),
        if "adu/sec" in  quick_fields and expt > 0:
            print "{:>8s}".format("adu/sec"),
        if "eper:s-cte" in  quick_fields:
            print "{:>9s}".format("s-cte"),
        if "eper:p-cte" in  quick_fields:
            print "{:>9s}".format("p-cte"),
        if "tearing" in  quick_fields:
            print "{:>15s}".format("tearing: L  R"),
        if "dipoles" in  quick_fields:
            print "{:>8s}".format("%dipoles"),
        print #-- newline

    if not optlist.noheadings:
        print " {:3d} {:>9s}".format(sid, name),

    if "mean" in  quick_fields:
        sig_mean = np.median(sig_buf)
        print " {:>6g}".format(sig_mean),
    if "bias" in  quick_fields:
        bias_mean = np.median(bias_buf)
        print " {:>5g}".format(bias_mean),
    if "signal" in  quick_fields:
        signal = sig_mean - bias_mean
        print " {:>6g}".format(signal),
    if "noise" in  quick_fields:
        (nrows, ncols) = np.shape(bias_buf)
        print " {:>7.4g}".format(np.std(
            bias_buf[int(nrows/2-nrows/20):int(nrows/2+nrows/20), 3:ncols-2])),
    if "adu/sec" in  quick_fields and expt > 0:
        print "{:>8.2f}".format(float(signal)/expt),
    if "eper:s-cte" in  quick_fields:
        debugmsg = "s-cte------------------"
        logging.debug(debugmsg)
        (nrows, ncols) = np.shape(sig_buf)
        nsig_cols = ncols
        #- define region to measure signal used in cte calc
        y1 = int(nrows/2-nrows/10)
        y2 = int(nrows/2+nrows/10)
        x0 = ncols-ncols/20
        x1 = ncols-1
        debugmsg = "s_n=sig_buf[{}:{},{}:{}]".format(y1,y2,x0,x1)
        logging.debug(debugmsg)
        s_n = sig_buf[y1:y2, x0:x1]
        (nrows, ncols) = np.shape(bias_buf)
        l_ncols = 3
        bias_mean = np.mean(bias_buf[y1:y2,l_ncols:ncols])
        debugmsg = "using bias_buf[{}:{},{}:{}]".format(y1,y2,l_ncols,ncols)
        logging.debug(debugmsg)
        l_n = np.mean(s_n) - bias_mean
        debugmsg = "l_n={:>10.6g}".format(l_n)
        logging.debug(debugmsg)
        y1 = int(nrows/2-nrows/10)
        y2 = int(nrows/2+nrows/10)
        x0 = 0
        x1 = l_ncols
        debugmsg = "b_n=bias_buf[{}:{},{}:{}]".format(y1,y2,x0,x1)
        logging.debug(debugmsg)
        b_n = bias_buf[y1:y2, x0:x1]
        l_nn = np.mean((b_n - bias_mean).sum(axis=1))
        debugmsg = "l_nn={:>10.6g}".format(l_nn)
        logging.debug(debugmsg)
        if l_n > 0.0:
            eper = 1 - (l_nn / (nsig_cols * l_n))
            print " {:>8.6g}".format(eper),
    if "eper:p-cte" in  quick_fields:
        debugmsg = "p-cte------------------"
        logging.debug(debugmsg)
        (nrows, ncols) = np.shape(p_bias_buf)
        l_nrows = 3 # number of overscan rows used to determing cte
        #- define region to measure bias used in cte calc
        x1 = int(ncols/2-ncols/10)
        x2 = int(ncols/2+ncols/10)
        y0 = l_nrows
        y1 = nrows
        p_bias_mean = np.mean(p_bias_buf[y0:y1,x1:x2])
        debugmsg = "p_bias_mean=mean(p_bias_buf[{}:{},{}:{}])".format(y0,y1,x1,x2)
        logging.debug(debugmsg)
        debugmsg = "p_bias_mean={:>10.6g}".format(p_bias_mean)
        logging.debug(debugmsg)
        #- define region to measure signal used in cte calc
        (nrows, ncols) = np.shape(sig_buf)
        nsig_rows = nrows
        x1 = int(ncols/2-ncols/10)
        x2 = int(ncols/2+ncols/10)
        y0 = nrows-100
        y1 = nrows-1
        debugmsg = "s_n=sig_buf[{}:{},{}:{}]".format(y0,y1,x1,x2)
        logging.debug(debugmsg)
        s_n = sig_buf[y0:y1, x1:x2]
        l_n = np.mean(s_n) - p_bias_mean
        debugmsg = "l_n={:>10.6g}".format(l_n)
        logging.debug(debugmsg)
        x1 = int(ncols/2-ncols/10)
        x2 = int(ncols/2+ncols/10)
        y0 = 0
        y1 = l_nrows
        debugmsg = "shape(p_bias_buf)={}".format(np.shape(p_bias_buf))
        logging.debug(debugmsg)
        debugmsg = "b_n=p_bias_buf[{}:{},{}:{}]".format(x1,x2,y0,y1)
        logging.debug(debugmsg)
        b_n = p_bias_buf[y0:y1, x1:x2]
        debugmsg = "shape(b_n)={}".format(np.shape(b_n))
        logging.debug(debugmsg)
        #l_nn = np.median(b_n) - p_bias_mean
        l_nn = np.mean((b_n - bias_mean).sum(axis=0))
        debugmsg = "l_nn={:>10.6g}".format(l_nn)
        logging.debug(debugmsg)
        (nrows, ncols) = np.shape(sig_buf)
        if l_n > 0.0:
            eper = 1 - (l_nn / (nrows * l_n))
            print " {:>8.6g}".format(eper),
    #---------
    if "tearing" in  quick_fields:
        debugmsg = "tearing check----------"
        logging.debug(debugmsg)
    #- column-1 into an array arr1
    #- column-2..N into 2nd array with dimension N-1 x Ncols arr2
    #- take median of 2nd array to make 1-D: arr3
    #- find stddev of arr3
    #- form (arr3 - arr1)/stddev as arr4
    #- find the first value of index "j" in sorted(arr4) less than 1.0
    #- report out (len(arr4)-j)/len(arr4) to 1 digit as tearing where
        #- this represents the fraction of the array less than 1.0
        #- left side
        arr3 = np.median(sig_buf[:,2:40], axis=1)
        arr4 = (arr3 - sig_buf[:,0])/np.std(arr3)
        tm = (1.0*np.size(arr4) - np.searchsorted(arr4, 1.0))/np.size(arr4)
        print "{:>4.1f}".format(tm),
        #- right side
        arr3 = np.median(sig_buf[:,-40:-2], axis=1)
        arr4 = (arr3 - sig_buf[:,-0])/np.std(arr3)
        tm = (1.0*np.size(arr4) - np.searchsorted(arr4, 1.0))/np.size(arr4)
        print "{:>4.1f}".format(tm),
    #---------
    if "dipoles" in  quick_fields:
        debugmsg = "dipoles check----------"
        logging.debug(debugmsg)
    #- region to work on is sig_buf, say 200 rows near top
    #- transpose to column order
    #- find sigma-clipped mean, median and stdev
    #- subtract mean from array
    #- divide the array by sigma
    #- go through array finding pixel pairs of differing sign
    #- and where |A(n)-A(n+1)| > 6
    #- add one to counter each time such a pair is found
    #- print out the % of pixels occupied by dipoles
        (nrows, ncols) = np.shape(sig_buf)
        arr1 = sig_buf[-nrows/10:-1,:] #- use top 10% of array
        debugmsg = "using subarray [{}:{},:]".format(-nrows/10,-1)
        logging.debug(debugmsg)
        arr2 = arr1.flatten('F') #- flatten to 1d in column order
        avg2, med2, std2 = stats.sigma_clipped_stats(arr2)
        debugmsg = "clipped stats: avg:{:>.3g} med:{} stdev:{:>.3g}".format(avg2, med2, std2)
        logging.debug(debugmsg)
        arr3 = (arr2 - avg2)/std2
        ndipole = 0
        for i in range(0, np.size(arr3) - 1):
            if (np.sign(arr3[i+1] * arr3[i]) == -1) and abs(arr3[i+1] - arr3[i]) > 5:
                ndipole += 1
        debugmsg = "dipole count = {}".format(ndipole)
        logging.debug(debugmsg)
        print "{:>8.2f}".format(100.0*float(2*ndipole)/(np.size(arr1)))
    print #-- newline
    ncalls() #-- track call count, acts like static variable

def ncalls():
    """maintain a counter
    """
    ncalls.counter += 1

if __name__ == '__main__':
    main()
