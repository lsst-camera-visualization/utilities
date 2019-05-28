#!/usr/bin/env python
"""
Calculate statistical results for FITS images
"""

import re
import argparse
import logging
import textwrap
import os.path
from astropy.io import fits
from astropy import stats
import numpy as np
import imutils as iu


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            Calculate statistical quantities for image")
                                    '''),
        epilog=textwrap.dedent('''\
                               '''))
    parser.add_argument("fitsfile", nargs='+',
                        metavar="file", help="input fits file(s)")
    parser.add_argument("--quicklook", action='store_true',
                        help="estimate signal, noise, counts/sec in adus")
    sgroup = parser.add_argument_group("stats", "select statistics and regions"
                                       " (exclusive of quicklook)")
    sgroup.add_argument("--region", nargs='+', metavar='reg',
                        help="region fmt: \"x1:x2,y1:y2\",")
    sgroup.add_argument("--stats", nargs='+', metavar='stat',
                        help="select: mean median stddev min max")
    sgroup.add_argument("--bias", nargs='?', metavar='cols', const='overscan',
                        help="subtract bias, fmt: \"x1:x2\"")
    sgroup.add_argument("--btype", choices=['mean', 'median', 'byrow'],
                        help="bias subtract by-row (def) or constant",
                        default='byrow')
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument("--hduname", nargs='+',
                        metavar='idn', help="process HDU list by names")
    hgroup.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="process HDU list by ids")
    parser.add_argument("--tearing", action='store_true',
                        help="add tearing metric to quicklook output")
    parser.add_argument("--dipoles", action='store_true',
                        help="add dipole metric to quicklook output")
    parser.add_argument("--threshold", nargs=1, metavar='thresh',
                        help="count number of pixels above threshold")
    parser.add_argument("--info", action='store_true',
                        help="print the info() table summarizing file")
    parser.add_argument("--debug", action='store_true',
                        help="print additional debugging messages")
    parser.add_argument("--noheadings", action='store_true',
                        default=False,
                        help="Don't print column heads for stats")
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
        if optlist.info:  # just print the image info and exit
            hdulist.info()
            continue
        if not optlist.noheadings:  # print filename
            print("#")
            print("# {}".format(os.path.basename(ffile)))
        # Construct a list of the HDU's to work on
        hduids = iu.get_hduids(optlist, hdulist)
        if optlist.quicklook:
            quicklook(optlist, hduids, hdulist)
        else:
            if optlist.bias:
                iu.subtract_bias(optlist, hduids, hdulist)
            stats_proc(optlist, hduids, hdulist)
        ncalls.counter = 0  # reset per file, triggers headers


def stats_proc(optlist, hduids, hdulist):
    """print statistics for region according to options
    """
    # Process each HDU in the list "hduids"
    for hduid in hduids:
        pix = hdulist[hduid].data
        name = hdulist[hduid].name
        if optlist.region:
            for reg in optlist.region:  # if there are regions
                logging.debug('processing %s', reg)
                slice_spec = iu.parse_region(reg)
                if slice_spec:
                    stats_print(optlist, hduid, name,
                                pix[slice_spec], reg)
                else:
                    logging.error('skipping')
        else:
            stats_print(optlist, hduid, name, pix, None)


def stats_print(optlist, sid, name, buf, reg):
    """perform and print the given statistics quantities
    """
    if not optlist.stats:
        optlist.stats = ["mean", "median", "stddev", "min", "max"]
    if not optlist.noheadings and ncalls.counter == 0:
        print("#{:>3s} {:>9s}".format("id", "HDUname"), end="")
        if "mean" in optlist.stats:
            print(" {:>8s}".format("mean"), end="")
        if "median" in optlist.stats:
            print(" {:>8s}".format("median"), end="")
        if "stddev" in optlist.stats:
            print(" {:>7s}".format("stddev"), end="")
        if "min" in optlist.stats:
            print(" {:>7s}".format("min"), end="")
        if "max" in optlist.stats:
            print(" {:>7s}".format("max"), end="")
        if reg:
            print(" {:20s}".format("region"), end="")
        print("")  # newline)

    if not optlist.noheadings:
        print(" {:3d} {:>9s}".format(sid, name), end="")

    if "mean" in optlist.stats:
        print(" {:>8g}".format(np.mean(buf)), end="")
    if "median" in optlist.stats:
        print(" {:>8g}".format(np.median(buf)), end="")
    if "stddev" in optlist.stats:
        print(" {:>7.2g}".format(np.std(buf)), end="")
    if "min" in optlist.stats:
        print(" {:>7g}".format(np.min(buf)), end="")
    if "max" in optlist.stats:
        print(" {:>7g}".format(np.max(buf)), end="")
    if reg:
        reg = re.sub(r"^\[*([^\]]*)\]*$", r"\1", reg)
        print(" {:20s}".format(reg), end="")
    print("")  # newline)
    ncalls()  # track call count, acts like static variable)


def quicklook(optlist, hduids, hdulist):
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

    for hduid in hduids:
        # get header details to extract signal, bias regions)
        pix = hdulist[hduid].data
        logging.debug('shape(pix)=%s', np.shape(pix))
        name = hdulist[hduid].name
        hdr = hdulist[hduid].header
        try:
            dstr = hdr['DATASEC']
        except KeyError as ke:
            logging.error('KeyError: %s, required for quicklook mode', ke)
            exit(1)
        logging.debug('DATASEC=%s', dstr)
        res = re.match(r"\[?([0-9]*):([0-9]+),([0-9]+):([0-9]+)\]?",
                       dstr)
        if res:
            datasec = res.groups()
        else:
            emsg = "DATASEC:{} parsing failed".format(dstr)
            logging.error(emsg)
            exit(1)
        naxis1 = int(hdr['NAXIS1'])
        naxis2 = int(hdr['NAXIS2'])
        # define region to measure signal level)
        x1 = int(datasec[0]) - 1
        x2 = int(datasec[1])
        y1 = int(datasec[2]) - 1
        y2 = int(datasec[3])
        sig_buf = pix[y1:y2, x1:x2]
        logging.debug('sig_buf = pix[%s:%s, %s:%s]', y1, y2, x1, x2)
        logging.debug('shape(sig_buf)=%s', np.shape(sig_buf))

        # Serial bias_region = "[y1:y2,x1:x2]"
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
        logging.debug('bias_buf = pix[%s:%s, %s:%s]', y1, y2, x1, x2)
        logging.debug('shape(bias_buf)=%s', np.shape(bias_buf))

        # parallel bias_region = "[y1:y2,x1:x2]"
        x1 = int(datasec[0]) - 1
        x2 = int(datasec[1])
        y1 = int(datasec[3])
        y2 = naxis2
        if y1 > y2 or x1 > x2:
            logging.error('No bias region available for datasec=%s',
                          datasec)
            logging.error(' with naxis1=%s, naxis2=%s', naxis1, naxis2)
            exit(1)
        p_bias_buf = pix[y1:y2, x1:x2]
        logging.debug('p_bias_buf = pix[%s:%s, %s:%s]', y1, y2, x1, x2)
        logging.debug('shape(p_bias_buf)=%s', np.shape(p_bias_buf))
        quicklook_print(optlist, hduid, name, pix, sig_buf,
                        bias_buf, p_bias_buf, expt)


def quicklook_print(optlist, sid, name, pix, sig_buf,
                    bias_buf, p_bias_buf, expt):
    """perform and print the given statistics quantities
       fields are: mean, bias, signal, noise, adu/s
    """
    quick_fields = ["mean", "bias", "signal",
                    "noise", "adu/sec", "eper:s-cte", "eper:p-cte"]
    if optlist.tearing:
        quick_fields.append("tearing")
    if optlist.dipoles:
        quick_fields.append("dipoles")
    if optlist.threshold:
        quick_fields.append("threshold")

    if not optlist.noheadings and ncalls.counter == 0:
        print("#{:>3s} {:>9s}".format("id", "HDUname"), end="")
        if "mean" in quick_fields:
            print(" {:>6s}".format("median"), end="")
        if "bias" in quick_fields:
            print(" {:>5s}".format("bias"), end="")
        if "signal" in quick_fields:
            print(" {:>6s}".format("signal"), end="")
        if "noise" in quick_fields:
            print(" {:>7s}".format("noise"), end="")
        if "adu/sec" in quick_fields and expt > 0:
            print("{:>8s}".format("adu/sec"), end="")
        if "eper:s-cte" in quick_fields:
            print("{:>9s}".format("s-cte"), end="")
        if "eper:p-cte" in quick_fields:
            print("{:>9s}".format("p-cte"), end="")
        if "tearing" in quick_fields:
            print("{:>15s}".format("tearing: L  R"), end="")
        if "dipoles" in quick_fields:
            print("{:>9s}".format("%dipoles"), end="")
        if "threshold" in quick_fields:
            print("{:>9s}".format("N>thresh"), end="")
        print("")  # newline)

    if not optlist.noheadings:
        print(" {:3d} {:>9s}".format(sid, name), end="")

    if "mean" in quick_fields:
        sig_mean = np.median(sig_buf)
        print(" {:>6g}".format(sig_mean), end="")
    if "bias" in quick_fields:
        bias_mean = np.median(bias_buf)
        print(" {:>5g}".format(bias_mean), end="")
    if "signal" in quick_fields:
        signal = sig_mean - bias_mean
        print(" {:>6g}".format(signal), end="")
    if "noise" in quick_fields:
        (nrows, ncols) = np.shape(bias_buf)
        print(" {:>7.4g}".format(
            np.std(bias_buf[int(nrows/2-nrows/20):int(nrows/2+nrows/20),
                            3:ncols-2])), end="")
    if "adu/sec" in quick_fields and expt > 0:
        print("{:>8.2f}".format(float(signal)/expt), end="")
    if "eper:s-cte" in quick_fields:
        logging.debug('s-cte------------------')
        (nrows, ncols) = np.shape(sig_buf)
        nsig_cols = ncols
        # define region to measure signal used in cte calc)
        y1 = int(nrows/2-nrows/10)
        y2 = int(nrows/2+nrows/10)
        x0 = ncols-int(ncols/20)
        x1 = ncols-1
        logging.debug('s_n=sig_buf[%s:%s,%s:%s]', y1, y2, x0, x1)
        s_n = sig_buf[y1:y2, x0:x1]
        # define region to measure bias
        (nrows, ncols) = np.shape(bias_buf)
        l_ncols = 3
        bias_mean = np.mean(bias_buf[y1:y2, l_ncols:ncols])
        logging.debug('using bias_buf[%d:%d,%d:%d]',
                      y1, y2, l_ncols, ncols)
        l_n = np.mean(s_n) - bias_mean
        logging.debug('l_n=%10.6g', l_n)
        y1 = int(nrows/2-nrows/10)
        y2 = int(nrows/2+nrows/10)
        x0 = 0
        x1 = l_ncols
        logging.debug('b_n=bias_buf[%s:%s,%s:%s]', y1, y2, x0, x1)
        b_n = bias_buf[y1:y2, x0:x1]
        l_nn = np.mean((b_n - bias_mean).sum(axis=1))
        logging.debug('l_nn=%10.6g', l_nn)
        if l_n > 0.0:
            eper = 1 - (l_nn / (nsig_cols * l_n))
            print(" {:>8.6g}".format(eper), end="")
    if "eper:p-cte" in quick_fields:
        logging.debug('p-cte------------------')
        (nrows, ncols) = np.shape(p_bias_buf)
        l_nrows = 3  # number of overscan rows used to determing cte
        # define region to measure bias used in cte calc)
        x1 = int(ncols/2-ncols/10)
        x2 = int(ncols/2+ncols/10)
        y0 = l_nrows
        y1 = nrows
        p_bias_mean = np.mean(p_bias_buf[y0:y1, x1:x2])
        logging.debug('p_bias_mean=mean(p_bias_buf[%s:%s, %s:%s])',
                      y0, y1, x1, x2)
        logging.debug('p_bias_mean=%10.6g', p_bias_mean)
        # define region to measure signal used in cte calc)
        (nrows, ncols) = np.shape(sig_buf)
        x1 = int(ncols/2-ncols/10)
        x2 = int(ncols/2+ncols/10)
        y0 = nrows-100
        y1 = nrows-1
        logging.debug('s_n=sig_buf[%s:%s,%s:%s]', y0, y1, x1, x2)
        s_n = sig_buf[y0:y1, x1:x2]
        l_n = np.mean(s_n) - p_bias_mean
        logging.debug('l_n=%10.6g', l_n)
        x1 = int(ncols/2-ncols/10)
        x2 = int(ncols/2+ncols/10)
        y0 = 0
        y1 = l_nrows
        logging.debug('shape(p_bias_buf)=%s', np.shape(p_bias_buf))
        logging.debug('b_n=p_bias_buf[%s:%s,%s:%s]', x1, x2, y0, y1)
        b_n = p_bias_buf[y0:y1, x1:x2]
        logging.debug('shape(b_n)=%s', np.shape(b_n))
        # l_nn = np.median(b_n) - p_bias_mean
        l_nn = np.mean((b_n - bias_mean).sum(axis=0))
        logging.debug('l_nn=%10.6g', l_nn)
        (nrows, ncols) = np.shape(sig_buf)
        if l_n > 0.0:
            eper = 1 - (l_nn / (nrows * l_n))
            print(" {:>8.6g}".format(eper), end="")
    # ---------
    if "tearing" in quick_fields:
        logging.debug('tearing check----------')
    # column-1 into an array arr1)
    # column-2..N into 2nd array with dimension N-1 x Ncols arr2)
    # take median of 2nd array to make 1-D: arr3)
    # find stddev of arr3)
    # form (arr3 - arr1)/stddev as arr4)
    # find the first value of index "j" in sorted(arr4) less than 1.0)
    # report out (len(arr4)-j)/len(arr4) to 1 digit as tearing where)
        # this represents the fraction of the array less than 1.0)
        # left side)
        arr3 = np.median(sig_buf[:, 2:40], axis=1)
        arr4 = (arr3 - sig_buf[:, 0])/np.std(arr3)
        tm = (1.0*np.size(arr4) - np.searchsorted(arr4, 1.0))/np.size(arr4)
        print("{:>4.1f}".format(tm), end="")
        # right side)
        arr3 = np.median(sig_buf[:, -40:-2], axis=1)
        arr4 = (arr3 - sig_buf[:, -0])/np.std(arr3)
        tm = (1.0*np.size(arr4) - np.searchsorted(arr4, 1.0))/np.size(arr4)
        print("{:>4.1f}".format(tm), end="")
    # ---------
    if "dipoles" in quick_fields:
        logging.debug('dipoles check----------')
        # region to work on is sig_buf, say 200 rows near top)
        # transpose to column order)
        # find sigma-clipped mean, median and stdev)
        # subtract mean from array)
        # divide the array by sigma)
        # go through array finding pixel pairs of differing sign)
        # and where |A(n)-A(n+1)| > 6)
        # add one to counter each time such a pair is found)
        # print out the % of pixels occupied by dipoles)
        (nrows, ncols) = np.shape(sig_buf)
        arr1 = sig_buf[-int(nrows/10):-1, :]  # use top 10% of array)
        logging.debug('using subarray [%s:%s,:]', -int(nrows/10), -1)
        arr2 = arr1.flatten('F')  # flatten to 1d in column order)
        avg2, med2, std2 = stats.sigma_clipped_stats(arr2)
        logging.debug('clipped stats: avg:%.3g med:%s stdev:%.3g',
                      avg2, med2, std2)
        arr3 = (arr2 - avg2)/std2
        ndipole = 0
        for i in range(0, np.size(arr3) - 1):
            if (np.sign(arr3[i+1] * arr3[i]) == -1)\
                    and abs(arr3[i+1] - arr3[i]) > 5:
                ndipole += 1
        logging.debug('dipole count = %s', ndipole)
        print("{:>9.2f}".format(
            100.0*float(2*ndipole)/(np.size(arr1))), end="")
    print("")  # newline)
    ncalls()  # track call count, acts like static variable)
    # ---------
    if "threshold" in quick_fields:
        logging.debug('threshold check----------')
        # region to work on is sig_buf
        print("{:>9d}".format(
            np.size(np.where(np.reshape(pix, -1) >= optlist.threshold))),
              end="")


def ncalls():
    """maintain a counter
    """
    ncalls.counter += 1


if __name__ == '__main__':
    main()
