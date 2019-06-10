#!/usr/bin/env python
"""
Compute noise power spectrum for fits file images
"""

import re
import argparse
import logging
import textwrap
import numpy as np
from scipy import signal
from astropy.io import fits
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
import imutils as iu


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
                                    Compute Image Noise Power Spectrum
                                    '''),
        epilog=textwrap.dedent('''\
                               '''))
    parser.add_argument("fitsfile", help="fits_image",
                        metavar="file", nargs='+',)
    parser.add_argument("--rt", help="pixel read time (ns)",
                        type=int, default=1800)
    parser.add_argument("--row", help="choose starting image row",
                        default=10, type=int)
    parser.add_argument("--nrows", help="number of rows to average",
                        default=1, type=int)
    # hdu name|index exclusive
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument("--hduname", nargs='+',
                        metavar='idn', help="process HDU list by names")
    hgroup.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="process HDU list by ids")
    parser.add_argument("--scaling", choices=['density', 'spectrum'],
                        default='density', help="")
    parser.add_argument("--clip", action='store_true',
                        help="apply sigma (3) clip to rows used")
    parser.add_argument("--info", action='store_true',
                        help="print the info() table summarizing file")
    parser.add_argument("--debug", action='store_true',
                        help="print additional debugging messages")
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
    if optlist.scaling == 'density':
        window = 'boxcar'
    else:
        window = 'flattop'

    # Open files
    fileno = 0
    hduids = []
    # loop over files
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
        # Construct a list of the HDU's to work on
        hduids = iu.get_requested_image_hduids(optlist, hdulist)

        # loop over hdu's
        hducnt = 0
        for hduid in hduids:
            hdr = hdulist[hduid].header
            try:
                dstr = hdr['DATASEC']
            except KeyError as ke:
                emsg = "KeyError: {}, required".format(ke)
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

            # define region to measure
            x1 = int(datasec[0]) - 1
            x2 = int(datasec[1])
            # y1 = int(datasec[2]) - 1
            # y2 = int(datasec[3])

            # naxis1 = int(hdr['NAXIS1'])
            # naxis2 = int(hdr['NAXIS2'])
            stddev = float(hdr['STDVBIAS'])
            pix = hdulist[hduid].data
            fs = 1.0/(optlist.rt*1e-9)
            # measure the size needed
            arr = pix[optlist.row, x1:x2]
            x, p = signal.periodogram(arr, fs, window,
                                      scaling=optlist.scaling)
            flen = x.size
            plen = p.size
            if flen != plen:
                emsg = "flen({}) != plen({})".format(flen, plen)
                emsg = "DATASEC:{} parsing failed".format(dstr)
                logging.error(emsg)
                exit(1)
            # now do the real calculation
            f = np.empty((optlist.nrows, flen))
            Pxx_den = np.empty((optlist.nrows, plen))
            for rr in range(0, optlist.nrows):
                arr = pix[rr + optlist.row, x1:x2]
                if optlist.clip:
                    amed = np.median(arr)
                    farr = sigma_clip(arr)
                    x, p = signal.periodogram(farr.filled(amed),
                                              fs, window,
                                              scaling=optlist.scaling)
                else:
                    x, p = signal.periodogram(arr, fs, window,
                                              scaling=optlist.scaling)
                f[rr] = x
                Pxx_den[rr] = p

            f_avg = np.average(f, axis=0)
            Pxx_den_avg = np.average(Pxx_den, axis=0)
            if (fileno + hducnt) == 0:
                pmin = Pxx_den_avg.min()
                pmax = Pxx_den_avg.max()
                pavg = 0.5*(pmin + pmax)
                debugmsg = "pmin0={:>g}".format(pmin)
                logging.debug(debugmsg)
                debugmsg = "pmax0={:>g}".format(pmax)
                logging.debug(debugmsg)
            else:
                if pmin > Pxx_den_avg.min():
                    pmin = Pxx_den_avg.min()
                    debugmsg = "pmin={:>g}".format(pmin)
                    logging.debug(debugmsg)
                if pmax < Pxx_den_avg.max():
                    pmax = Pxx_den_avg.max()
                    debugmsg = "pmax={:>g}".format(pmax)
                    logging.debug(debugmsg)

            plt.semilogy(f_avg, Pxx_den_avg,
                         label="{}:{:>02d}:{:>7.2f}".format(
                             fileno, hduid, stddev))
            hducnt += 1
            # end loop over hdui's
        fileno += 1
        # end loop over files
    #
    plt.ylim([0.8*pmin, 1.2*pmax])
    plt.xlabel('freqquency [Hz]')
    if optlist.scaling == 'density':
        plt.ylabel('PSD [V**2/Hz]')
    else:
        plt.ylabel('Linear spectrum [V RMS]')
    plt.grid(True)
    plt.legend(fontsize='xx-small', title='File:HDUi RN')
    plt.show()


if __name__ == '__main__':
    main()
