#!/usr/bin/env python
"""
Perform image arithmetic
"""

import sys
import re
import argparse
import logging
import textwrap
import datetime
from string import join
import os.path
from astropy.io import fits
from astropy import wcs
import numpy as np

#- module scope objects
operandtypes = {'error': 0, 'file': 1, 'number': 2}
optypes = ['+', '-', '*', '/']


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        Simple Image Arithmetic
                                    '''),
        epilog=textwrap.dedent('''\
        Output fits_image is created (can't already exist)
        In case of two images, first is primary with header copied
        to output with comments added
                               '''))
    parser.add_argument("operand1", help="fits_image or number")
    parser.add_argument("op", help="operator: +,-,*,/",
                        choices=['+', '-', '*', '/'])
    parser.add_argument("operand2", help="fits_image or number")
    parser.add_argument("result", nargs='?',
                        help="output fits_image")
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument("--hduname", nargs='+',
                        metavar='idn', help="process HDU list by names")
    hgroup.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="process HDU list by ids")
    parser.add_argument("--region", help="region fmt: \"x1:x2,y1:y2\"")
    return parser.parse_args()


def main():
    """main logic:"""
    logging.basicConfig(format='%(levelname)s: %(message)s',
                        level=logging.DEBUG)
    optlist = parse_args()
    # evaluate operands as either a filename, float or error
    op1_type = get_operand_type(optlist.operand1)
    if op1_type == operandtypes['error']:
        logging.error("Invalid Operand1: must be file or number")
        exit(1)
    op2_type = get_operand_type(optlist.operand2)
    if op2_type == operandtypes['error']:
        logging.error("Invalid Operand2: must be file or number")
        exit(1)
        if ((op1_type == operandtypes['file'] or
             op2_type == operandtypes['file']) and
                not optlist.result):
            logging.error("Output file name required")
            exit(1)

    #- 4 cases
    if (op1_type == operandtypes['file'] and
        op2_type == operandtypes['file']):
        # get hdulist1
        # get hdulist2
        # create output hdulist/image -- verify filename etc.
        # do the pixel manipulations
        # load the output headers and update accordingly
        # write the output file
        #----------------------------------
        try:
            hdulist1 = fits.open(optlist.operand1, memmap=True)
        except IOError as ioerr:
            emsg = "IOError: {}".format(ioerr)
            logging.error(emsg)
            exit(1)
        try:
            hdulist2 = fits.open(optlist.operand2, memmap=True)
        except IOError as ioerr:
            emsg = "IOError: {}".format(ioerr)
            logging.error(emsg)
            exit(1)
        # create the output image and append the primary hdu
        hdulisto = fits.HDUList()
        hdulisto.append(hdulist1[0])
        # add/modify keys to document the changes
        #- update DATE keyword and preserve old one
        hdr = hdulisto[0].header
        hdr.rename_keyword('DATE','DATE0')
        cstr = hdr.comments['DATE0']
        hdr.comments['DATE0'] = "Previous file date/time"
        # FITS date format: 'yyyy-mm-ddTHH:MM:SS[.sss]'
        dt = datetime.datetime.utcnow()
        dtstr = ("{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}".
                 format(dt.year, dt.month, dt.day, dt.hour, dt.minute,
                        dt.second, int(dt.microsecond/1e3)))
        hdr.insert('DATE0',('DATE',dtstr,cstr))
        #- add HISTORY lines
        hdr.add_history("Header written by {} at: {}Z".
                        format(os.path.basename(sys.argv[0]), dtstr))
        hdr.add_history("CMD: {} {}".format(
            os.path.basename(sys.argv[0]), join(sys.argv[1:])))
        print "{}".format(hdr.tostring('\n'))
        #- loop over HDUs
        for hduid in range(1, len(hdulist1)):
            if (optlist.hduindex and str(hduid) not in optlist.hduindex):
                continue
            if (optlist.hduname and
                    hdulist1[hduid].name not in optlist.hduname):
                continue
            if not isinstance(hdulist1[hduid], fits.ImageHDU):
                # copy non-image HDUs to output
                hdulisto.append(hdulist1[hduid])
                continue
            # verify operand2 is valid
            if not (isinstance(hdulist2[hduid], fits.ImageHDU) and
                    (hdulist1[hduid].header['NAXIS1'] ==
                     hdulist2[hduid].header['NAXIS1']) and
                    (hdulist1[hduid].header['NAXIS2'] ==
                     hdulist2[hduid].header['NAXIS2'])):
                logging.error("Images are not comensurate")
                exit(1)

            #- create the hdu
            wcses = wcs.find_all_wcs(hdulist1[hduid].header)
            hdulisto.append(fits.ImageHDU(None, None,
            hduoid = len(hdulisto) - 1
            print "hduoid={}".format(hduoid)
            hdr1 = hdulist1[hduid].header,
            #- parse region to create slice
            reg = ""
            if optlist.region:
                (yslice, xslice) = parse_region(optlist.region)
                reg = (yslice, xslice)
                naxis2 = yslice.stop - yslice.start
                naxis1 = xslice.stop - xslice.start
                print "naxis1={} naxis2={}".format(naxis1, naxis2)
                hdr = hdulisto[hduoid].header
                hdr['NAXIS'] = 2
                hdr.set('NAXIS1', naxis1,
                        "size of the n'th axis", after='NAXIS')
                hdr.set('NAXIS2', naxis2,
                        "size of the n'th axis", after='NAXIS1')
                #- remove any wcs info?
                new_wcses = []
                for w in wcses:
                    wreg = w.slice(reg)
                    new_wcses.append(wreg)
                for w in new_wcses:
                    print w.wcs.name
                    print w.to_header_string()
            #- do the arithmetic
            ffcalc(hdulist1[hduid], hdulist1[hduid], optlist.op,
                   hdulisto[hduoid], reg)
            #print "{}".format(hdr.tostring('\n'))
            hdulisto[hduoid].update_header()
            #print "{}".format(hdr.tostring('\n'))

        #- write the output file
        hdulisto.writeto(optlist.result, overwrite=True)
        hdulisto.close()

        #- update checksum?
        exit(0)

    elif (op1_type == operandtypes['file'] and
          op2_type == operandtypes['number']):
        print "file and number"
        # get hdulist1
        #
        # create output hdulist/image -- verify filename etc.
        # do the pixel manipulations
        # load the output headers and update accordingly
        # write the output file
    elif (op1_type == operandtypes['number'] and
          op2_type == operandtypes['file']):
        print "number and file"
        # same as previous case, just swap operands
    elif (op1_type == operandtypes['number'] and
          op2_type == operandtypes['number']):
        # no putput file, just print the number
        print "{}".format(scalc(float(optlist.operand1),
                                float(optlist.operand2), optlist.op))
    else:
        logging.error("Invalid operands")
        exit(1)

    exit(0)


def ffcalc(hdu1, hdu2, op, hduout, reg):
    """
    """
    print "ffcalc() called"
    if reg:
        yslice = reg[0]
        xslice = reg[1]
        pix1 = hdu1.data[yslice, xslice]
        pix2 = hdu2.data[yslice, xslice]
    else:
        pix1 = hdu1.data
        pix2 = hdu2.data
    fmap = {
        '+': (ffaddition, pix1, pix2),
        '-': (ffsubtract, pix1, pix2),
        '*': (ffmultiply, pix1, pix2),
        '/': (ffdivision, pix1, pix2)}
    func, arr1, arr2 = fmap[op]
    hduout.data = func(arr1, arr2)

def ffaddition(arg1, arg2):
    """add the args"""
    return arg1 + arg2

def ffsubtract(arg1, arg2):
    """subtract the args"""
    return arg1 - arg2

def ffmultiply(arg1, arg2):
    """multiply the args"""
    return arg1 * arg2

def ffdivision(arg1, arg2):
    """divide arg1 by arg2 """
    return arg1 / arg2


def scalc(opd1, opd2, op):
    """
    """
    fmap = {
        '+': (saddition, opd1, opd2),
        '-': (ssubtract, opd1, opd2),
        '*': (smultiply, opd1, opd2),
        '/': (sdivision, opd1, opd2)}
    func, arg1, arg2 = fmap[op]
    return func(arg1, arg2)

def saddition(arg1, arg2):
    """add the args"""
    return arg1 + arg2

def ssubtract(arg1, arg2):
    """subtract the args"""
    return arg1 - arg2

def smultiply(arg1, arg2):
    """multiply the args"""
    return arg1 * arg2

def sdivision(arg1, arg2):
    """divide arg1 by arg2 """
    return arg1 / arg2

def get_operand_type(opstring):
    """validates operand string is either a filename or float(str)
       returns 0=file, 1=number, -1 otherwise
    """
    if os.path.isfile(opstring):
        return operandtypes['file']
    try:
        float(opstring)
    except ValueError as verr:
        emsg = "ValueError: {}".format(verr)
        logging.error(emsg)
        return operandtypes['error']
    return operandtypes['number']

def parse_region(reg):
    """return list with 2 slices for the region
    """
    #- region = [x1:x2,y1,y2]
    res = re.match(r"\[*([0-9]*):([0-9]+),([0-9]+):([0-9]+)\]*", reg)
    if res:
        (x1, x2, y1, y2) = res.groups()
        return (slice(int(y1),int(y2)), slice(int(x1), int(x2)))
    #- region = [*,y1:y2]
    res = re.match(r"\[*(\*),([0-9]+):([0-9]+)\]*", reg)
    if res:
        (x, y1, y2) = res.groups()
        return (slice(int(y1),int(y2)), slice(None, None))
    # reg = [x1:x2,*]
    res = re.match(r"\[*([0-9]+):([0-9]+),(\*)\]*", reg)
    if res:
        (x1, x2, y) = res.groups()
        return (slice(None, None), slice(int(x1), int(x2)))
    # reg = [*,*] #- redundant, but for completeness
    res = re.match(r"\[*(\*),(\*)\]*", reg)
    if res:
        (x, y) = res.groups()
        return (slice(None, None), slice(None, None))
    else:
        emsg = "bad region spec {}, no match produced".\
            format(reg)
        logging.error(emsg)
        exit(1)




if __name__ == '__main__':
    main()

