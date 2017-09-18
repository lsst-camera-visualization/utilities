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
        to output with comments added and adjustments
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
        #- Open file used as the Master being modified
        try:
            hdulist1 = fits.open(optlist.operand1, memmap=True)
        except IOError as ioerr:
            emsg = "IOError: {}".format(ioerr)
            logging.error(emsg)
            exit(1)
        #- Open file used to modify the Master
        try:
            hdulist2 = fits.open(optlist.operand2, memmap=True)
        except IOError as ioerr:
            emsg = "IOError: {}".format(ioerr)
            logging.error(emsg)
            exit(1)
        # Create the output image and append the primary hdu
        hdulisto = fits.HDUList()
        hdulisto.append(hdulist1[0])
        # add/modify keys to document the changes
        #- update DATE keyword and preserve old one
        #- Note checksums will be invalid
        hdr = hdulisto[0].header
        hdr.rename_keyword('DATE', 'DATE0')
        cstr = hdr.comments['DATE0']
        hdr.comments['DATE0'] = "Previous file date/time"
        # FITS date format: 'yyyy-mm-ddTHH:MM:SS[.sss]'
        dt = datetime.datetime.utcnow()
        dtstr = ("{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}".
                 format(dt.year, dt.month, dt.day, dt.hour, dt.minute,
                        dt.second, int(dt.microsecond/1e3)))
        hdr.insert('DATE0',('DATE', dtstr, cstr))
        #- add HISTORY lines
        hdr.add_history("Header written by {} at: {}Z".
                        format(os.path.basename(sys.argv[0]), dtstr))
        hdr.add_history("CMD: {} {}".format(
            os.path.basename(sys.argv[0]), join(sys.argv[1:])))
        hdulisto[0].add_checksum(dtstr, True)
        #- loop over HDU id's from Master file, copy non-image HDUs
        #- and process the image HDUs accordingly
        for hduid in range(1, len(hdulist1)):
            if (optlist.hduindex and
                    hduid not in optlist.hduindex):
                continue
            if (optlist.hduname and
                    hdulist1[hduid].name not in optlist.hduname):
                continue
            if not isinstance(hdulist1[hduid], fits.ImageHDU):
                # not an Image so just copy to output
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
            reg = ""
            if not optlist.region:
                hdr1 = hdulist1[hduid].header
                hdulisto.append(
                    fits.ImageHDU(None, hdr1, hdr1['EXTNAME']))
                hduoid = len(hdulisto) - 1
                hdr = hdulisto[hduoid].header
                logging.debug("Before Update header:")
                logging.debug("\n{}".format(hdr.tostring('\n')))
                #- Update header
                #- set dimensions
                naxis1 = hdr1['NAXIS1']
                naxis2 = hdr1['NAXIS2']
                hdr['NAXIS'] = 2
                hdr.set('NAXIS1', naxis1,
                        "size of the n'th axis", after='NAXIS')
                hdr.set('NAXIS2', naxis2,
                        "size of the n'th axis", after='NAXIS1')
                hdr['BITPIX'] = -32
                logging.debug("After Update header:")
                logging.debug("\n{}".format(hdr.tostring('\n')))
            else: #- handle regions
                hdulisto.append( #- no header template
                    fits.ImageHDU(None, None, None))
                hduoid = len(hdulisto) - 1
                logging.debug("hduoid={}".format(hduoid))
                hdr1 = hdulist1[hduid].header.copy()
                wcses = wcs.find_all_wcs(hdr1)
                (yslice, xslice) = parse_region(optlist.region)
                reg = (yslice, xslice)
                naxis2 = yslice.stop - yslice.start
                naxis1 = xslice.stop - xslice.start
                hdr = hdulisto[hduoid].header
                logging.debug("Before Update header:")
                logging.debug("\n{}".format(hdr.tostring('\n')))
                hdr['NAXIS'] = 2
                hdr.set('NAXIS1', naxis1,
                        "size of the n'th axis", after='NAXIS')
                hdr.set('NAXIS2', naxis2,
                        "size of the n'th axis", after='NAXIS1')
                hdr['EXTNAME'] = hdr1['EXTNAME']
                if hdr1.count('CHANNEL'):
                    hdr['CHANNEL'] = hdr1['CHANNEL']
                #- update any wcses
                new_wcses = []
                for w in wcses:
                    wreg = w.slice(reg)
                    wreghdr = wreg.to_header()
                    for card in wreghdr.cards:
                        hdr.append(card)
                #- add HISTORY lines
                hdr.add_history("Header written by {}".
                           format(os.path.basename(sys.argv[0])))
                logging.debug("After Update header:")
                logging.debug("\n{}".format(hdr.tostring('\n')))
            #- do the arithmetic
            hdulisto[hduoid].data = ffcalc(hdulist1[hduid].data,
                                           hdulist2[hduid].data,
                                           optlist.op, reg)
            hdulisto[hduoid].update_header()
            dt = datetime.datetime.utcnow()
            dtstr = ("{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}".
                    format(dt.year, dt.month, dt.day, dt.hour, dt.minute,
                            dt.second, int(dt.microsecond/1e3)))
            hdulisto[hduoid].add_checksum(dtstr)

        #- write the output file
        hdulisto.writeto(optlist.result, overwrite=True)
        hdulisto.close()
        hdulist1.close()
        hdulist2.close()
        exit(0)

    elif (op1_type == operandtypes['file'] and
          op2_type == operandtypes['number']):
        print ""
        # get hdulist1
        #
        # create output hdulist/image -- verify filename etc.
        # do the pixel manipulations
        # load the output headers and update accordingly
        # write the output file
    elif (op1_type == operandtypes['number'] and
          op2_type == operandtypes['file']):
        print ""
        # same as previous case, just swap operands
    elif (op1_type == operandtypes['number'] and
          op2_type == operandtypes['number']):
        # no putput file, just print the number
        print "{}".format(scalc(float(optlist.operand1),
                                float(optlist.operand2), optlist.op))
        exit(0)
    else:
        logging.error("Invalid operands")
        exit(1)



def ffcalc(arr1, arr2, op, reg):
    """
    """
    if reg:
        yslice = reg[0]
        xslice = reg[1]
        pix1 = arr1[yslice, xslice].astype('float64')
        pix2 = arr2[yslice, xslice].astype('float64')
    else:
        pix1 = arr1.astype('float64')
        pix2 = arr2.astype('float64')
    fmap = {
        '+': (ffaddition, pix1, pix2),
        '-': (ffsubtract, pix1, pix2),
        '*': (ffmultiply, pix1, pix2),
        '/': (ffdivision, pix1, pix2)}
    func, arr1, arr2 = fmap[op]
    arr3 = func(arr1, arr2)
    return arr3.astype('float32')

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

