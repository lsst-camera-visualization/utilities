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
import os.path
from astropy.io import fits
from astropy import wcs
import numpy as np
import imutils as iu

# module scope objects
operandtypes = {'error': 0, 'file': 1, 'number': 2, 'list': 3}
optypes = ['+', '-', '*', '/']


def parse_args():
    """handle command line"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        Simple Image Arithmetic
                                    '''),
        epilog=textwrap.dedent('''\
        Output fits_image is overwritten if it already exists.
        First image is primary providing output header to be updated.
        At least one operand must be a fits_image. The other operand
        can be a fits image, a single scalar, or a quoted set of scalars
        as in "1.34 1.29" where the number of values must match the desired
        number of requested output data HDUs.
        The "--bias" option provides bias subtraction and is applied to
        all images.
        Note a "--" is often needed to indicate the end of options.
                               '''))
    parser.add_argument("operand1", help="fits_image or number(s) string")
    parser.add_argument("op", help="operator: +,-,*,/",
                        choices=['+', '-', '*', '/'])
    parser.add_argument("operand2", help="fits_image or number(s) string")
    parser.add_argument("result", nargs='?',
                        help="output fits_image")
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument("--hduname", nargs='+',
                        metavar='idn', help="process HDU list by names")
    hgroup.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="process HDU list by ids")
    parser.add_argument("--region", help="region fmt: \"x1:x2, y1:y2\"")
    parser.add_argument("--bias", nargs='?', metavar='cols', const='overscan',
                        help="subtract bias, fmt: \"x1:x2\"")
    parser.add_argument("--btype", default='byrow',
                        choices=['mean', 'median', 'byrow', 'byrowcol'],
                        help="bias subtraction method, default is byrow")
    parser.add_argument("--debug", action='store_true',
                        help="print additional debugging messages")
    return parser.parse_args()


def main():
    """main logic:"""
    optlist = parse_args()
    iu.init_logging(optlist.debug)
    iu.init_warnings()

    # evaluate operands as either a filename, float, floats or error
    op1_type = get_operand_type(optlist.operand1)
    op2_type = get_operand_type(optlist.operand2)
    if op1_type == operandtypes['error'] or op2_type == operandtypes['error']:
        logging.error("operands must be file or float(s) (string)")
        exit(1)
    if not optlist.result:
        logging.error("Output file name required")
        exit(1)

    # arrange operands so file is first
    if (op1_type == operandtypes['number'] and
            op2_type == operandtypes['file']):
        optlist.operand1, optlist.operand2 = optlist.operand2, optlist.operand1

    # Open files, throws exception on error
    hdulist1 = fits.open(optlist.operand1, mode='readonly')
    if op2_type == operandtypes['file']:
        hdulist2 = fits.open(optlist.operand2, mode='readonly')
    else:
        hdulist2 = None
        operand2 = optlist.operand2.split()  # list of floats as strings

    # create output image with primary header and updates
    hdulisto = iu.create_output_hdulist(hdulist1, sys.argv)

    # loop over HDU id's from Master file, copy non-image HDUs
    # and process the image HDUs accordingly
    hduids = iu.get_requested_image_hduids(optlist, hdulist1)
    if hduids is None:
        logging.info('No data HDUs found or requested')
        exit(1)
    for hduid in hduids:  # process these images
        #
        hdu2id = hdulist2.index_of(hdulist1[hduid].name)
        if not hdu2id:
            logging.error('HDU %s does not exist in %s',
                          hdulist1[hduid].name, hdulist2.filename())
        #
        if hdulist2 and np.shape(hdulist1[hduid].data) != \
                np.shape(hdulist2[hdu2id].data):
            logging.error("Images are not comensurate")
            exit(1)
        #
        # prepare the output hdu
        if not optlist.region:
            hduoid = iu.init_hdu(hdulist1[hduid], hdulisto)
            reg = None
        else:
            reg = iu.parse_region(optlist.region)
            if not reg[0] or not reg[1]:
                logging.error('parse_region(%s) failed', optlist.region)
                exit(1)
            hduoid = iu.init_reg(hdulist1[hduid], hdulisto, reg)
        #
        # optionally subtract bias
        if optlist.bias:
            iu.subtract_bias(optlist.bias, optlist.btype, hdulist1[hduid])
            if hdulist2:
                iu.subtract_bias(optlist.bias, optlist.btype, hdulist2[hdu2id])
        #
        # do the arithmetic
        if hdulist2:
            hdulisto[hduoid].data = ffcalc(hdulist1[hduid].data,
                                           hdulist2[hdu2id].data,
                                           optlist.op, reg)
        else:  # scalar or list of scalars
            if len(operand2) == 1:
                arg2 = float(operand2[0])
            else:
                arg2 = float(operand2.pop(0))
            hdulisto[hduoid].data = fscalc(hdulist1[hduid].data,
                                           arg2, optlist.op, reg)
        # finish up this hdu
        hdulisto[hduoid].update_header()
        dtstr = datetime.datetime.utcnow().isoformat(
            timespec='milliseconds')
        hdulisto[hduoid].add_checksum(dtstr)

    for hdu in hdulist1:
        # append to output if it does not contain image data
        if not isinstance(hdu, (fits.ImageHDU,
                                fits.CompImageHDU, fits.PrimaryHDU)):
            hdulisto.append(hdu)

    # write the output file
    hdulisto.writeto(optlist.result, overwrite=True)
    exit(0)


def ffcalc(arr1, arr2, op, reg):
    """
    returns (arr1 op arr2)
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
        '+': faddition,
        '-': fsubtract,
        '*': fmultiply,
        '/': fdivision}
    func = fmap[op]
    arr3 = func(pix1, pix2)
    return arr3.astype('float32')


def fscalc(arr1, s, op, reg):
    """
    returns (array op scalar)
    """
    if reg:
        yslice = reg[0]
        xslice = reg[1]
        pix1 = arr1[yslice, xslice].astype('float64')
    else:
        pix1 = arr1.astype('float64')
    fmap = {
        '+': faddition,
        '-': fsubtract,
        '*': fmultiply,
        '/': fdivision}
    func = fmap[op]
    arr3 = func(pix1, s)
    return arr3.astype('float32')


def faddition(arg1, arg2):
    """add the args"""
    return arg1 + arg2


def fsubtract(arg1, arg2):
    """subtract the args"""
    return arg1 - arg2


def fmultiply(arg1, arg2):
    """multiply the args"""
    return arg1 * arg2


def fdivision(arg1, arg2):
    """divide arg1 by arg2 """
    return arg1 / arg2


def get_operand_type(operand):
    """validates operand string is either a filename or float(str)
       returns 0=file, 1=number, 2=list, -1 otherwise
    """
    logging.debug("operand = {}".format(operand))
    # Check for file first
    if os.path.isfile(operand):
        return operandtypes['file']
    # Check for scalars (single or list)
    values = operand.split()
    for val in values:
        logging.debug("scalar = {}".format(val))
        try:
            float(val)
        except ValueError as verr:
            emsg = "ValueError: {}".format(verr)
            logging.error(emsg)
            return operandtypes['error']
    return operandtypes['list']


if __name__ == '__main__':
    main()
