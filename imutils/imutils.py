"""
Common imutil functions
"""
import re
import logging
from astropy.io import fits
import numpy as np


def parse_region(reg):
    """
    Return a pair of slices (slice1, slice2) corresponding
    to the region give as input in ~IRAF format
    If the region can't be parsed None is returned
    """
    # peel off any outer brackets
    reg = re.sub(r"^\[*([^\]]*)\]*$", r"\1", reg)
    #
    # reg = [x1:x2,y1:y2] -- standard rectangle)
    if re.match(r"([0-9]*):([0-9]+),\s*([0-9]+):([0-9]+)$", reg):
        (x1, x2, y1, y2) = re.match(
            r"([0-9]*):([0-9]+),\s*([0-9]+):([0-9]+)$", reg).groups()
        retval = (slice(int(y1)-1, int(y2)), slice(int(x1)-1, int(x2)))
    #
    # reg = [x0,y1:y2] -- single column section)
    elif re.match(r"([0-9]+),\s*([0-9]+):([0-9]+)$", reg):
        (x0, y1, y2) = re.match(
            r"([0-9]+),\s*([0-9]+):([0-9]+)$", reg).groups()
        retval = (slice(int(y1)-1, int(y2)), slice(int(x0)-1))
    #
    # reg = [*,y1:y2]) -- row selection
    elif re.match(r"(\*),\s*([0-9]+):([0-9]+)$", reg):
        (x, y1, y2) = re.match(r"(\*),\s*([0-9]+):([0-9]+)$", reg).groups()
        retval = (slice(int(y1)-1, int(y2)), slice(None, None))
    #
    # reg = [x0,y0] -- single pixel
    elif re.match(r"([0-9]+),\s*([0-9]+)$", reg):
        (x0, y0) = re.match(r"([0-9]+),\s*([0-9]+)$", reg).groups()
        retval = (slice(int(y0)), slice(int(x0)))
    #
    # reg = [x1:x2,y0] -- single row section
    elif re.match(r"([0-9]+):([0-9]+),\s*([0-9]+)$", reg):
        (x1, x2, y0) = re.match(
            r"([0-9]+):([0-9]+),\s*([0-9]+)$", reg).groups()
        retval = (slice(int(y0)-1), slice(int(x1)-1, int(x2)))
    #
    # reg = [x1:x2,*] -- column selection
    elif re.match(r"([0-9]+):([0-9]+),\s*(\*)$", reg):
        (x1, x2, y) = re.match(r"([0-9]+):([0-9]+),\s*(\*)$", reg).groups()
        retval = (slice(None, None), slice(int(x1)-1, int(x2)))
    #
    # reg = [*,*] # redundant, for completeness)
    elif re.match(r"(\*),\s*(\*)$", reg):
        (x, y) = re.match(r"(\*),\s*(\*)$", reg).groups()
        retval = (slice(None, None), slice(None, None))
    #
    # no match found, bad spec
    else:
        logging.error('bad region spec: \'%s\' no match produced', reg)
        retval = None
    #
    return retval


def get_hduids(optlist, hdulist):
    """ return a list of hduids requested in optlist or all by default
    """
    # Construct a list, "hduids" of the HDU's to work on
    hduids = []
    if optlist.hduname:
        for hduname in optlist.hduname:
            try:
                hduids.append(hdulist.index_of(hduname))
            except KeyError as ke:
                emsg = "KeyError: {}".format(ke)
                logging.error(emsg)
                return None

    elif optlist.hduindex:
        hduids = optlist.hduindex
    else:  # all segments with pixel data
        for hdu in hdulist:
            if isinstance(hdu, fits.PrimaryHDU):
                if np.shape(hdu.data):
                    logging.debug('adding %s with index %d to hduid list',
                                  hdu.name, hdulist.index_of(hdu.name))
                    hduids.append(hdulist.index_of(hdu.name))
            elif isinstance(hdu, (fits.ImageHDU, fits.CompImageHDU)):
                logging.debug('adding %s with index %d to hduid list',
                              hdu.name, hdulist.index_of(hdu.name))
                hduids.append(hdulist.index_of(hdu.name))
            else:
                logging.debug('%s with index %d is not type (Comp)ImageHDU',
                              hdu.name, hdulist.index_of(hdu.name))
    return hduids
