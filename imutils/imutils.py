"""
Common imutil functions
"""
import re
import logging
from astropy.io import fits
from astropy import stats
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


def subtract_bias(optlist, hduids, hdulist):
    """
    """
    logging.debug('processing bias subtraction')
    for hduid in hduids:  # process each with regions
        if optlist.bias == 'overscan':
            logging.debug('use overscan subset in hdu:%s',
                          hduid)
            hdr = hdulist[hduid].header
            try:
                dstr = hdr['DATASEC']
            except KeyError as ke:
                logging.error('KeyError: %s required', ke)
                exit(1)
            # logging.debug('DATASEC=%s', dstr)
            try:
                naxis1 = hdr['NAXIS1']
            except KeyError as ke:
                logging.error('KeyError: %s required', ke)
                exit(1)
            bias_spec = parse_region(dstr)
            (x1, x2) = (bias_spec[1].start or 0,
                        bias_spec[1].stop or len(hdulist[hduid].data))
            if naxis1 - x2 > x1:  # area to right of DATASEC
                o_width = naxis1 - x2
                # logging.debug(
                #   'overscan= %d:%d', x2 + 1, naxis1)
                if o_width > 10:
                    (b1, b2) = (x2 + 4, naxis1)
                else:
                    (b1, b2) = (x2, naxis1)
            else:  # area to left of DATASEC
                o_width = x1
                # logging.debug('overscan= 1:%d', x1)
                if o_width > 10:
                    (b1, b2) = (0, x1 - 3)
                else:
                    (b1, b2) = (0, x1)
            reg = "[{}:{},*]".format(b1, b2)
        elif len(optlist.bias) == 1:
            # peel off outer brackets
            reg = re.sub(r"^\[*([^\]]*)\]*$", r"\1", optlist.bias)
            # recast as a column selection region
            reg = "[{},*]".format(reg)
        else:
            logging.error('bad bias column selection %s', optlist.bias)
            exit(1)

        # parse the bias spec and process
        logging.debug('bias column selection: %s', reg)
        bias_spec = parse_region(reg)
        # logging.debug('bias_spec: %s', bias_spec)
        (y1, y2) = (bias_spec[0].start or 0,
                    bias_spec[0].stop or len(hdulist[hduid].data[:, 0]))
        (x1, x2) = (bias_spec[1].start or 0,
                    bias_spec[1].stop or len(hdulist[hduid].data[0, :]))
        #  extract the bias region data
        biasbuf = hdulist[hduid].data[y1:y2, x1:x2]
        #  convert biasbuf into a scalar or a 1d array
        logging.debug('biastype=%s', optlist.btype)
        if not optlist.btype or optlist.btype == 'byrow':
            # get sigma clipped median per row, deal with mask also
            b1_avg, bias, b1_std = stats.sigma_clipped_stats(
                biasbuf, axis=1)
            fbias = bias.data[~bias.mask]
            hdulist[hduid].data = \
                hdulist[hduid].data - fbias[:, None]
        elif optlist.btype == 'mean':
            bias = np.mean(biasbuf)
            hdulist[hduid].data = \
                hdulist[hduid].data - bias
        elif optlist.btype == 'median':
            bias = np.median(biasbuf)
            hdulist[hduid].data = \
                hdulist[hduid].data - bias
        else:
            logging.error('btype: %s not valid', optlist.btype)
            exit(1)
