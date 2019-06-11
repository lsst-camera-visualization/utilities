"""
Common imutil functions
"""
import re
import logging
import warnings
import datetime
import os.path
from astropy.io import fits
from astropy import stats
from astropy import wcs
from astropy.utils.exceptions import AstropyWarning
import numpy as np


def init_logging(debug):
    """
    """
    if debug:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.INFO)


def init_warnings():
    """
    """
    warnings.simplefilter('ignore', category=AstropyWarning)


def create_output_hdulist(hdulisti, argv):
    """
    """
    logging.debug("creating output hdulist")
    # Create the output image, copy and update header comments, history
    hdulisto = fits.HDUList()
    hdulisto.append(fits.ImageHDU(None, hdulisti[0].header))
    hdr = hdulisto[0].header
    cstr = hdr.comments['DATE']  # copy comment
    hdr.rename_keyword('DATE', 'DATEORIG', force=True)
    hdr.comments['DATEORIG'] = "Previous file date/time"
    # FITS date format: 'yyyy-mm-ddTHH:MM:SS[.sss]'
    dtstr = datetime.datetime.utcnow().isoformat(timespec='milliseconds')
    hdr.insert('DATEORIG', ('DATE', dtstr, cstr))
    # add HISTORY lines
    hdr.add_history("Header written by {} at: {}".
                    format(os.path.basename(argv[0]), dtstr))
    hdr.add_history("CMD: {} {}".format(
        os.path.basename(argv[0]), ' '.join(argv[1:])))
    return hdulisto


def init_reg(hdui, hdulisto, region):
    """
    Use hdui as a template to append a new hdu to hdulisto
    copy the header and set the size in preparation defined by
    the region for data to be added later. Only sets sizes for PrimaryHDU.
    """
    # create the output hdu region from the master
    logging.debug('region = {}'.format(region))
    if not isinstance(hdui, fits.PrimaryHDU):
        hdri = hdui.header.copy()
        hdulisto.append(fits.ImageHDU(None, hdri, hdri['EXTNAME']))
    hduoid = len(hdulisto) - 1
    logging.debug('hduoid=%s', hduoid)
    naxis2 = (region[0].stop or len(hdui.data[:, 0])) - (region[0].start or 0)
    naxis1 = (region[1].stop or len(hdui.data[0, :])) - (region[1].start or 0)
    hdro = hdulisto[hduoid].header
    hdro['NAXIS'] = 2
    hdro.set('NAXIS1', naxis1,
             "size of the n'th axis", after='NAXIS')
    hdro.set('NAXIS2', naxis2,
             "size of the n'th axis", after='NAXIS1')
    hdro['EXTNAME'] = hdri['EXTNAME']
    if hdri.count('CHANNEL'):
        hdro['CHANNEL'] = hdri['CHANNEL']
    # update any wcses
    wcses = wcs.find_all_wcs(hdro, fix=False)
    for w in wcses:
        wreg = w.slice(region)
        wreghdro = wreg.to_header()
        for card in wreghdro.cards:
            key = card.keyword
            value = card.value
            comment = card.comment
            hdro.set(key, value, comment)
    return hduoid


def init_hdu(hdui, hdulisto):
    """
    Use hdui as a template to append a new hdu to hdulisto
    copy the header and set the size in preparation for data
    to be added later.  Only sets sizes for PrimaryHDU.
    """
    # create the output hdu from the master
    if not isinstance(hdui, fits.PrimaryHDU):
        hdri = hdui.header.copy()
        hdulisto.append(fits.ImageHDU(None, hdri, hdri['EXTNAME']))
    hduoid = len(hdulisto) - 1
    hdro = hdulisto[hduoid].header
    # logging.debug('output header before:\n%s\n', hdro.tostring())
    hdro['NAXIS'] = 2
    hdro.set('NAXIS1', hdri['NAXIS1'],
             "size of the n'th axis", after='NAXIS')
    hdro.set('NAXIS2', hdri['NAXIS2'],
             "size of the n'th axis", after='NAXIS1')
    hdro['BITPIX'] = -32
    # logging.debug('output header after:\n%s\n', hdro.tostring())
    return hduoid


def parse_region(reg):
    """
    Return a pair of slices (slice1, slice2) corresponding
    to the region give as input in ~IRAF format
    If the region can't be parsed (None, None) is returned
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
        retval = (None, None)
    #
    return retval


def get_requested_image_hduids(optlist, hdulist):
    """
    Return a list of imaage hduids requested in optlist or all by default.
    Check that they exist in hdulist and have data.  Requested hduids that
    don't exist are skipped.
    """
    # Construct a list of candidate HDUs to work on
    chduids = []  # list of candidate hduids
    if optlist.hduname:
        for hduname in optlist.hduname:
            try:
                chduids.append(hdulist.index_of(hduname))
            except KeyError as ke:
                logging.error('KeyError: %s', ke)
                logging.error('HDU[%s] not found, skipping', hduname)
    if optlist.hduindex:
        for hduid in optlist.hduindex:
            try:
                hdu = hdulist[hduid]
                chduids.append(hduid)
            except IndexError:
                logging.error('HDU[%d] not found, skipping', hduid)
    if not optlist.hduindex and not optlist.hduname:
        for hdu in hdulist:
            chduids.append(hdulist.index(hdu))

    # Validate the list of candidate HDUs, keep those with pixels
    hduids = []
    for hduid in chduids:
        hdu = hdulist[hduid]
        if isinstance(hdu, fits.PrimaryHDU):  # check for data
            if np.shape(hdu.data):
                logging.debug('adding %s with index %d to hduid list',
                              hdu.name, hdulist.index_of(hdu.name))
                hduids.append(hduid)
        elif isinstance(hdu, (fits.ImageHDU, fits.CompImageHDU)):
            logging.debug('adding %s with index %d to hduid list',
                          hdu.name, hdulist.index_of(hdu.name))
            hduids.append(hdulist.index_of(hdu.name))
        else:
            logging.debug('%s with index %d is not type (Comp)ImageHDU',
                          hdu.name, hdulist.index_of(hdu.name))
    if hduids:
        return hduids
    return None


def get_overscan_spec(hdu):
    """
    Given an hdu, return serial and parallel overscan regions as slices
    Returns a pair of regions (sbias, pbias) where each are in slice notation
    """
    # first get serial and parallel overscan region defs
    hdr = hdu.header
    try:
        dstr = hdr['DATASEC']
    except KeyError as ke:
        logging.error('KeyError: %s required', ke)
        return (None, None)
    logging.debug('DATASEC=%s', dstr)
    try:
        naxis1 = hdr['NAXIS1']
    except KeyError as ke:
        logging.error('KeyError: %s required', ke)
        return (None, None)
    try:
        naxis2 = hdr['NAXIS2']
    except KeyError as ke:
        logging.error('KeyError: %s required', ke)
        return (None, None)
    # get DATASEC region
    data_spec = parse_region(dstr)  # this can return (None, None)
    (p1, p2) = (data_spec[0].start or 0,
                data_spec[0].stop or len(hdu.data[:, 0]))
    (s1, s2) = (data_spec[1].start or 0,
                data_spec[1].stop or len(hdu.data[0, :]))
    if naxis1 - s2 > s1:  # area to right of DATASEC
        (b1, b2) = (s2, naxis1)
    else:  # area to left of DATASEC
        (b1, b2) = (0, s1)
    soscan = (slice(0, naxis2), slice(b1, b2))
    poscan = (slice(p2, naxis2), slice(0, naxis1))

    return (soscan, poscan)


def subtract_bias(bstring, btype, hdu):
    """
    Subtract a bias (constant or row base or special) from an hdu.
    Choices are mean, median or byrow subtraction of a bias calculated
    in either a given set of columns or using DATASEC to infer the
    overscan region.  Using DATASEC skips the 1st 4 columns of overscan.
    The rows-and-columns version uses both the serial and parallel oscan
    to account for the shape of the bias (particularly for ITL sensors)
    """
    (soscan, poscan) = get_overscan_spec(hdu)
    #  [sp]oscan is ((y1, y2), (x1, x2)) format

    if bstring == 'overscan':
        logging.debug('use overscan subregion in hdu:%s', hdu.name)
    elif len(bstring) >= 3:  # user supplies column selection
        # eg. a:b is minimum spec size
        # peel off outer brackets
        reg = re.sub(r"^\[*([0-9]+:[0-9]+)\]*$", r"\1", bstring)
        (x1, x2) = re.match(r"([0-9]+):([0-9]+)$", reg).groups()
        soscan[1] = (x1, x2)
    else:
        logging.error('bad bias selection %s', bstring)
        exit(1)

    logging.debug('biastype=%s', btype)
    # default method is "byrow"
    if not btype or btype == 'byrow':
        # get sigma clipped median per row
        so_avg, so_bias, so_std = stats.sigma_clipped_stats(
            hdu.data[soscan], axis=1)
        so_bias = so_bias.reshape(np.shape(so_bias)[0], 1)
        logging.debug("np.shape(so_bias)={}".format(np.shape(so_bias)))
        hdu.data = hdu.data - so_bias.data
    elif btype == 'mean':
        bias = np.mean(hdu.data[soscan])
        hdu.data = hdu.data - bias
    elif btype == 'median':
        bias = np.median(hdu.data[soscan])
        hdu.data = hdu.data - bias
    elif btype == 'byrowcol':
        # get sigma clipped median per row
        so_avg, so_bias, so_std = stats.sigma_clipped_stats(
            hdu.data[soscan], axis=1)
        so_bias = so_bias.reshape(np.shape(so_bias)[0], 1)
        logging.debug("np.shape(so_bias)={}".format(np.shape(so_bias)))
        hdu.data = hdu.data - so_bias.data
        # get sigma clipped median per column
        logging.debug('poscan=((%d, %d), (%d, %d))',
                      poscan[0].start, poscan[0].stop,
                      poscan[1].start, poscan[1].stop)
        po_avg, po_bias, po_std = stats.sigma_clipped_stats(
            hdu.data[poscan], axis=0)
        po_bias = po_bias.reshape(1, np.shape(po_bias)[0])
        logging.debug("np.shape(po_bias)={}".format(np.shape(po_bias)))
        hdu.data = hdu.data - po_bias.data
    else:
        logging.error('btype: %s not valid', btype)
        exit(1)
