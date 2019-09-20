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
    """ Set up debug and info level logging
    """
    if debug:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s',
                            level=logging.INFO)


def init_warnings():
    """ Block warnings from Astropy
    """
    warnings.simplefilter('ignore', category=AstropyWarning)


def create_output_hdulist(hdulisti: fits.HDUList, argv: list) -> fits.HDUList:
    """
    Create output HDUList from input HDUList for building new image

    DATE and an HISTORY header cards added to record what was done
    This is generally the first step before subsequent ops to modify
    data arrays and changing additional header keys.
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


def init_hdu(hdui: fits.ImageHDU, hdulisto: fits.HDUList,
             region: tuple = None) -> fits.ImageHDU:
    """
    Append a new image HDU to output image using input HDU as a template.

    Copy the header and set the size/region specs in preparation for data
    to be added later.

    Returns
    -------
    hduo: fits.ImageHDU That was created during the call.

    """
    # create the output hdu from the master, (primary already exists)
    if not isinstance(hdui, fits.PrimaryHDU):
        hdri = hdui.header.copy()
        hdulisto.append(fits.ImageHDU(None, hdri, hdri['EXTNAME']))
    hduo = hdulisto[len(hdulisto) - 1]
    hdro = hduo.header
    hdro['NAXIS'] = 2
    hdro.set('NAXIS1', hdri['NAXIS1'], "size of the n'th axis", after='NAXIS')
    hdro.set('NAXIS2', hdri['NAXIS2'], "size of the n'th axis", after='NAXIS1')
    hdro['BITPIX'] = -32
    # make changes to account for region of interest subimage
    if region and region is not (None, None):
        logging.debug('region = {}'.format(region))
        naxis2 = ((region[0].stop or len(hdui.data[:, 0]))
                  - (region[0].start or 0))
        naxis1 = ((region[1].stop or len(hdui.data[0, :]))
                  - (region[1].start or 0))
        hdro.set('NAXIS1', naxis1,
                 "size of the n'th axis", after='NAXIS')
        hdro.set('NAXIS2', naxis2,
                 "size of the n'th axis", after='NAXIS1')
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
    # logging.debug('output header:\n%s\n', hdro.tostring())
    return hduo


def parse_region(reg: str) -> tuple:
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
            r"([0-9]+):([0-9]+),\s*([0-9]+):([0-9]+)$", reg).groups()
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
    # reg = [x1:*,y1:y2]) -- row selection w/cols to end
    elif re.match(r"([0-9]+):\s*(\*),\s*([0-9]+):([0-9]+)$", reg):
        (x1, x2, y1, y2) = re.match(
            r"([0-9]+):\s*(\*),\s*([0-9]+):([0-9]+)$", reg).groups()
        retval = (slice(int(y1)-1, int(y2)), slice(int(x1)-1, None))
    #
    # reg = [*:x1,y1:y2]) -- row selection w/cols from beginning
    elif re.match(r"(\*):\s*([0-9]+),\s*([0-9]+):([0-9]+)$", reg):
        (x1, x2, y1, y2) = re.match(
            r"(\*):\s*([0-9]+),\s*([0-9]+):([0-9]+)$", reg).groups()
        retval = (slice(int(y1)-1, int(y2)), slice(None, int(x2)-1))
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


def get_requested_image_hduids(hdulist: fits.HDUList,
                               hdunames:list, hduindices:list) -> list:
    """
    Return a list of imaage hduids requested in optlist or all by default.

    Check that they exist in hdulist and have data.  Requested hduids that
    don't exist are skipped.  Redundant values are dropped.
    """
    chduids = []  # list of candidate hduids
    logging.debug('hdunames= %s', hdunames)
    for name in hdunames or []:
        try:
            hduid = hdulist.index_of(name)
            if hduid not in chduids:
                chduids.append(hduid)
        except KeyError as ke:
            logging.error('KeyError: %s', ke)
            logging.error('HDU[%s] not found, skipping', name)
    for hduid in hduindices or []:
        try:
            hdu = hdulist[hduid]
            if hduid not in chduids:
                chduids.append(hduid)
        except IndexError:
            logging.error('HDU[%d] not found, skipping', hduid)
    if not hduindices and not hdunames:
        for hdu in hdulist:
            chduids.append(hdulist.index(hdu))

    # Validate the list of candidate HDUs, keep those with pixels
    hduids = []
    for hduid in chduids:
        hdu = hdulist[hduid]
        if isinstance(hdu, fits.PrimaryHDU):  # check for data
            if np.shape(hdu.data):
                logging.debug('adding %s with index %d to hduid list',
                              hdu.name, hduid)
                hduids.append(hduid)
        elif isinstance(hdu, (fits.ImageHDU, fits.CompImageHDU)):
            logging.debug('adding %s with index %d to hduid list',
                          hdu.name, hduid)
            hduids.append(hduid)
        else:
            logging.debug('%s with index %d is not type (Comp)ImageHDU',
                          hdu.name, hduid)
    if hduids:
        return hduids
    return None


def get_data_oscan_slices(hdu: fits.FitsHDU) -> tuple:
    """
    Get datasec, serial/parallel overscan as slice specifications.

    Given an hdu, uses header keys to infer slice specs.  If a particular
    region cannot be obtained a spec of (None, None) is returned for that region.
    Returns a tuple of slice definitions (datasec, soscan, poscan).
    The serial overscan is assumed to be at the end of each row is present.
    """
    # first get serial and parallel overscan region defs
    hdr = hdu.header
    try:
        dstr = hdr['DATASEC']
    except KeyError as ke:
        logging.error('KeyError: %s required', ke)
        return (None, None, None)
    logging.debug('DATASEC=%s', dstr)
    try:
        n1 = hdr['NAXIS1']
    except KeyError as ke:
        logging.error('KeyError: %s required', ke)
        return (None, None, None)
    try:
        n2 = hdr['NAXIS2']
    except KeyError as ke:
        logging.error('KeyError: %s required', ke)
        return (None, None, None)
    # get DATASEC region
    datasec = parse_region(dstr)
    if datasec == (None, None):
        return (None, None, None)
    (p1, p2) = (datasec[0].start or 0,
                datasec[0].stop or len(hdu.data[:, 0]))
    (s1, s2) = (datasec[1].start or 0,
                datasec[1].stop or len(hdu.data[0, :]))
    if n1 > s2:
        soscan = (slice(0, n2), slice(s2, n1))
    else:  # no serial overscan
        soscan = (slice(None), slice(None))
    if n2 > p2:
        poscan = (slice(p2, n2), slice(0, n1))
    else:
        poscan = (slice(None), slice(None))

    return (datasec, soscan, poscan)


def subtract_bias(bstring: str, btype: str, hdu: fits.ImageHDU) -> fits.ImageHDU:
    """
    Subtract a bias estimate (from overscans) from an hdu.

    Choices are mean, median, byrow or byrowcol subtraction of a bias
    calculated in either a given set of columns or using DATASEC to infer
    the overscan region. The rows-and-columns method uses both the serial
    and parallel oscan to account for the shape of the bias
    (particularly for ITL sensors)
    """
    (datasec, soscan, poscan) = get_data_oscan_slices(hdu)

    if bstring == 'overscan':
        logging.debug('use overscan region in hdu:%s', hdu.name)
    elif len(bstring) >= 3:  # user supplies column selection
        # eg. a:b is minimum spec size, peel off outer brackets
        try:
            logging.debug('use %s cols in hdu:%s bias subtraction',
                          bstring, hdu.name)
            reg = re.sub(r"^\[*([0-9]+:[0-9]+)\]*$", r"\1", bstring)
            (x1, x2) = re.match(r"([0-9]+):([0-9]+)$", reg).groups()
            soscan = (soscan[0], slice(int(x1), int(x2)))
        except SyntaxError as se:
            logging.error('SyntaxError: %s', se)
            logging.error('bad bias selection %s', bstring)
            exit(1)
    else:
        logging.error('bad bias selection %s', bstring)
        exit(1)

    logging.debug('biastype=%s', btype)
    # default method is "byrow"
    if btype in ('byrow', 'byrowcol'):
        # get sigma clipped median per row
        so_avg, so_bias, so_std = stats.sigma_clipped_stats(
            hdu.data[soscan], axis=1)
        so_bias = so_bias.reshape(np.shape(so_bias)[0], 1)
        logging.debug("np.shape(so_bias)={}".format(np.shape(so_bias)))
        hdu.data = hdu.data - so_bias.data
        if btype == 'byrowcol':
            # get sigma clipped median per column
            logging.debug('poscan=((%d, %d), (%d, %d))',
                          poscan[0].start, poscan[0].stop,
                          poscan[1].start, poscan[1].stop)
            po_avg, po_bias, po_std = stats.sigma_clipped_stats(
                hdu.data[poscan], axis=0)
            po_bias = po_bias.reshape(1, np.shape(po_bias)[0])
            logging.debug("np.shape(po_bias)={}".format(np.shape(po_bias)))
            hdu.data = hdu.data - po_bias.data
    elif btype == 'mean':
        bias = np.mean(hdu.data[soscan])
        hdu.data = hdu.data - bias
    elif btype == 'median':
        bias = np.median(hdu.data[soscan])
        hdu.data = hdu.data - bias
    else:
        logging.error('btype: %s not valid', btype)
        exit(1)


def long_substr(ilist: list) -> str:
    """
    Find the longest common substring from a list of strings.

    This is from:
    https://stackoverflow.com/questions/2892931/\
        longest-common-substring-from-more-than-two-strings-python#
    """
    substr = ''
    if len(ilist) > 1 and len(ilist[0]) > 0:
        for i in range(len(ilist[0])):
            for j in range(len(ilist[0])-i+1):
                if j > len(substr) and all(ilist[0][i:i+j] in x for x in ilist):
                    substr = ilist[0][i:i+j]
    return substr


def cleave(arr: list, substr: str) -> (list, list):
    """
    Split each element in array of strings using a substring.

    Returns
    -------
    pre_arr, post_arr :
        Lists of strings with preceeding/following segments after splitting.
        If either preceeding/following segment is empty, (None, None) is returned.
    """

    cleav_pat = re.compile("^(.*){}(.*)$".format(substr))
    pre_arr = []
    post_arr = []
    for istr in arr:
        tl, tr = re.match(cleav_pat, istr).groups()
        if tl and tr:
            pre_arr.append(tl)
            post_arr.append(tr)
        else:  # if any pre/post string is 0 length this branch is done
            return (None, None)

    return (pre_arr, post_arr)


def get_lcs_array(s_arr: list,  ss_arr: list, index: int,
                  order: str, minsz: int) -> list:
    """
    Create ordered list of common substrings for list of strings.

    Given a list of filenames or paths or similar, find the set of all
    common substrings longer than 'minsz'.  The result can be used to create
    a 'glob-like' pattern to summarize the original list of strings.

    Parameters
    ----------
    s_arr: list Array of strings to process to find lcs.

    ss_arr: list Substring array -- this is the product.

    index: int Index in ss_arr of parent lcs.

    order: str Insert new lcs before|after index.

    minsz: int Shortest allowed common substring.
    """

    lcstr = long_substr(s_arr)  # get longest common substring
    if len(lcstr) > minsz:
        if len(ss_arr) == 0:
            logging.debug('lcstr=%s', lcstr)
            logging.debug("first: ss_arr={}".format(ss_arr))
            logging.debug("first: index={}".format(index))
            ss_arr.append(lcstr)
            logging.debug("first: ss_arr={}".format(ss_arr))
            new_index = 0
        else:
            logging.debug('lcstr=%s', lcstr)
            if order == "before":
                logging.debug("before: ss_arr={}".format(ss_arr))
                logging.debug("before: index={}".format(index))
                ss_arr.insert(index, lcstr)
                logging.debug("before: ss_arr={}".format(ss_arr))
                new_index = index
            elif order == "after":
                # treat index as offset from end
                logging.debug("after: ss_arr={}".format(ss_arr))
                logging.debug("after: index={}".format(index))
                if index == 0:
                    ss_arr.append(lcstr)
                else:
                    ss_arr.insert(-1 * index, lcstr)
                logging.debug("after: ss_arr={}".format(ss_arr))
                new_index = index  # stay at same relative place from end
            else:
                logging.error('invalid order %s', order)
                exit(1)
    else:
        return

    pre_arr, post_arr = cleave(s_arr, lcstr)
    if pre_arr:
        if order == "after":
            new_index = len(ss_arr) - new_index - 1
        get_lcs_array(pre_arr, ss_arr, new_index, 'before', minsz)
    if post_arr:
        if order == "before":
            new_index = len(ss_arr) - new_index - 1
        get_lcs_array(post_arr, ss_arr, new_index, 'after', minsz)
