#!/usr/bin/env python
"""
Display plots from fits images
"""

import re
import argparse
import logging
import math
from astropy.io import fits
from astropy import stats
import numpy as np
import matplotlib.pyplot as plt
import imutils as iu


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
    parser.add_argument("--noheadings", action='store_true', default=False,
                        help="Don't print column heads for stats")
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument("--hduname", nargs='+',
                        metavar='idn', help="process HDU list by names")
    hgroup.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="process HDU list by ids")
    pgroup = parser.add_mutually_exclusive_group()
    pgroup.add_argument("--line", nargs='*', metavar='reg',
                        help="region fmt: \"x1:x2,y1:y2\"")
    pgroup.add_argument("--column", nargs='*', metavar='reg',
                        help="region fmt: \"x1:x2,y1:y2\"")
    parser.add_argument("--ltype", choices=['median', 'mean', 'clipped'],
                        default='median', help="line plot type")
    parser.add_argument("--ctype", choices=['median', 'mean', 'clipped'],
                        help="line plot type")
    parser.add_argument("--offset", choices=['mean', 'delta'],
                        help="plot diff from mean or diff w/offset")
    parser.add_argument("--title", nargs='?', help="specify TitleString")
    parser.add_argument("--layout", default='landscape', required=False,
                        choices=['landscape', 'portrait'])
    parser.add_argument("--overlay", action='store_true',
                        help="overlay all plots in one (default: per file)")
    return parser.parse_args()

# define a region for the plot to act on
#
# specify what type of plots to make:
#    -line [average|median|clipped]
#    -column [average|median|clipped]
#    -histogram [bin parameters?]
#    -surface [surface params, view angles, scaling]
#    -scatter [hmm?]
#
# style issues
#    -logy, -logx
#    -overlay hdu level, file level, all levels???
#    -offset [(5%),offsetval]


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

    # layout of plot view
    nfiles = len(optlist.fitsfile)
    if nfiles > 1:
        if optlist.layout == 'portrait':
            npcols = int(math.ceil(math.sqrt(nfiles)))
            nprows = int(math.ceil(float(nfiles)/npcols))
            if npcols*nprows == nfiles & nfiles > 1:
                nprows += 1
        elif optlist.layout == 'landscape':
            nprows = int(math.ceil(math.sqrt(nfiles)))
            npcols = int(math.ceil(float(nfiles)/nprows))
            if npcols*nprows == nfiles & nfiles > 1:
                npcols += 1
    else:
        npcols = nprows = 1

    logging.info('subplots layout for nfiles=%s is nprows=%s x npcols=%s',
                 nfiles, nprows, npcols)
    if optlist.overlay:
        fig, axes = plt.subplots(1, 1, sharey=True)
    else:
        fig, axes = plt.subplots(nprows, npcols, sharey=True)
    #title_string = long_substr(optlist.fitsfile)
    #fig.suptitle("*{}*...".format(title_string))

    # begin processing -- loop over files
    # for ffile in optlist.fitsfile:
    fcnt = 0
    for prow in range(nprows):
        for pcol in range(npcols):
            if fcnt < nfiles:
                ffile = optlist.fitsfile[fcnt]
                logging.debug('processing %s', ffile)
                logging.debug('prow=%s pcol=%s', prow, pcol)
                logging.debug('shape(axes)=%s', np.shape(axes))
                if np.shape(axes) == ():
                    ax = axes
                else:
                    ax = axes[prow, pcol]
                try:
                    hdulist = fits.open(ffile)
                except IOError as ioerr:
                    logging.error('IOError: %s', ioerr)
                    exit(1)
                if optlist.info:  # just print the image info and exit
                    hdulist.info()
                    continue
                # Construct a list of the HDU's to work on
                hduids = iu.get_hduids(optlist, hdulist)
                if optlist.line is not None:
                    lineplot(optlist, hduids, hdulist, ax)
                    ax.set_xlabel('column')
                    ax.set_ylabel('signal')
                    ax.grid(True)
                    ax.set_title("{}".format(ffile))
                    # ax.legend(fontsize='xx-small',title='HDUi:Region ')
                elif optlist.column is not None:
                    columnplot(optlist, hduids, hdulist, ax)
                    ax.set_xlabel('row')
                    ax.set_ylabel('signal')
                    ax.grid(True)
                    ax.set_title('xyz')
                    # ax.legend(fontsize='xx-small',title='HDUi:Region ')
                else:
                    exit(0)
            else:  # files are done, no more boxes or grids
                if not optlist.overlay:
                    ax = axes[prow, pcol]
                    ax.grid(False)
                    ax.set_frame_on(False)
                    ax.get_xaxis().set_visible(False)
                    ax.get_yaxis().set_visible(False)
                    ncalls.counter = 0
            fcnt += 1
            if fcnt == nfiles:
                handles, labels = ax.get_legend_handles_labels()
            if fcnt > nfiles and fcnt == nprows * npcols:
                ax.legend(handles, labels,
                          fontsize='xx-small', title='HDUi:Region ')

    if optlist.info:  # just print the image info and exit
        exit()
    # fig.legend(handles, labels, loc='lower right', fontsize='xx-small')
    # ax.legend(fontsize='xx-small',title='HDUi:Region ')
    # plt.plot(x, column1,
    #         label="{}:[{}:{},{}:{}]".format(name, x1, x2, y1, y2))
    # fig.set_size_inches()
    # ax.legend(fontsize='xx-small',title='HDUi:Region ')
    # fig.legend(loc=7, fontsize='xx-small', title='HDUi:Region ')
    # fig.legend(handles=(), labels=(), loc='center')
    # plt.legend( lines, labels, loc = 'lower center',
    #             bbox_to_anchor = (0,-0.1,1,1),
    #                       bbox_transform = plt.gcf().transFigure )
    # ---handles, labels = ax.get_legend_handles_labels()
    # ---fig.legend(handles, labels, loc='lower right', fontsize='xx-small')
    fig.tight_layout()
    # fig.tight_layout(rect=(0,0,1,0.9))
    plt.show()


def long_substr(data):
    """
    https://stackoverflow.com/questions/2892931/\
        longest-common-substring-from-more-than-two-strings-python#
    """
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                    substr = data[0][i:i+j]
    return substr


def columnplot(optlist, hduids, hdulist, ax):
    """ Plot column data as specified in optlist
    """
    # hdr = hdulist[0].header
    # try:
    #     fname = hdr['FILENAME']
    # except KeyError as ke:
    #     logging.error('KeyError: %s, required', ke)
    #     exit(1)

    for hduid in hduids:  # process each with regions
        name = hdulist[hduid].name
        hdr = hdulist[hduid].header
        reg = ""
        if optlist.column == 0:  # - region will be DATASEC
            logging.debug('shape(hdulist[%s].data)=%s',
                          hduid, np.shape(hdulist[hduid].data))
            try:
                dstr = hdr['DATASEC']
            except KeyError as ke:
                logging.error('KeyError: %s, required', ke)
                exit(1)
            logging.debug('DATASEC=%s', dstr)
            optlist.column.append(dstr)
        for reg in optlist.column:
            if reg.startswith("["):
                reg = reg[1:]
            if reg.endswith("]"):
                reg = reg[:-1]
            (nrows, ncols) = np.shape(hdulist[hduid].data)
            slice_spec = iu.parse_region(reg)
            if slice_spec:
                logging.debug('calling line_plot() %s[%s]', name, reg)
                line_plot(optlist, hduid, name,
                          hdulist[hduid].data, slice_spec, ax)
            else:
                logging.error('skipping')
            logging.debug('calling column_plot() %s[%s]', name, reg)
            column_plot(optlist, hduid, name,
                        hdulist[hduid].data, slice_spec)
    plt.xlabel('column')
    plt.ylabel('adu')
    plt.grid(True)
    # plt.title("{}".format(fname))
    plt.legend(fontsize='x-small', title='HDUi:Region ')
    plt.show()


def column_plot(optlist, hduid, name, hdudata, slice_spec):
    """plot column
    """
    ((y1, y2), (x1, x2)) = slice_spec
    logging.debug('slice_spec: x1=%s, x2=%s, y1=%s, y2=%s', x1, x2, y1, y2)
    imbuf = hdudata[y1:y2, x1:x2]
    column1 = np.median(imbuf, axis=1)
    x = np.arange(y1, y2)
    logging.debug('shape(x)=%s', np.shape(x))
    plt.plot(x, column1,
             label="{}:[{}:{},{}:{}]".format(name, x1, x2, y1, y2))


def lineplot(optlist, hduids, hdulist, ax):
    """
    make a line plot according to optlist for hduids from hdulist on ax
    """
    # hdr = hdulist[0].header
    # try:
    #     fname = hdr['FILENAME']
    # except KeyError as ke:
    #     logging.error('KeyError: %s, required', ke)
    #     exit(1)

    for hduid in hduids:  # process each with regions
        name = hdulist[hduid].name
        hdr = hdulist[hduid].header
        reg = ""
        if not optlist.line:  # region will be DATASEC
            logging.debug('shape(hdulist[%s].data)=%s',
                          hduid, np.shape(hdulist[hduid].data))
            try:
                dstr = hdr['DATASEC']
            except KeyError as ke:
                logging.error('KeyError: %s, required', ke)
                exit(1)
            logging.debug('DATASEC=%s', dstr)
            optlist.line.append(dstr)
        for reg in optlist.line:
            slice_spec = iu.parse_region(reg)
            if slice_spec:
                logging.debug('calling line_plot() %s[%s]', name, reg)
                line_plot(optlist, hduid, name,
                          hdulist[hduid].data, slice_spec, ax)
            else:
                logging.error('skipping')
            ncalls()


def line_plot(optlist, hduid, name, hdudata, slice_spec, ax):
    """
    """
    logging.debug('slice_spec: %s', slice_spec)
    (y1, y2, x1, x2) = (slice_spec[0].start, slice_spec[0].stop,
                        slice_spec[1].start, slice_spec[1].stop)
    logging.debug('slice_spec: x1=%d, x2=%d, y1=%d, y2=%d', x1, x2, y1, y2)
    imbuf = hdudata[y1:y2, x1:x2]
    logging.debug('optlist.ltype=%s', optlist.ltype)
    if not optlist.ltype:
        line1 = np.median(imbuf, axis=0)
    elif optlist.ltype == 'median':
        line1 = np.median(imbuf, axis=0)
    elif optlist.ltype == 'mean':
        line1 = np.mean(imbuf, axis=0)
    elif optlist.ltype == 'clipped':
        line1 = np.mean(stats.sigma_clip(imbuf, axis=0), axis=0)
    else:
        logging.error('ltype incorrect or programmer error')
        exit(1)
    if optlist.offset is not None:
        logging.debug('optlist.offset=%s', optlist.offset)
        avg, med, std = stats.sigma_clipped_stats(line1)
        if optlist.offset == 'mean':
            line1 = line1 - avg
        elif optlist.offset == 'delta':
            logging.debug('ncalls.counter=%s', ncalls.counter)
            line1 = line1 - avg + ncalls.counter * 5 * std
        else:
            logging.error('--offset arg incorrect or programmer error')
    x = np.arange(x1, x2)
    logging.debug('shape(x)=%s', np.shape(x))
    logging.debug('shape(line1)=%s', np.shape(line1))
    ax.plot(x, line1, label="{}:[{}:{},{}:{}]".format(name, x1, x2, y1, y2))


def ncalls():
    """maintain a counter
    """
    ncalls.counter += 1


if __name__ == '__main__':
    main()
