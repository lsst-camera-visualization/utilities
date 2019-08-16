#!/usr/bin/env python
"""
Display line plots from fits images
"""

import re
import argparse
import logging
import textwrap
import math
from astropy.io import fits
from astropy import stats
import numpy as np
import matplotlib.pyplot as plt
import imutils as iu


def parse_args():
    """handle command line"""
    # style_list = ['default', 'classic'] + sorted(
    #      style for style in plt.style.available if style != 'classic')
    style_list = ['default'] + sorted(
        ['fast', 'ggplot', 'seaborn-poster', 'seaborn-notebook'])
    #
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        Line plots of row/column regions
                                    '''),
        epilog=textwrap.dedent('''\
        For each file a plot is produced yielding a grid of plots
        unless the "--overlay" option is used which makes just one plot.
        Within each file the specified hdu's (or all) are plotted for
        each region given.
        Each region is collapsed to 1-d by taking the mean, median or
        clipped median on the other axis from the plot.  That is for a
        row plot, each column of the region is mapped to a scalar.
        The "--bias" and "btype" options provide bias subtraction.
        N.B. a "--" is often needed before the file(s) to indicate
        the end of options.
                                '''))
    parser.add_argument("fitsfile", nargs='+',
                        metavar="file", help="input fits file(s)")
    # row or column plot stuff affecting drawing the plots
    pgroup = parser.add_mutually_exclusive_group()
    pgroup.add_argument("--row", nargs='*', metavar='reg',
                        help="row plot, region fmt: \"x1:x2,y1:y2\"")
    pgroup.add_argument("--col", nargs='*', metavar='reg',
                        help="column plot, region fmt: \"x1:x2,y1:y2\"")
    parser.add_argument("--ltype", choices=['median', 'mean', 'clipped'],
                        default='median', help="line type: 2d->1d method")
    parser.add_argument("--offset", choices=['mean', 'median', 'delta'],
                        help="plot relative from mean, median or w/offset")
    parser.add_argument("--overlay", action='store_true',
                        help="all lines in one plot")
    parser.add_argument("--bias", nargs='?', metavar='cols', const='overscan',
                        help="subtract bias, fmt: \"x1:x2\"")
    parser.add_argument("--btype", default='byrow',
                        choices=['mean', 'median', 'byrow', 'byrowcol'],
                        help="bias subtraction method, default is byrow")
    ygroup = parser.add_mutually_exclusive_group()
    ygroup.add_argument('--sharey', dest='sharey', action='store_true',
                        help="Share Y-axis range, this is default")
    ygroup.add_argument('--no-sharey', dest='sharey', action='store_false',
                        help="Don't share Y-axis range")
    ygroup.set_defaults(sharey=True)
    # hdu name|index exclusive
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument("--hduname", nargs='+',
                        metavar='idn', help="process HDU list by names")
    hgroup.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="process HDU list by ids")
    parser.add_argument("--info", action='store_true',
                        help="print the info() table summarizing file")
    # title
    parser.add_argument("--title", nargs='?', metavar='Title',
                        const='auto', help="specify Title String")
    # fluff
    parser.add_argument("--layout", default='landscape',
                        help="\"landscape\"|\"portrait\"|\"nxm\"")
    parser.add_argument("--style", default='ggplot', required=False,
                        choices=style_list)
    parser.add_argument("--steps", action='store_const', required=False,
                        const='steps-mid', default='default',
                        help='Use \'steps-mid\' style, best for short lines')
    parser.add_argument("--debug", action='store_true',
                        help="print additional debugging messages")
    parser.add_argument("--nomemmap", action='store_true', default=False,
                        help="for case where memmap doesn't work")
    return parser.parse_args()


def main():
    """main logic:"""
    optlist = parse_args()
    iu.init_logging(optlist.debug)
    ncalls.counter = 0

    # layout of plot view for array of plots:
    #     landscape has columns > rows
    #     portrait has more rows > columns w/room for legend
    plt.style.use(optlist.style)
    plt.rcParams['figure.constrained_layout.use'] = True

    fig, axes = get_fig_and_axis(len(optlist.fitsfile), optlist.layout,
                                 optlist.overlay, optlist.sharey)

    fsize = fig.get_size_inches()
    logging.debug("width= %5.2f, height= %5.2f", fsize[0], fsize[1])
    logging.debug("len(axes)={}".format(len(axes)))
    logging.debug('axes.shape= %s', axes.shape)
    nprows, npcols = (axes.shape[0], axes.shape[1])
    logging.debug('nprows= %d, npcols= %d', nprows, npcols)

    set_title(optlist.title, optlist.fitsfile, fig)
    nfiles = len(optlist.fitsfile)

    # findex = 0
    # for ffile in optlist.fitsfile:
    for findex in range(0, nfiles):
        try:
            hdulist = fits.open(optlist.fitsfile[findex],
                                memmap=optlist.nomemmap)
        except IOError as ioerr:
            logging.error('IOError: %s', ioerr)
            exit(1)
        # info option
        if optlist.info:
            hdulist.info()
            continue

        if optlist.overlay:
            ax = axes[0, 0]
        else:
            logging.debug('ax = np.ravel(axes)[%d + %d]',
                          int(findex / npcols) * npcols, findex % npcols)
            ax = np.ravel(axes)[
                int(findex / npcols) * npcols + findex % npcols]

        # construct a list of the HDU's to work on
        hduids = iu.get_requested_image_hduids(
            hdulist, optlist.hduname, optlist.hduindex)
        if hduids is None:
            logging.error('No valid HDUs found in %s',
                          optlist.hduname or optlist.hduindex)
            exit(1)

        # title string is truncated filename (w/out path or .fit(s))
        title_str = re.sub(r"^.*/(.*)$", r"\1", optlist.fitsfile[findex])
        title_str = re.sub(r"^(.*)\.fits?(\.fz)*$", r"\1", title_str)
        title_str = "{}...".format(title_str[:27])
        if optlist.overlay:
            if nfiles > 1:  # empty line=> title_str in legend
                ax.plot([], [], ' ', label=title_str)
        elif nfiles == 1 and not optlist.title:
            ax.set_title("{}".format(title_str))
        else:
            ax.set_title("{}".format(title_str))

        # y label depends on offset type
        if not optlist.offset:
            ax.set_ylabel('signal')
        elif optlist.offset == 'mean':
            ax.set_ylabel('signal - mean')
        elif optlist.offset == 'median':
            ax.set_ylabel('signal - median')
        elif optlist.offset == 'delta':
            ax.set_ylabel('signal - mean + 5*j*stdev, j=0,1,..')
        else:
            logging.error('invalid --offset choice')
            exit(1)
        # x label
        ax.grid(True)
        if optlist.row is not None:
            ax.set_xlabel('column')
        elif optlist.col is not None:
            ax.set_xlabel('row')
        else:
            logging.error('must have one of --row or --col')
            exit(1)

        # do the plotting
        plot_hdus(optlist, hduids, hdulist, ax)

        #  done with file, close it
        hdulist.close()

        # reset counter for each file
        if not optlist.overlay:
            ncalls.counter = 0

        # end of loop over files

    if optlist.info:  # just print the image info and exit
        exit()

    # Deal with the legend (ugly)
    if optlist.overlay:
        ax = np.ravel(axes)[0]
    else:
        ax = np.ravel(axes)[
            int((nfiles - 1) / npcols) * npcols + (nfiles - 1) % npcols]
    handles, labels = ax.get_legend_handles_labels()
    if nfiles == 1 or optlist.overlay:
        ax = axes[0, 0]
        if len(handles) < 6:  # inside the plot box
            ax.legend(handles, labels, loc='best',
                      fontsize='small', title='HDU:Region')
        else:  # outside at upper right
            fsize = fig.get_size_inches()
            fig.set_size_inches(fsize[0]*1.2, fsize[1])
            ax.legend(handles, labels, loc="upper left",
                      bbox_to_anchor=(1, 1), fontsize='small',
                      title='HDU:Region')
    else:
        ax = np.ravel(axes)[-1]   # use the last slot
        if len(handles) < 6:
            fsize = 'small'
        elif len(handles) < 12:
            fsize = 'x-small'
        elif len(handles) < 25:
            fsize = 'xx-small'
        else:
            labels[24] = "truncated labels[{}:{}]".format(
                24, len(handles))
        ax.legend(handles[:25], labels[:25],
                  fontsize=fsize, title='HDU:Region ')

    if not optlist.overlay:
        for gidx in range(nfiles, nprows * npcols):
            ax = np.ravel(axes)[
                int(gidx / npcols) * npcols + gidx % npcols]
            ax.grid(False)
            ax.set_frame_on(False)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

    plt.show()


def plot_hdus(optlist, hduids, hdulist, ax):
    """
    make a plot according to optlist for hduids from hdulist on ax
    """
    if optlist.row:
        map_axis = 0  # first axis (y) converts to scalar
    elif optlist.col:
        map_axis = 1  # second axis (x) converts to scalar
    elif optlist.gradient:
        logging.debug('not implemented yet')
    else:
        exit(1)
    # Process each HDU in the list "hduids"
    for hduid in hduids:
        hdu = hdulist[hduid]
        try:
            name = hdu.name
        except IndexError as ierr:
            logging.debug('IndexError: %s', ierr)
            logging.debug('using name=%s', hduid)
            name = "{}".format(hduid)
        if optlist.bias:
            iu.subtract_bias(optlist.bias, optlist.btype, hdu)
        (datasec, soscan, poscan) = iu.get_data_oscan_slices(hdu)
        slices = []  # define regions to plot
        for reg in optlist.row or optlist.col:
            logging.debug('processing %s', reg)
            if re.match(r"data", reg):
                slice_spec = datasec
            elif re.match(r"over", reg):
                slice_spec = soscan
            elif re.match(r"pover", reg):
                slice_spec = poscan
            else:
                slice_spec = iu.parse_region(reg)
            if slice_spec is not (None, None):
                slices.append(slice_spec)
            else:
                logging.error('skipping region %s', reg)
        for slice_spec in slices:
            logging.debug('calling line_plot() %s[%s]', name, reg)
            line_plot(slice_spec, hdu.data, optlist.ltype, optlist.steps,
                      optlist.offset, map_axis, name, ax)
            ncalls()


def line_plot(slice_spec, pix, plot_type, steps,
              plot_offset, map_axis, hduname, ax):
    """
    make a line plot on the axis given (ax) using the region etc.
    """
    logging.debug('slice_spec: %s', slice_spec)
    logging.debug('shape(pix[slice_spec])=%s', np.shape(pix[slice_spec]))

    #  transform region to a row or column (average, median, ...)
    logging.debug('plot_type=%s', plot_type)
    if not plot_type:
        line1 = np.median(pix[slice_spec], axis=map_axis)
    elif plot_type == 'median':
        line1 = np.median(pix[slice_spec], axis=map_axis)
    elif plot_type == 'mean':
        line1 = np.mean(pix[slice_spec], axis=map_axis)
    elif plot_type == 'clipped':
        l1_avg, line1, l1_std = stats.sigma_clipped_stats(
            pix[slice_spec], axis=map_axis)
        logging.debug(
            'shape(l1_avg)=%s, shape(pix[slice_spec])=%s, shape(l1_std)=%s',
            np.shape(l1_avg), np.shape(pix[slice_spec]), np.shape(l1_std))
    else:
        logging.error('ltype incorrect or programmer error')
        exit(1)

    #  shift the line (row|column) if needed for plotting
    if plot_offset:
        logging.debug('plot_offset=%s', plot_offset)
        avg, med, std = stats.sigma_clipped_stats(line1)
        if plot_offset == 'mean':
            line1 = line1 - avg
        elif plot_offset == 'median':
            line1 = line1 - med
        elif plot_offset == 'delta':
            # offset by 5 std's per line
            line1 = line1 - med + ncalls.counter * 5 * std
        else:
            logging.error('invalid offset or programmer error')
            exit(1)

    #  plot the line: N.B. arrays are 0 indexed, fits is 1 indexed
    x = np.arange((slice_spec[1].start or 0) + 1,
                  (slice_spec[1].stop or len(pix[0, :])) + 1)
    y = np.arange((slice_spec[0].start or 0) + 1,
                  (slice_spec[0].stop or len(pix[:, 0])) + 1)
    s = x if map_axis == 0 else y
    slabel = "{}:[{}:{},{}:{}]".format(hduname, x[0], x[-1], y[0], y[-1])
    ax.plot(s, line1, drawstyle="{}".format(steps), label=slabel)


def set_title(title, fitsfile, fig):
    """
    Given the title option string, the file list and the figure,
    generate the suptitle.
    """
    if title:
        if title == "auto":
            filenames = []
            for fname in fitsfile:
                fname = re.sub(r"^.*/(.*)$", r"\1", fname)
                fname = re.sub(r"^(.*)\.fits?(\.fz)*$", r"\1", fname)
                filenames.append(fname)
            logging.debug("using autogenerated title...")
            logging.debug("filenames[]={}".format(filenames))
            if len(filenames) > 1:
                # title is longest common substring array
                # joined with *'s to look like a glob pattern
                ss_arr = []
                iu.get_lcs_array(filenames, ss_arr, 0, '', 2)
                if ss_arr:
                    fig.suptitle("*{}*".format('*'.join(ss_arr)))
            else:
                fig.suptitle("{}".format(filenames[0]))
        else:
            logging.debug("using title=%s", title)
            fig.suptitle("{}".format(title))


def get_fig_and_axis(nfiles, layout, overlay, sharey):
    """
    Create the figure and subplot axes based on number of input files
    and command line options.
    Note plt.subplot() is called with squeeze=False which forces axes
    to be a proper array of Axes objects even if there is only one.
    """
    fsize = plt.rcParams["figure.figsize"]
    if nfiles == 1 or overlay:
        npcols = nprows = 1
    else:
        # default is landscape
        npcols = int(math.ceil(math.sqrt(nfiles)))
        nprows = int(math.ceil(float(nfiles)/npcols))
        while float(npcols)/float(nprows) < 1.3:
            npcols += 1
            nprows = int(math.ceil(float(nfiles)/npcols))

        # add extra for legend if needed
        if npcols*nprows < (nfiles + 1):
            npcols += 1

        # portrait
        if layout == 'portrait':
            nprows, npcols = npcols, nprows
            fsize[0], fsize[1] = fsize[1], fsize[0]

        if re.match(r"([1-9]+)x([1-9]+)$", layout):
            (rstr, cstr) = re.match(
                r"([1-9]+)x([1-9]+)$", layout).groups()
            nprows, npcols = int(rstr), int(cstr)
            if nprows > npcols:
                fsize[0], fsize[1] = fsize[1], fsize[0]

    logging.debug('subplots layout for nfiles=%s is nprows=%s x npcols=%s',
                  nfiles, nprows, npcols)

    fig, axes = plt.subplots(nprows, npcols, sharey=sharey, squeeze=False,
                             figsize=(fsize[0], fsize[1]))
    return (fig, axes)


def ncalls():
    """maintain a counter
    """
    ncalls.counter += 1


if __name__ == '__main__':
    main()
