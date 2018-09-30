#!/usr/bin/env python
"""
Display plots from fits images
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
        region given.
        Each region is collapsed to 1-d by taking the mean, median or
        clipped median on the other axis from the plot.  That is for a
        row plot, columns in the region are mapped to a scalar.
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
                        help="overlay all plots in one (default: per file)")
    # hdu name|index exclusive
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument("--hduname", nargs='+',
                        metavar='idn', help="process HDU list by names")
    hgroup.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="process HDU list by ids")
    parser.add_argument("--info", action='store_true',
                        help="print the info() table summarizing file")
    # title
    parser.add_argument("--title", nargs='?', help="specify TitleString")
    # fluff
    parser.add_argument("--layout", default='landscape', required=False,
                        choices=['landscape', 'portrait'])
    parser.add_argument("--style", default='default', required=False,
                        choices=style_list)
    parser.add_argument("--debug", action='store_true',
                        help="print additional debugging messages")
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

    # layout of plot view for array of plots:
    #     landscape has columns > rows, target ratio is sqrt(2)
    #     portrait has more rows > columns, leave room for legend
    # plt.style.use('fivethirtyeight')
    plt.style.use(optlist.style)
    nfiles = len(optlist.fitsfile)
    fig_width = plt.rcParams["figure.figsize"][0]
    fig_height = plt.rcParams["figure.figsize"][1]
    if nfiles == 1:
        npcols = nprows = 1
    else:
        # default is landscape
        npcols = int(math.ceil(math.sqrt(nfiles)))
        nprows = int(math.ceil(float(nfiles)/npcols))
        # while float(npcols)/float(nprows) < 0.9*math.sqrt(2):
        while float(npcols)/float(nprows) < 1.3:
            npcols += 1
            nprows = int(math.ceil(float(nfiles)/npcols))

        # add extra for legend if needed
        if npcols*nprows < (nfiles + 1):
            npcols += 1

        # portrait
        if optlist.layout == 'portrait':
            nprows, npcols = npcols, nprows
            fig_width, fig_height = fig_height, fig_width

    logging.debug('subplots layout for nfiles=%s is nprows=%s x npcols=%s',
                  nfiles, nprows, npcols)

    if optlist.overlay:
        fig, axes = plt.subplots(1, 1)
    else:
        fig, axes = plt.subplots(nprows, npcols, sharey=True,
                                 figsize=(fig_width, fig_height))
        logging.debug("width= %5.2f, height= %5.2f", fig_width, fig_height)
    # title_string = long_substr(optlist.fitsfile)
    # fig.suptitle("*{}*...".format(title_string))
    if optlist.title:
        fig.suptitle("{}".format(optlist.title))

    # begin processing -- loop over rows, cols
    fcnt = 0
    for prow in range(nprows):
        for pcol in range(npcols):
            if fcnt < nfiles:
                #
                ffile = optlist.fitsfile[fcnt]
                logging.debug('processing %s', ffile)
                logging.debug('prow=%s pcol=%s', prow, pcol)

                # choose which subplot
                if np.shape(axes):
                    logging.debug('np.shape(axes)=%s', np.shape(axes))
                    # ax = axes[prow, pcol]
                    ax = np.ravel(axes)[prow*npcols + pcol]
                else:
                    ax = axes

                # open file, no memmap for now...
                try:
                    hdulist = fits.open(ffile, memmap=False)
                except IOError as ioerr:
                    logging.error('IOError: %s', ioerr)
                    exit(1)

                # info option
                if optlist.info:
                    hdulist.info()
                    continue

                # construct a list of the HDU's to work on
                hduids = iu.get_hduids(optlist, hdulist)

                # title string is filename (w/out path or .fit(s))
                if not optlist.overlay:
                    title_str = re.sub(r"^.*/(.*)$", r"\1", ffile)
                    title_str = re.sub(r"^(.*)\.fits?$", r"\1", title_str)

                # labels, grid, title for ax
                ax.set_ylabel('signal')
                ax.grid(True)
                ax.set_title("{}".format(title_str))
                if optlist.row is not None:
                    ax.set_xlabel('column')
                elif optlist.col is not None:
                    ax.set_xlabel('row')
                else:
                    exit(1)

                # do the plotting
                plot_hdus(optlist, hduids, hdulist, ax)

                #  done with file, close it
                hdulist.close()

            else:  # files are done, no more boxes or grids
                if not optlist.overlay:
                    # ax = axes[prow, pcol]
                    ax = np.ravel(axes)[prow*npcols + pcol]
                    ax.grid(False)
                    ax.set_frame_on(False)
                    ax.get_xaxis().set_visible(False)
                    ax.get_yaxis().set_visible(False)
                    ncalls.counter = 0
            fcnt += 1
            if not optlist.overlay:  # reset counter for each file
                ncalls.counter = 0
            if fcnt == nfiles:      # last file, get legend stuff
                handles, labels = ax.get_legend_handles_labels()
                if nfiles == 1:
                    if np.size(hduids) < 4:
                        ax.legend(handles, labels, loc='best',
                                  fontsize='xx-small', title='HDUi:Region')
                    else:
                        fig.set_size_inches(fig_width*1.2, fig_height)
                        ax.legend(handles, labels, loc="upper left",
                                  bbox_to_anchor=(1, 1), fontsize='x-small',
                                  title='HDUi:Region ')

            # place the legend in the last slot (bottom right)
            if fcnt > nfiles and fcnt == nprows * npcols:
                ax.legend(handles, labels,
                          fontsize='xx-small', title='HDUi:Region ')

    if optlist.info:  # just print the image info and exit
        exit()
    if optlist.title:
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    else:
        fig.tight_layout()
    plt.show()


def plot_hdus(optlist, hduids, hdulist, ax):
    """
    make a plot according to optlist for hduids from hdulist on ax
    """

    for hduid in hduids:  # process each with regions
        name = hdulist[hduid].name
        reg = ""
        if optlist.row:
            plot_axis = 0
        elif optlist.col:
            plot_axis = 1
        else:
            exit(1)

        for reg in optlist.row or optlist.col:
            slice_spec = iu.parse_region(reg)
            if slice_spec:
                logging.debug('calling line_plot() %s[%s]', name, reg)
                line_plot(slice_spec, hdulist[hduid].data, optlist.ltype,
                          plot_axis, optlist.offset, hdulist[hduid].name,
                          ax)
            else:
                logging.error('skipping')
            ncalls()


def line_plot(slice_spec, hdudata, plot_type,
              plot_axis, plot_offset, hduname, ax):
    """
    make a line plot on the axis given (ax) using the region etc.
    """
    #  extract the region data
    logging.debug('slice_spec: %s', slice_spec)
    (y1, y2, x1, x2) = (slice_spec[0].start, slice_spec[0].stop,
                        slice_spec[1].start, slice_spec[1].stop)
    logging.debug('slice_spec: x1=%d, x2=%d, y1=%d, y2=%d', x1, x2, y1, y2)
    imbuf = hdudata[y1:y2, x1:x2]

    #  transform region to a row or column (average, median, ...)
    logging.debug('plot_type=%s', plot_type)
    if not plot_type:
        line1 = np.median(imbuf, axis=plot_axis)
    elif plot_type == 'median':
        line1 = np.median(imbuf, axis=plot_axis)
    elif plot_type == 'mean':
        line1 = np.mean(imbuf, axis=plot_axis)
    elif plot_type == 'clipped':
        l1_avg, line1, l1_std = stats.sigma_clipped_stats(
            imbuf, axis=plot_axis)
        logging.debug('shape(l1_avg)=%s, shape(l1_std)=%s',
                      np.shape(l1_avg), np.shape(l1_std))
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
            logging.debug('ncalls.counter=%s', ncalls.counter)
            line1 = line1 - med + ncalls.counter * 5 * std
        else:
            logging.error('invalid offset or programmer error')
            exit(1)

    #  plot the line
    if plot_axis == 0:
        x = np.arange(x1, x2)
    elif plot_axis == 1:
        x = np.arange(y1, y2)
    logging.debug('shape(x)=%s', np.shape(x))
    logging.debug('shape(line1)=%s', np.shape(line1))
    ax.plot(x, line1, label="{}:[{}:{},{}:{}]".format(hduname, x1, x2, y1, y2))


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


def ncalls():
    """maintain a counter
    """
    ncalls.counter += 1


if __name__ == '__main__':
    main()
