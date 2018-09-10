#!/usr/bin/env python
"""
Display plots from fits images
"""

import re
import argparse
import logging
import os.path
from astropy.io import fits
from astropy import stats
import math
import numpy as np
import matplotlib.pyplot as plt

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
    parser.add_argument("--noheadings", action='store_true',
                        default=False, help="Don't print column heads for stats")
    hgroup = parser.add_mutually_exclusive_group()
    hgroup.add_argument("--hduname", nargs='+',
                        metavar='idn', help="process HDU list by names")
    hgroup.add_argument("--hduindex", nargs='+', type=int,
                        metavar='idx', help="process HDU list by ids")
    pgroup = parser.add_mutually_exclusive_group()
    pgroup.add_argument("--line", nargs='*', metavar='reg',
                        help="region fmt: \"x1:x2,y1:y2\",")
    pgroup.add_argument("--column", nargs='*', metavar='reg',
                        help="region fmt: \"x1:x2,y1:y2\",")
    parser.add_argument("--ltype", choices=['median','mean','clipped'],
                        default='median', help="line plot type")
    parser.add_argument("--ctype", choices=['median','mean','clipped'],
                        help="line plot type")
    parser.add_argument("--offset", choices=['mean', 'delta'],
                        help="plot diff from mean or diff w/offset")
    return parser.parse_args()

#define a region for the plot to act on
#
#specify what type of plots to make:
#    -line [average|median|clipped]
#    -column [average|median|clipped]
#    -histogram [bin parameters?]
#    -surface [surface params, view angles, scaling]
#    -scatter [hmm?]
#
#style issues
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
       npcols = int(math.ceil(math.sqrt(nfiles)))
       nprows = int(math.ceil(float(nfiles)/npcols))
       if npcols*nprows == nfiles & nfiles > 1:
           npcols += 1
    else:
        npcols = nprows = 1
    logging.info("subplots layout for nfiles={} is nprows={} x npcols={}".format(nfiles, nprows, npcols))
    fig, axes = plt.subplots(nprows, npcols, sharey=True)
    titleString = long_substr(optlist.fitsfile)
    fig.suptitle("*{}*...".format(titleString))

    # begin processing -- loop over files
    #for ffile in optlist.fitsfile:
    fcnt = 0
    for row in range(nprows):
        for col in range(npcols):
            if fcnt < nfiles:
                ffile = optlist.fitsfile[fcnt]
                logging.debug("processing {}".format(ffile))
                logging.debug("row={} col={}".format(row, col))
                logging.info("row={} col={}".format(row, col))
                logging.debug("shape(axes)={}".format(np.shape(axes)))
                if np.shape(axes) == ():
                    ax = axes
                else:
                    #ax = axes.flatten()[row * npcols + col]
                    ax = axes[row, col]
                try:
                    hdulist = fits.open(ffile)
                except IOError as ioerr:
                    emsg = "IOError: {}".format(ioerr)
                    logging.error(emsg)
                    exit(1)
                if optlist.info: # just print the image info and exit
                    hdulist.info()
                    continue
                if optlist.line is not None:
                    lineplot(optlist, hdulist, ax)
                    ax.set_xlabel('column')
                    ax.set_ylabel('signal')
                    ax.grid(True)
                    ax.set_title('xyz')
                    #ax.legend(fontsize='xx-small',title='HDUi:Region ')
                elif optlist.column is not None:
                    columnplot(optlist, hdulist)
                    ax.set_xlabel('row')
                    ax.set_ylabel('signal')
                    ax.grid(True)
                    ax.set_title('xyz')
                    #ax.legend(fontsize='xx-small',title='HDUi:Region ')
                else:
                    exit(0)
            else:
                ax = axes.flatten()[fcnt]
                ax.grid(False)
                ax.set_frame_on(False)
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
            ncalls.counter = 0
            fcnt += 1
            if fcnt == nfiles:
                handles, labels = ax.get_legend_handles_labels()
            if fcnt == nprows * npcols:
                #ax.legend(handles, labels, fontsize='x-small',title='HDUi:Region ')
                ax.legend(handles, labels, title='HDUi:Region ')

    if optlist.info: # just print the image info and exit
        exit()
    #fig.legend(handles, labels, loc='lower right', fontsize='xx-small')
    #ax.legend(fontsize='xx-small',title='HDUi:Region ')
    #plt.plot(x, column1,
    #         label="{}:[{}:{},{}:{}]".format(name, x1, x2, y1, y2))
    #fig.set_size_inches()
    #ax.legend(fontsize='xx-small',title='HDUi:Region ')
    #fig.legend(loc=7, fontsize='xx-small', title='HDUi:Region ')
    #fig.legend(handles=(), labels=(), loc='center')
    #plt.legend( lines, labels, loc = 'lower center', bbox_to_anchor = (0,-0.1,1,1),
    #                       bbox_transform = plt.gcf().transFigure )
    #---handles, labels = ax.get_legend_handles_labels()
    #---fig.legend(handles, labels, loc='lower right', fontsize='xx-small')
    fig.tight_layout()
    #fig.tight_layout(rect=(0,0,1,0.9))
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

def columnplot(optlist, hdulist):
    """
    """
    hdr = hdulist[0].header
    try:
        fname = hdr['FILENAME']
    except KeyError as ke:
        emsg = "KeyError: {}, required".format(ke)
        logging.error(emsg)
        exit(1)
    hduids = []  #-- make a list of HDU ids to work on
    if optlist.hduname:
        for hduname in optlist.hduname:
            try:
                hduids.append(hdulist.index_of(hduname))
            except KeyError as ke:
                emsg = "KeyError: {}".format(ke)
                logging.error(emsg)
                exit(1)

    elif optlist.hduindex:
        hduids = optlist.hduindex
    else: #- all segments
        for hdu in hdulist:
            if isinstance(hdu, fits.ImageHDU):
                hduids.append(hdulist.index_of(hdu.name))

    for hduid in hduids: # process each with regions
        name = hdulist[hduid].name
        hdr = hdulist[hduid].header
        reg = ""
        if len(optlist.column) == 0: #- region will be DATASEC
            logging.debug("shape(hdulist[{}].data)={}".format(hduid, np.shape(hdulist[hduid].data)))
            try:
                dstr = hdr['DATASEC']
            except KeyError as ke:
                emsg = "KeyError: {}, required".format(ke)
                logging.error(emsg)
                exit(1)
            debugmsg = "DATASEC={}".format(dstr)
            logging.debug(debugmsg)
            optlist.column.append(dstr)
        for reg in optlist.column:
            if reg.startswith("["):
                reg = reg[1:]
            if reg.endswith("]"):
                reg = reg[:-1]
            (nrows, ncols) = np.shape(hdulist[hduid].data)
            slicespec = slice_region(reg, nrows, ncols)
            logging.debug("calling column_plot() {}[{}]".format(name, reg))
            column_plot(optlist, hduid, name,
                                           hdulist[hduid].data, slicespec)
    plt.xlabel('column')
    plt.ylabel('adu')
    plt.grid(True)
    plt.title("{}".format(fname))
    plt.legend(fontsize='x-small',title='HDUi:Region ')
    plt.show()

def column_plot(optlist, hduid, name, hdudata, slicespec):
    """
    """
    (y1, y2, x1, x2) = slicespec
    logging.debug("slicespec: x1={}, x2={}, y1={}, y2={}".format(x1,x2,y1,y2))
    imbuf = hdudata[y1:y2, x1:x2]
    column1 = np.median(imbuf, axis=1)
    lmin = column1.min() * 0.985
    lmax = column1.max() * 1.015
    x = np.arange(y1, y2)
    logging.debug("shape(x)={}".format(np.shape(x)))
    plt.plot(x, column1,
             label="{}:[{}:{},{}:{}]".format(name, x1, x2, y1, y2))


def lineplot(optlist, hdulist, ax):
    """print statistics for region according to options
    """
    hdr = hdulist[0].header
    try:
        fname = hdr['FILENAME']
    except KeyError as ke:
        emsg = "KeyError: {}, required".format(ke)
        logging.error(emsg)
        exit(1)
    hduids = []  #-- make a list of HDU ids to work on
    if optlist.hduname:
        for hduname in optlist.hduname:
            try:
                hduids.append(hdulist.index_of(hduname))
            except KeyError as ke:
                emsg = "KeyError: {}".format(ke)
                logging.error(emsg)
                exit(1)

    elif optlist.hduindex:
        hduids = optlist.hduindex
    else: #- all segments
        for hdu in hdulist:
            if isinstance(hdu, fits.ImageHDU):
                hduids.append(hdulist.index_of(hdu.name))

    for hduid in hduids: # process each with regions
        name = hdulist[hduid].name
        hdr = hdulist[hduid].header
        reg = ""
        if len(optlist.line) == 0: #- region will be DATASEC
            logging.debug("shape(hdulist[{}].data)={}".format(hduid, np.shape(hdulist[hduid].data)))
            try:
                dstr = hdr['DATASEC']
            except KeyError as ke:
                emsg = "KeyError: {}, required".format(ke)
                logging.error(emsg)
                exit(1)
            debugmsg = "DATASEC={}".format(dstr)
            logging.debug(debugmsg)
            optlist.line.append(dstr)
        for reg in optlist.line:
            if reg.startswith("["):
                reg = reg[1:]
            if reg.endswith("]"):
                reg = reg[:-1]
            (nrows, ncols) = np.shape(hdulist[hduid].data)
            slicespec = slice_region(reg, nrows, ncols)
            logging.debug("calling line_plot() {}[{}]".format(name, reg))
            line_plot(optlist, hduid, name,
                                           hdulist[hduid].data, slicespec, ax)
            ncalls()

def line_plot(optlist, hduid, name, hdudata, slicespec, ax):
    """
    """
    (y1, y2, x1, x2) = slicespec
    logging.debug("slicespec: x1={}, x2={}, y1={}, y2={}".format(x1,x2,y1,y2))
    imbuf = hdudata[y1:y2, x1:x2]
    logging.debug("optlist.ltype={}".format(optlist.ltype))
    if not optlist.ltype:
        line1 = np.median(imbuf, axis=0)
    elif optlist.ltype == 'median':
        line1 = np.median(imbuf, axis=0)
    elif optlist.ltype == 'mean':
        line1 = np.mean(imbuf, axis=0)
    elif optlist.ltype == 'clipped':
        line1 = np.mean(stats.sigma_clip(imbuf, axis=0), axis=0)
    else:
        logging.error("ltype incorrect or programmer error")
        exit(1)
    if optlist.offset is not None:
        logging.debug("optlist.offset={}".format(optlist.offset))
        avg, med, std = stats.sigma_clipped_stats(line1)
        if optlist.offset == 'mean':
            line1 = line1 - avg
        elif optlist.offset == 'delta':
            logging.debug("ncalls.counter={}".format(ncalls.counter))
            line1 = line1 - avg + ncalls.counter * 5 * std
        else:
            logging.error("--offset arg incorrect or programmer error")
    x = np.arange(x1, x2)
    logging.debug("shape(x)={}".format(np.shape(x)))
    logging.debug("shape(line1)={}".format(np.shape(line1)))
    ax.plot(x, line1,
             label="{}:[{}:{},{}:{}]".format(name, x1, x2, y1, y2))

def slice_region(reg, nrows, ncols):
    """ return slice defined by the region spec for input pix array
    """
    #- measure pix array
    #- reg = "x0:x1,y1:y2" -- regular region
    res = re.match(r"([0-9]*):([0-9]+),([0-9]+):([0-9]+)", reg)
    if res:
        (x1, x2, y1, y2) = res.groups()
        return (int(y1)-1, int(y2), int(x1)-1, int(x2))
    #- reg = "x0,y1:y2" -- single column, row selection
    res = re.match(r"([0-9]+),([0-9]+):([0-9]+)", reg)
    if res:
        (x0, y1, y2) = res.groups()
        return (int(y1)-1, int(y2), int(x0)-1, int(x0))
    #- reg = "*,y1:y2" -- row selection, all columns
    res = re.match(r"(\*),([0-9]+):([0-9]+)", reg)
    if res:
        (x, y1, y2) = res.groups()
        return (int(y1)-1, int(y2), 0, ncols)
    # reg = "x0,y0" -- single pixel
    res = re.match(r"([0-9]+),([0-9]+)", reg)
    if res:
        (x0, y0) = res.groups()
        return (int(y0)-1,int(y0), int(x0)-1, int(x0))
    # reg = "x1:x2,y0" -- single row
    res = re.match(r"([0-9]+):([0-9]+),([0-9]+)", reg)
    if res:
        (x1, x2, y0) = res.groups()
        return  (int(y0)-1, int(y0), int(x1)-1, int(x2))
    # reg = "x1:x2,*" -- column selection, all rows
    res = re.match(r"([0-9]+):([0-9]+),(\*)", reg)
    if res:
        (x1, x2, y) = res.groups()
        return (0, nrows, int(x1)-1, int(x2))
    # reg = "*,*" #- redundant, but for completeness
    res = re.match(r"(\*),(\*)", reg)
    if res:
        (x, y) = res.groups()
        return (0, nrows, 0, ncols)
    emsg = "bad region spec {}, no match produced".\
        format(reg)
    logging.error(emsg)
    exit(1)



def quicklook(optlist, hdulist):
    """print quicklook for hdu according to options
    """
    hdr = hdulist[0].header
    try:
        expt = float(hdr['EXPTIME'])
    except KeyError as ke:
        emsg = "KeyError: {}".format(ke)
        logging.warn(emsg)
        emsg = "adu/sec won't be available"
        logging.warn(emsg)
        expt = 0.0
    hduids = []  #-- ids of hdus to operate on
    if optlist.hduname:
        for hduname in optlist.hduname:
            hduids.append(hdulist.index_of(hduname))
    elif optlist.hduindex:
        hduids = optlist.hduindex
    else: #- all segments
        for hdu in hdulist:
            if isinstance(hdu, fits.ImageHDU):
                hduids.append(hdulist.index_of(hdu.name))

    for hduid in hduids:
        #- get header details to extract signal, bias regions
        pix = hdulist[hduid].data
        logging.debug("shape(pix)={}".format(np.shape(pix)))
        name = hdulist[hduid].name
        hdr = hdulist[hduid].header
        try:
            dstr = hdr['DATASEC']
        except KeyError as ke:
            emsg = "KeyError: {}, required for quicklook mode".format(ke)
            logging.error(emsg)
            exit(1)
        debugmsg = "DATASEC={}".format(dstr)
        logging.debug(debugmsg)
        res = re.match(r"\[*([0-9]*):([0-9]+),([0-9]+):([0-9]+)\]*",
                       dstr)
        if res:
            datasec = res.groups()
        else:
            emsg = "DATASEC:{} parsing failed".format(dstr)
            logging.error(emsg)
            exit(1)
        naxis1 = int(hdr['NAXIS1'])
        naxis2 = int(hdr['NAXIS2'])
        #- define region to measure signal level
        x1 = int(datasec[0]) - 1
        x2 = int(datasec[1])
        y1 = int(datasec[2]) - 1
        y2 = int(datasec[3])
        sig_buf = pix[y1:y2, x1:x2]
        debugmsg = "sig_buf = pix[{}:{}, {}:{}]".format(y1, y2, x1, x2)
        logging.debug(debugmsg)
        debugmsg =  "shape(sig_buf)={}".format(np.shape(sig_buf))
        logging.debug(debugmsg)

        #Serial bias_region = "[y1:y2,x1:x2]"
        x1 = int(datasec[1])
        x2 = naxis1
        y1 = int(datasec[2]) - 1
        y2 = int(datasec[3])
        if y1 > y2 or x1 > x2:
            emsg = "No bias region available for datasec={}"\
                " with naxis1={}, naxis2={}".\
                format(datasec, naxis1, naxis2)
            logging.error(emsg)
            exit(1)
        bias_buf = pix[y1:y2, x1:x2]
        debugmsg = "bias_buf = pix[{}:{}, {}:{}]".format(y1, y2, x1, x2)
        logging.debug(debugmsg)
        debugmsg = "shape(bias_buf)={}".format(np.shape(bias_buf))
        logging.debug(debugmsg)

        #parallel bias_region = "[y1:y2,x1:x2]"
        x1 = int(datasec[0]) - 1
        x2 = int(datasec[1])
        y1 = int(datasec[3])
        y2 = naxis2
        if y1 > y2 or x1 > x2:
            emsg = "No bias region available for datasec={}"\
                " with naxis1={}, naxis2={}".\
                format(datasec, naxis1, naxis2)
            logging.error(emsg)
            exit(1)
        p_bias_buf = pix[y1:y2, x1:x2]
        debugmsg = "p_bias_buf = pix[{}:{}, {}:{}]".format(y1, y2, x1, x2)
        logging.debug(debugmsg)
        debugmsg = "shape(p_bias_buf)={}".format(np.shape(p_bias_buf))
        logging.debug(debugmsg)
        quicklook_print(optlist, hduid, name, sig_buf,
                        bias_buf, p_bias_buf, expt)

def quicklook_print(optlist, sid, name, sig_buf,
                    bias_buf, p_bias_buf, expt):
    """perform and print the given statistics quantities
       fields are: mean, bias, signal, noise, adu/s
    """
    #quick_fields = ["mean", "bias", "signal", "noise", "adu/sec"]
    quick_fields = ["mean", "bias", "signal",
                    "noise", "adu/sec", "eper:s-cte", "eper:p-cte"]
    if optlist.tearing:
        quick_fields.append("tearing")
    if optlist.dipoles:
        quick_fields.append("dipoles")
    if not optlist.noheadings and ncalls.counter == 0:
        print("#{:>3s} {:>9s}".format("id", "HDUname"), end="")
        if "mean" in  quick_fields:
            print(" {:>6s}".format("median"), end="")
        if "bias" in  quick_fields:
            print(" {:>5s}".format("bias"), end="")
        if "signal" in  quick_fields:
            print(" {:>6s}".format("signal"), end="")
        if "noise" in  quick_fields:
            print(" {:>7s}".format("noise"), end="")
        if "adu/sec" in  quick_fields and expt > 0:
            print("{:>8s}".format("adu/sec"), end="")
        if "eper:s-cte" in  quick_fields:
            print("{:>9s}".format("s-cte"), end="")
        if "eper:p-cte" in  quick_fields:
            print("{:>9s}".format("p-cte"), end="")
        if "tearing" in  quick_fields:
            print("{:>15s}".format("tearing: L  R"), end="")
        if "dipoles" in  quick_fields:
            print("{:>8s}".format("%dipoles"), end="")
        print("") #-- newline

    if not optlist.noheadings:
        print(" {:3d} {:>9s}".format(sid, name), end="")

    if "mean" in  quick_fields:
        sig_mean = np.median(sig_buf)
        print(" {:>6g}".format(sig_mean), end="")
    if "bias" in  quick_fields:
        bias_mean = np.median(bias_buf)
        print(" {:>5g}".format(bias_mean), end="")
    if "signal" in  quick_fields:
        signal = sig_mean - bias_mean
        print(" {:>6g}".format(signal), end="")
    if "noise" in  quick_fields:
        (nrows, ncols) = np.shape(bias_buf)
        print(" {:>7.4g}".format(np.std(
            bias_buf[int(nrows/2-nrows/20):int(nrows/2+nrows/20), 3:ncols-2])), end="")
    if "adu/sec" in  quick_fields and expt > 0:
        print("{:>8.2f}".format(float(signal)/expt), end="")
    if "eper:s-cte" in  quick_fields:
        debugmsg = "s-cte------------------"
        logging.debug(debugmsg)
        (nrows, ncols) = np.shape(sig_buf)
        nsig_cols = ncols
        #- define region to measure signal used in cte calc
        y1 = int(nrows/2-nrows/10)
        y2 = int(nrows/2+nrows/10)
        x0 = ncols-ncols/20
        x1 = ncols-1
        debugmsg = "s_n=sig_buf[{}:{},{}:{}]".format(y1,y2,x0,x1)
        logging.debug(debugmsg)
        s_n = sig_buf[y1:y2, x0:x1]
        (nrows, ncols) = np.shape(bias_buf)
        l_ncols = 3
        bias_mean = np.mean(bias_buf[y1:y2,l_ncols:ncols])
        debugmsg = "using bias_buf[{}:{},{}:{}]".format(y1,y2,l_ncols,ncols)
        logging.debug(debugmsg)
        l_n = np.mean(s_n) - bias_mean
        debugmsg = "l_n={:>10.6g}".format(l_n)
        logging.debug(debugmsg)
        y1 = int(nrows/2-nrows/10)
        y2 = int(nrows/2+nrows/10)
        x0 = 0
        x1 = l_ncols
        debugmsg = "b_n=bias_buf[{}:{},{}:{}]".format(y1,y2,x0,x1)
        logging.debug(debugmsg)
        b_n = bias_buf[y1:y2, x0:x1]
        l_nn = np.mean((b_n - bias_mean).sum(axis=1))
        debugmsg = "l_nn={:>10.6g}".format(l_nn)
        logging.debug(debugmsg)
        if l_n > 0.0:
            eper = 1 - (l_nn / (nsig_cols * l_n))
            print(" {:>8.6g}".format(eper), end="")
    if "eper:p-cte" in  quick_fields:
        debugmsg = "p-cte------------------"
        logging.debug(debugmsg)
        (nrows, ncols) = np.shape(p_bias_buf)
        l_nrows = 3 # number of overscan rows used to determing cte
        #- define region to measure bias used in cte calc
        x1 = int(ncols/2-ncols/10)
        x2 = int(ncols/2+ncols/10)
        y0 = l_nrows
        y1 = nrows
        p_bias_mean = np.mean(p_bias_buf[y0:y1,x1:x2])
        debugmsg = "p_bias_mean=mean(p_bias_buf[{}:{},{}:{}])".format(y0,y1,x1,x2)
        logging.debug(debugmsg)
        debugmsg = "p_bias_mean={:>10.6g}".format(p_bias_mean)
        logging.debug(debugmsg)
        #- define region to measure signal used in cte calc
        (nrows, ncols) = np.shape(sig_buf)
        nsig_rows = nrows
        x1 = int(ncols/2-ncols/10)
        x2 = int(ncols/2+ncols/10)
        y0 = nrows-100
        y1 = nrows-1
        debugmsg = "s_n=sig_buf[{}:{},{}:{}]".format(y0,y1,x1,x2)
        logging.debug(debugmsg)
        s_n = sig_buf[y0:y1, x1:x2]
        l_n = np.mean(s_n) - p_bias_mean
        debugmsg = "l_n={:>10.6g}".format(l_n)
        logging.debug(debugmsg)
        x1 = int(ncols/2-ncols/10)
        x2 = int(ncols/2+ncols/10)
        y0 = 0
        y1 = l_nrows
        debugmsg = "shape(p_bias_buf)={}".format(np.shape(p_bias_buf))
        logging.debug(debugmsg)
        debugmsg = "b_n=p_bias_buf[{}:{},{}:{}]".format(x1,x2,y0,y1)
        logging.debug(debugmsg)
        b_n = p_bias_buf[y0:y1, x1:x2]
        debugmsg = "shape(b_n)={}".format(np.shape(b_n))
        logging.debug(debugmsg)
        #l_nn = np.median(b_n) - p_bias_mean
        l_nn = np.mean((b_n - bias_mean).sum(axis=0))
        debugmsg = "l_nn={:>10.6g}".format(l_nn)
        logging.debug(debugmsg)
        (nrows, ncols) = np.shape(sig_buf)
        if l_n > 0.0:
            eper = 1 - (l_nn / (nrows * l_n))
            print(" {:>8.6g}".format(eper), end="")
    #---------
    if "tearing" in  quick_fields:
        debugmsg = "tearing check----------"
        logging.debug(debugmsg)
    #- column-1 into an array arr1
    #- column-2..N into 2nd array with dimension N-1 x Ncols arr2
    #- take median of 2nd array to make 1-D: arr3
    #- find stddev of arr3
    #- form (arr3 - arr1)/stddev as arr4
    #- find the first value of index "j" in sorted(arr4) less than 1.0
    #- report out (len(arr4)-j)/len(arr4) to 1 digit as tearing where
        #- this represents the fraction of the array less than 1.0
        #- left side
        arr3 = np.median(sig_buf[:,2:40], axis=1)
        arr4 = (arr3 - sig_buf[:,0])/np.std(arr3)
        tm = (1.0*np.size(arr4) - np.searchsorted(arr4, 1.0))/np.size(arr4)
        print("{:>4.1f}".format(tm), end="")
        #- right side
        arr3 = np.median(sig_buf[:,-40:-2], axis=1)
        arr4 = (arr3 - sig_buf[:,-0])/np.std(arr3)
        tm = (1.0*np.size(arr4) - np.searchsorted(arr4, 1.0))/np.size(arr4)
        print("{:>4.1f}".format(tm), end="")
    #---------
    if "dipoles" in  quick_fields:
        debugmsg = "dipoles check----------"
        logging.debug(debugmsg)
    #- region to work on is sig_buf, say 200 rows near top
    #- transpose to column order
    #- find sigma-clipped mean, median and stdev
    #- subtract mean from array
    #- divide the array by sigma
    #- go through array finding pixel pairs of differing sign
    #- and where |A(n)-A(n+1)| > 6
    #- add one to counter each time such a pair is found
    #- print out the % of pixels occupied by dipoles
        (nrows, ncols) = np.shape(sig_buf)
        arr1 = sig_buf[-nrows/10:-1,:] #- use top 10% of array
        debugmsg = "using subarray [{}:{},:]".format(-nrows/10,-1)
        logging.debug(debugmsg)
        arr2 = arr1.flatten('F') #- flatten to 1d in column order
        avg2, med2, std2 = stats.sigma_clipped_stats(arr2)
        debugmsg = "clipped stats: avg:{:>.3g} med:{} stdev:{:>.3g}".format(avg2, med2, std2)
        logging.debug(debugmsg)
        arr3 = (arr2 - avg2)/std2
        ndipole = 0
        for i in range(0, np.size(arr3) - 1):
            if (np.sign(arr3[i+1] * arr3[i]) == -1) and abs(arr3[i+1] - arr3[i]) > 5:
                ndipole += 1
        debugmsg = "dipole count = {}".format(ndipole)
        logging.debug(debugmsg)
        print("{:>8.2f}".format(100.0*float(2*ndipole)/(np.size(arr1))))
    print #-- newline
    ncalls() #-- track call count, acts like static variable

def ncalls():
    """maintain a counter
    """
    ncalls.counter += 1

if __name__ == '__main__':
    main()
