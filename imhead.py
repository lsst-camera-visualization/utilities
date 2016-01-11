#!/usr/bin/env python
#
#
#------------------------
import sys
import argparse
import numpy as np
from astropy.io import fits

#- process arguments
parser = argparse.ArgumentParser(description="Print headers (Primary+Image)")
parser.add_argument("fitsfile", help="input fits file")
parser.add_argument("--info", action='store_true', 
        help="print the info() table summarizing file contents")
parser.add_argument("--hduname", help="dump only the named extension header")
parser.add_argument("--hduindex", type=int, help="dump only id'd extension header")
parser.add_argument("--all", type=int, help="print all extensions")
args = parser.parse_args()


seglist = {}        # dict to hold segment names old and new hduList indexes
otherlist = {}      # dict to hold non-image HDU's (except for primary)
hduList = fits.open(args.fitsfile)
#- just print the image info with indexes, names, sizes
if(args.info):
    hduList.info()
#- print single hdu by name
elif(args.hduname):
    index = hduList.index_of(args.hduname)
    hdu = fits.ImageHDU(header=hduList[index].header)
    print "#--------{}---------".format(args.hduname)
    hdu.header.totextfile(sys.stdout)
#- print single hdu by index
elif(args.hduindex):
    hdu = fits.ImageHDU(header=hduList[args.hduindex].header)
    print "#--------extension {}---------".format(args.hduindex)
    hdu.header.totextfile(sys.stdout)
#- print primary and Image headers, others optionally
else:
    #- build dicts
    for hdu in hduList:
        index = hduList.index_of(hdu.name)
        if (isinstance(hdu, fits.ImageHDU)):
            seglist[hdu.name] = index
        elif (isinstance(hdu, fits.PrimaryHDU)):
            pindex = index
        else:
            otherlist[hdu.name] = index
    #- print the primary
    hdu = fits.ImageHDU(header=hduList[pindex].header)
    print "#--------{}---------".format(hdu.name)
    hdu.header.totextfile(sys.stdout)
    #- print the Image headers
    for name, index in sorted(seglist.iteritems()):
        hdu = fits.ImageHDU(header=hduList[index].header)
        print "#--------{}---------".format(hdu.name)
        hdu.header.totextfile(sys.stdout)
    #- print the other headers
    if(args.all):
        ids = list(otherlist.values())
        for index in sorted(ids):
            hdu = fits.ImageHDU(header=hduList[index].header)
            print "#--------{}---------".format(hdu.name)
            hdu.header.totextfile(sys.stdout)


hduList.close()
sys.exit()
