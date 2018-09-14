#!/bin/bash
#
# Get the default pump-plate snapshot report
#-------------------------------------------------------------------------
rp=$(/usr/bin/realpath $0)
bindir=${rp%/*}
basedir=${bindir%/*}
echo \# ${basedir}/trendutils/trendmonitor.py ${basedir}/data/pump-plate.report
${basedir}/trendutils/trendmonitor.py ${basedir}/data/pump-plate.report
