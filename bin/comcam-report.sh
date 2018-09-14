#!/bin/bash
#
# Get the default comcam snapshot report
#-------------------------------------------------------------------------
rp=$(/usr/bin/realpath $0)
bindir=${rp%/*}
basedir=${bindir%/*}
echo \# ${basedir}/trendutils/trendmonitor.py ${basedir}/data/comcam.report
${basedir}/trendutils/trendmonitor.py ${basedir}/data/comcam.report
