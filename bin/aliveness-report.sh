#!/bin/bash
#
# Get the default aliveness-bench snapshot report
#-------------------------------------------------------------------------
rp=$(/usr/bin/realpath $0)
bindir=${rp%/*}
basedir=${bindir%/*}
echo \# ${basedir}/trendutils/trendmonitor.py ${basedir}/data/aliveness-bench.report
${basedir}/trendutils/trendmonitor.py ${basedir}/data/aliveness-bench.report
