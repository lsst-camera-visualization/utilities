#!/bin/bash
#
# Get the default ts71G snapshot report
#-------------------------------------------------------------------------
echo \$0=$0
rp=$(/usr/bin/realpath $0)
bindir=${rp%/*}
basedir=${bindir%/*}
echo \# ${basedir}/trendutils/trendmonitor.py ${basedir}/data/ts71g.report
${basedir}/trendutils/trendmonitor.py ${basedir}/data/ts71g.report
