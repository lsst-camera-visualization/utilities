#!/bin/bash
#
# script to dump out trending database channels for a given subsystem
#----------------------------

# Print out error message when called commandline error
function usage {
cat <<-EOM
    Usage: ${0##*/} [Options] SubSystemName
    where "SubSystemName" is like "ts5-raft"
    Options:
    -h  (print this help msg)
EOM
exit 1
}

#-- process commandline arguments
if [ $# -ne 1 ]; then
   usage
fi
ssname=$1

url="http://lsst-mcm.slac.stanford.edu:8080/rest/data/dataserver/listchannels"
/usr/bin/wget --quiet -O - ${url} |\
/lnfs/lsst/dh/software/centos7-gcc48/anaconda/py-2.7-4.3.14/bin/xmllint --format - |\
/usr/bin/gawk '\
   BEGIN {println=0;channelstart=0};\
   /<datachannel>/ {channelstart=1; next;};\
   /<path>/ {next;};\
   /<\/path>/ {next;};\
   /<pathelement>/ {if(channelstart==1 && index($0,"'${ssname}'")) {printf("\n"); println=1;}};\
   /<pathelement>/ {if(println) {gsub("<pathelement>","");gsub("</pathelement>","");print $0;}};\
   /<id>/ {if(println==1) {gsub("<id>","");gsub("</id>","");print $0}};\
   /<\/datachannel>/ {channelstart=0; println=0;};' |\
/usr/bin/gawk '\
   /^$/ {printf("\n");next;}; {printf("%s ",$0);};'
