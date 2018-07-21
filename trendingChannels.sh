#!/bin/bash
#
# script to dump out trending database channels for a given subsystem
#----------------------------
XMLLINT=$(which xmllint)
WGET=$(which wget)
GAWK=$(which gawk)

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

if $(hostname -I | grep -q 134.79) ; then
    trending_server=lsst-mcm.slac.stanford.edu
else
    #- expect to have an ssh tunnel set up
    trending_server=localhost
fi
BASE_URL="http://${trending_server}:8080/rest/data/dataserver"

url="${BASE_URL}/listchannels"
${WGET} --quiet -O - ${url} |\
${XMLLINT} --format - |\
${GAWK} '\
   BEGIN {println=0;channelstart=0};\
   /<datachannel>/ {channelstart=1; next;};\
   /<path>/ {next;};\
   /<\/path>/ {next;};\
   /<pathelement>/ {if(channelstart==1 && index($0,"'${ssname}'")) {printf("\n"); println=1;}};\
   /<pathelement>/ {if(println) {gsub("<pathelement>","");gsub("</pathelement>","");print $0;}};\
   /<id>/ {if(println==1) {gsub("<id>","");gsub("</id>","");print $0}};\
   /<\/datachannel>/ {channelstart=0; println=0;};' |\
${GAWK} '\
   /^$/ {printf("\n");next;}; {printf("%s ",$0);};'
