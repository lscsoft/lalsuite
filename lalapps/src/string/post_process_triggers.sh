#
# Bash script to do post processing on the output of the non-injection DAG.
#

. post_process_functions.sh

TRIGGERDIR=/scratch2/xavi/CosmicStrings2006October6

# Build the time slides.
ligolw_tisi --verbose --instrument H1=0:0:0 --instrument H2=0:0:0 --instrument L1=-89.26990816987241548050:89.26990816987241548050:1.78539816339744830961 >time_slides.xml

# Build a LAL cache of all trigger files and pipe to ligolw_cafe.  Need to
# use a glob pattern because there are too many files.
lalapps_ll2cache --verbose --description "-" "${TRIGGERDIR}/triggers/*.xml" | ligolw_cafe --verbose --time-slides time_slides.xml

# Merge trigger files according to ligolw_cafe's prescription (don't forget
# the time slides!).  Note that adding .gz to the end of the output name
# causes it to be gzip compressed.
for f in cafe_*.cache ; do ligolw_add --verbose --input-cache $f time_slides.xml --output ${f/.cache/.xml.gz} ; done

# Coincidence
condor_run_burca cafe_*.xml.gz
rm -f logs/*
