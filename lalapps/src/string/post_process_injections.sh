#
# Bash script to do post-processing on output of the injection DAG.
#

. post_process_functions.sh

INJECTIONDIR=/scratch2/xavi/CosmicStrings2006October6Injections

# Build the time slides (only zero-lag for injections).
ligolw_tisi --verbose --instrument H1=0:0:0 --instrument H2=0:0:0 --instrument L1=0:0:0 >time_slides.xml

# Build a LAL cache of all trigger files and pipe to ligolw_cafe.  Need to
# pass glob patterns because there are too many files for the command line.
# Include the injection files, too!  Fun trick:  because the injection
# description files are included in the cache, and each covers a whole
# segment, ligolw_cafe will naturally group files by segment.
lalapps_ll2cache --verbose --description "-" "${INJECTIONDIR}/triggers/*.xml" "${INJECTIONDIR}/HL-INJECTIONS-*.xml" | ligolw_cafe --verbose --time-slides time_slides.xml

# Remove any caches that don't include injection files to work around an
# issue with the DAG
rm -f $(grep -L "HL-" cafe_*.cache)

# Merge files according to ligolw_cafe's prescription.  Include the time
# slides, too.  Note that adding .gz to the end of the output name causes
# it to be gzip compressed.
for f in cafe_*.cache ; do ligolw_add --verbose --input-cache $f time_slides.xml --output ${f/.cache/.xml.gz} ; done

# Trim injections to just those that were made.  (commented out because DAG
# does this correctly)
#ligolw_bucut --verbose --program StringSearch --inj-made-only cafe_*.xml.gz

# Coincidence
condor_run_burca cafe_*.xml.gz
rm -f logs/*

# Injection finding
ligolw_binjfind --verbose --compare bytime cafe_*.xml.gz
