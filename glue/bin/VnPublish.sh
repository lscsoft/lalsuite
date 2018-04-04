#!/bin/bash

########################################################
# Usage:
# sh VnPublish.sh gps_end_time version_number
# for example:
# sh VnPublish.sh 956000000 2
#


# This script calls VnPublish.py to execute gap checking
#########################################################


clear
# Read in start time from file to use in cache file
. last_et
cache_start_time=$last_end_time
cache_end_time=$1
n=$2  # version to publish
echo "gps times to use in cache file:  ["$cache_start_time", "$cache_end_time")"

# make a temp dir
tmp_dir=`mktemp -d`
echo "temp dir is: "$tmp_dir
cd ${tmp_dir}

# retrieve version n data from frame file
# retrieve LHO version n data
ligo_data_find -r ldr.ligo.caltech.edu -s $cache_start_time -e $cache_end_time -o H -t H1_LDAS_C02_L2 -l -u file > H-H1_LDAS_C02_L2-${cache_start_time}-$((${cache_end_time} - ${cache_start_time})).cache || exit 1

ligolw_fr_to_science --cache-file=H-H1_LDAS_C02_L2-${cache_start_time}-$((${cache_end_time} - ${cache_start_time})).cache --ifo=H1 --segment-version=${n} --type=H1_LDAS_C02_L2 || exit 1


# retreive LLO version n data
ligo_data_find -r ldr.ligo.caltech.edu -s $cache_start_time -e $cache_end_time -o L -t L1_LDAS_C02_L2 -l -u file > L-L1_LDAS_C02_L2-${cache_start_time}-$((${cache_end_time} - ${cache_start_time})).cache || exit 1

ligolw_fr_to_science --cache-file=L-L1_LDAS_C02_L2-${cache_start_time}-$((${cache_end_time} - ${cache_start_time})).cache --ifo=L1 --segment-version=${n} --type=L1_LDAS_C02_L2 || exit 1




#--------------------------------------------------------------------------------
# Compare segment_summary table between H1:DMT-SCIENCE:1 and H1:DMT-SCIENCE:n
echo "Checking if H1:DMT-SCIENCE:1 segments are a subset of H1:DMT-SCIENCE:${n} segments ..."
# 4.1 Read in version 1 data from the segment database
ligolw_segment_query --segment-url https://segdb.ligo.caltech.edu --query-segments -s $cache_start_time -e $cache_end_time --include-segments 'H1:DMT-SCIENCE:1' > v1h1_sci_segs.xml || exit 1

ligolw_segment_intersect --segment -i H1:RESULT:1,H1:DMT-SCIENCE:$n v1h1_sci_segs.xml H-H1_LDAS_C02_L2_SEGMENTS*.xml > v1vn_intersect.xml

v1vn_diff=`ligolw_segment_diff --segment -i H1:RESULT:1,:RESULT:1 v1h1_sci_segs.xml v1vn_intersect.xml | ligolw_print -t segment -c start_time -c end_time`
if [ "$v1vn_diff" == "" ]; then
   echo "test passed"
else
   echo "test failed. Discrepancies are: $v1vn_diff"
   exit 1
fi


#--------------------------------------------------------------------------------
# check if gaps in Vn summary file equal to LDAS allowed gap list
echo
echo "Checking if gaps of summary intervals for H1:DMT-SCIENCE:$n are a subset of the LDAS allowed gaps ... "

# retrieve LDAS allowed gaps from the database:
ligolw_segment_query --segment-url https://segdb-dev.phy.syr.edu --query-segments --gps-start-time $cache_start_time --gps-end-time $cache_end_time --include-segments H1:DCH-MISSING_LDAS_C02_L2:1 > ldas_gaps.xml || exit 1

result=($(python -c "import sys; sys.path.append('/archive/home/piwei/hoft/test'); import VnPublish; VnPublish.vn_ldas_gap_check('`echo H-H1_LDAS_C02_L2_SEGMENTS_V2-*.xml`','H1:DMT-SCIENCE:$n','ldas_gaps.xml')")) || exit 1

if [ "$result" != "" ]; then
   echo $result;
   exit 1
else 
     echo "test passed"
fi



#--------------------------------------------------------------------------------
# if files passed all the checking, insert Vn data into the segment database
echo
echo "Inserting V$n data into the segment database ..."
ldbdc --segment-url https://segdb-dev.phy.syr.edu --dmt-segments --insert H-H1_LDAS_C02_L2_SEGMENTS_V2-*.xml
error=`echo $?`
if [ $error -ne 0 ]; then
   echo "Error inserting LHO data into the segment database"
   exit 1
fi

ldbdc --segment-url https://segdb-dev.phy.syr.edu --dmt-segments --insert L-L1_LDAS_C02_L2_SEGMENTS_V2-*.xml
error=`echo $?`
if [ $error -ne 0 ]; then
   echo "Error inserting LLO data into the segment database"
   exit 1
fi


#  update end_time to the .ini file for next time use
`echo "last_end_time=$cache_end_time" > /archive/home/piwei/hoft/H-H1_LDAS_C02/last_et` || exit 1

#  copy cache file and xml file to central location
cp *.cache /archive/home/piwei/hoft/H-H1_LDAS_C02 || exit 1
cp H-H1_LDAS_C02_L2_SEGMENTS_V2-*.xml /archive/home/piwei/hoft/H-H1_LDAS_C02 || exit 1
cp L-L1_LDAS_C02_L2_SEGMENTS_V2-*.xml /archive/home/piwei/hoft/H-H1_LDAS_C02 || exit 1



# Remove the tmp dir
cd
rm -rf ${tmp_dir}


echo "Done!"
exit 0
