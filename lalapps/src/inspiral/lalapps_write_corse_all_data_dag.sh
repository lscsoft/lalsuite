#!/bin/bash 

################################################################################
# edit these appropriately

month_gps_time='847555570'
month_duration='2419200'
cat='CAT_3'

coire_path='/home/cdcapano/local/s5_2yr_lowcbc_20080829/bin/lalapps_coire'
corse_path='/home/cdcapano/local/s5_2yr_lowcbc_20080829/bin/lalapps_corse'

log_path='/usr1/cdcapano/log'
condor_priority='20'

# don't touch anything below here
################################################################################

#get septime txt files organized by combo and run add_septime to add up times
echo "Getting SEPTIME_SLIDE*.txt file names and writing to cache files..."
pushd septime_files/${cat}/ > /dev/null
for combo in H1H2L1 H1H2 H1L1 H2L1; do
  echo -n "  Getting ${cat}/${combo}-SEPTIME_SLIDE*.txt files..."
  for file in ${combo}-SEPTIME_SLIDE_H*.txt; do
    echo ${file}
  done > ../septime_${cat}_${combo}_times.cache
  echo " done."
done
popd > /dev/null
echo "...done."

#generate corse dag
/bin/echo -n "Generating corse_all_data.dag and .sub files..."
if [ 1 ]; then
  #write add_septime job
  for combo in H1L1 H2L1 H1H2L1; do
    infile=septime_${cat}_${combo}_times.cache
    outfile=${combo}_V3_${cat}.txt
    echo "JOB $outfile corse_all_data.add_septime.sub"
    echo "RETRY $outfile 0"
    echo "VARS $outfile macrooutfile=\"$outfile\" macroinfile=\"$infile\""
    echo "CATEGORY $outfile add_septime"
    echo
  done
  for mass in mchirp_2_8 mchirp_8_17 mchirp_17_35; do
    #write corse jobs for double-in_double files
    for combo in H1L1 H2L1; do
      zerofile="second_coire_files/${mass}/${combo}-SECOND_COIRE_${cat}_${combo}-${month_gps_time}-${month_duration}.xml.gz"
      slidefile="second_coire_files/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_${combo}-${month_gps_time}-${month_duration}.xml.gz"
      outfile="corse_all_data_files/${mass}/${combo}-CORSE_ALL_DATA_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      timeanalyzedfile="septime_files/${combo}_V3_${cat}.txt"
      echo "JOB $outfile corse_all_data.corse.sub"
      echo "RETRY $outfile 1"
      echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\""
      echo "CATEGORY $outfile corse"
      parent_file="${combo}_V3_${cat}.txt"
      echo "PARENT $parent_file CHILD $outfile"
      echo
    done
    #write corse jobs for double-in_triple and triple-in_triple files
    for combo in H1L1 H2L1 H1H2L1; do
      zerofile="second_coire_files/${mass}/${combo}-SECOND_COIRE_H1H2L1-${month_gps_time}-${month_duration}.xml.gz"
      slidefile="second_coire_files/${mass}/${combo}-SECOND_COIRE_SLIDE_H1H2L1-${month_gps_time}-${month_duration}.xml.gz"
      outfile="corse_all_data_files/${mass}/H1H2L1_${combo}-CORSE_ALL_DATA_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      timeanalyzedfile="septime_files/H1H2L1_V3_${cat}.txt"
      echo "JOB $outfile corse_all_data.corse.sub"
      echo "RETRY $outfile 1"
      echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\""
      echo "CATEGORY $outfile corse"
      parent_file="H1H2L1_V3_${cat}.txt"
      echo "PARENT $parent_file CHILD $outfile"
      echo
    done
  done
  echo "MAXJOBS corse 20"
  echo "MAXJOBS add_septime 3"
fi > corse_all_data.dag

#write add_septime.sub file
if [ 1 ] ; then
  echo "universe = vanilla"
  echo "executable = add_septime.py"
  echo "arguments = --input-file septime_files/\$(macroinfile) --output-file septime_files/\$(macrooutfile) --num-slides 50"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/add_septime-\$(cluster)-\$(process).err"
  echo "output = logs/add_septime-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > corse_all_data.add_septime.sub
 
#write corse.sub file
if [ 1 ] ; then
  echo "universe = standard"
  echo "executable = ${corse_path}"
  echo "arguments = --glob-zero \$(macrozerofile) --glob-slide \$(macroslidefile) --output \$(macrooutfile) --data-type all_data --coinc-stat effective_snrsq --num-slides 50 --time-analyzed-file  \$(macrotimeanalyzedfile)"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/corse-\$(cluster)-\$(process).err"
  echo "output = logs/corse-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > corse_all_data.corse.sub
echo " done."

#setup directory structure
/bin/echo -n "Creating corse_all_data_files directory structure..."
if [ ! -d corse_all_data_files ] ; then
  mkdir corse_all_data_files
fi
for mass in mchirp_2_8 mchirp_8_17 mchirp_17_35; do
  if [ ! -d corse_all_data_files/${mass} ] ; then
    mkdir corse_all_data_files/${mass}
  fi
done
echo " done."
echo "******************************************************"
echo "  Now run: condor_submit_dag corse_all_data.dag"
echo "******************************************************"

