#!/bin/bash

################################################################################
# get needed options from ini file

data_type=`cat write_ifar_scripts.ini | grep 'data_type' | awk '{print $3}'`

month_gps_time=`cat write_ifar_scripts.ini | grep 'month_gps_time' | awk '{print $3}'`
month_duration=`cat write_ifar_scripts.ini | grep 'month_duration' | awk '{print $3}'`
cat=`cat write_ifar_scripts.ini | grep 'cat' | awk '{print $3}'`

coire_path=`cat write_ifar_scripts.ini | grep 'coire_path' | awk '{print $3}'`
plotifar_path=`cat write_ifar_scripts.ini | grep 'plotifar_path' | awk '{print $3}'`

log_path=`cat write_ifar_scripts.ini | grep 'log_path' | awk '{print $3}'`
condor_priority=`cat write_ifar_scripts.ini | grep 'condor_priority' | awk '{print $3}'`

#Print options out to screen for verification
echo "Options used are:"
echo "  data_type = ${data_type}"
echo "  month_gps_time = ${month_gps_time}"
echo "  month_duration = ${month_duration}"
echo "  cat = ${cat}"
echo "  coire_path = ${coire_path}"
echo "  plotifar_path = ${plotifar_path}"
echo "  log_path = ${log_path}"
echo "  condor_priority = ${condor_priority}"
echo
################################################################################

#gps_end_time is needed for plotifar (don't need to edit)
gps_end_time=$(( ${month_gps_time} + ${month_duration} ))

#generate dag
/bin/echo -n "Generating ifar_result_${data_type}.dag and .sub files..."
if [ 1 ]; then
  #run plotifar on the individual mass bins
  for mass in mchirp_2_8 mchirp_8_17 mchirp_17_35; do
    for combo in H1H2L1V1 H1H2L1 H1L1V1 H2L1V1 H1H2V1 H1H2 H1L1 H2L1 H1V1 H2V1 L1V1; do
      job_name="${combo}-plotifar_${mass}_${data_type}"
      glob_files="corse_all_data_files/${data_type}/${combo}*CORSE_`echo ${data_type} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
      outpath="ifar_result_files/${data_type}/${mass}"
      time_correct_file="second_coire_files/summ_files_all_data/${combo}-SECOND_COIRE_${cat}_${combo}-${month_gps_time}-${month_duration}.txt"
      user_tag="${mass}-${data_type}"
      echo "JOB $job_name ifar_result_${data_type}.plotifar.sub"
      echo "RETRY $job_name 0"
      echo "VARS $job_name macroglob=\"$glob_files\" macrooutpath=\"$outpath\" macroifotimes=\"$combo\"  macrotcorrfile=\"$time_correct_file\" macrousertag=\"$user_tag\"" 
      echo "CATEGORY $job_name plotifar"
      echo
    done
  done
  # run plotifar on all the mass bins globbed together 
  for combo in H1H2L1V1 H1H2L1 H1L1V1 H2L1V1 H1H2V1 H1H2 H1L1 H2L1 H1V1 H2V1 L1V1; do
    job_name="${combo}-plotifar_ALL_MASSES_${data_type}"
    glob_files="corse_all_data_files/${data_type}/${combo}*CORSE_`echo ${data_type} | tr '[a-z]' '[A-Z]'`_*_${cat}-${month_gps_time}-${month_duration}.xml.gz"
    outpath="ifar_result_files/${data_type}/"
    time_correct_file="second_coire_files/summ_files_all_data/${combo}-SECOND_COIRE_${cat}_${combo}-${month_gps_time}-${month_duration}.txt"
    user_tag="ALL_MASSES-${data_type}"
    echo "JOB $job_name ifar_result_${data_type}.plotifar.sub"
    echo "RETRY $job_name 0"
    echo "VARS $job_name macroglob=\"$glob_files\" macrooutpath=\"$outpath\" macroifotimes=\"$combo\" macrotcorrfile=\"$time_correct_file\" macrousertag=\"$user_tag\""
    echo "CATEGORY $job_name plotifar"
    echo
  done
  echo "MAXJOBS filter_coire 20"
  echo "MAXJOBS coire 20"
  echo "MAXJOBS plotifar 20"
fi > ifar_result_${data_type}.dag

if [ 1 ]; then
  echo "universe = vanilla"
  echo "executable = ${plotifar_path}"
  echo "arguments = --glob \$(macroglob) --output-path \$(macrooutpath) --enable-output --ifo-times \$(macroifotimes) --gps-start-time ${month_gps_time} --gps-end-time ${gps_end_time} --ifan-dist --ifar-dist --plot-uncombined --plot-combined --show-min-bkg --show-max-bkg --show-two-sigma-error --time-correct-file \$(macrotcorrfile) --user-tag \$(macrousertag)"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/plotifar-\$(cluster)-\$(process).err"
  echo "output = logs/plotifar-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > ifar_result_${data_type}.plotifar.sub

#make directory structure
if [ ! -d ifar_result_files ] ; then
 mkdir ifar_result_files
fi
if [ ! -d ifar_result_files/summ_files ] ; then
  mkdir ifar_result_files/summ_files
fi
if [ ! -d ifar_result_files/${data_type} ] ; then
  mkdir ifar_result_files/${data_type}
fi
for mass in mchirp_2_8 mchirp_8_17 mchirp_17_35; do
  if [ ! -d ifar_result_files/${data_type}/${mass} ] ; then
    mkdir ifar_result_files/${data_type}/${mass}
  fi
done
echo " done."
echo "*******************************************************************"
echo "  Now run: condor_submit_dag ifar_result_${data_type}.dag"
echo "*******************************************************************"

