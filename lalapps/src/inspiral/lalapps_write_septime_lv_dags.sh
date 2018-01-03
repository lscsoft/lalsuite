#!/bin/bash  

################################################################################
# get needed options from ini file

month_gps_time=`cat write_ifar_scripts_lv.ini | grep 'month_gps_time' | awk '{print $3}'`
month_duration=`cat write_ifar_scripts_lv.ini | grep 'month_duration' | awk '{print $3}'`
cat=`cat write_ifar_scripts_lv.ini | grep 'cat' | awk '{print $3}'`

septime_path=`cat write_ifar_scripts_lv.ini | grep 'septime_path' | awk '{print $3}'`

log_path=`cat write_ifar_scripts_lv.ini | grep 'log_path' | awk '{print $3}'`
condor_priority=`cat write_ifar_scripts_lv.ini | grep 'condor_priority' | awk '{print $3}'`

hipe_cache=`cat write_ifar_scripts_lv.ini | grep 'hipe_cache' | awk '{print $3}'`

#Print options out to screen for verification
echo "Options used are:"
echo "  month_gps_time = ${month_gps_time}"
echo "  month_duration = ${month_duration}"
echo "  cat = ${cat}"
echo "  septime_path = ${septime_path}"
echo "  log_path = ${log_path}"
echo "  condor_priority = ${condor_priority}"
echo "  hipe_cache = ${hipe_cache}"
echo

# These shouldn't need changing

h1_veto_file='/scratch2/jclayton/runlvtag5/thirdstage/866088014-868721414/segments/H1-COMBINED_'${cat}'_VETO_SEGS-'${month_gps_time}'-'${month_duration}'.txt'
h2_veto_file='/scratch2/jclayton/runlvtag5/thirdstage/866088014-868721414/segments/H2-COMBINED_'${cat}'_VETO_SEGS-'${month_gps_time}'-'${month_duration}'.txt'
l1_veto_file='/scratch2/jclayton/runlvtag5/thirdstage/866088014-868721414/segments/L1-COMBINED_'${cat}'_VETO_SEGS-'${month_gps_time}'-'${month_duration}'.txt'
v1_veto_file='/scratch2/jclayton/runlvtag5/thirdstage/866088014-868721414/segments/V1-COMBINED_'${cat}'_VETO_SEGS-'${month_gps_time}'-'${month_duration}'.txt'
# don't touch anything below here
################################################################################

#get veto segment list
/bin/echo -n "Generating veto segment list... "
ligolw_segments --name vetoes --output vetoes_${cat}.xml.gz \
  --insert-from-segwizard H1=${h1_veto_file} \
  --insert-from-segwizard H2=${h2_veto_file} \
  --insert-from-segwizard L1=${l1_veto_file} \
  --insert-from-segwizard V1=${v1_veto_file}
echo "done."

#generate zero-lag dag
echo "Generating septime_zero_lag.dag and .sub files... "
num_thincas=`grep THINCA_SECOND.*FULL_DATA_${cat} ${hipe_cache} | awk '{print $5}' | sed s+file://localhost++g | wc -l` 
thinca_idx=1

if [ 1 ]; then
  for file in `grep THINCA_SECOND.*FULL_DATA_${cat} ${hipe_cache} | awk '{print $5}' | sed s+file://localhost++g`; do
    echo -ne "processing ${thinca_idx} / ${num_thincas}\r" >&2
    thinca_idx=$(( ${thinca_idx} + 1 ))
    infile=`basename $file`
    job_name=`echo $infile | awk 'gsub("THINCA","SEPTIME")'`
    starttime=`echo $infile | awk 'gsub("-"," ") {print $3}'`
    duration=`echo $infile | sed 's/\./ /g' | awk 'gsub("-"," ") {print $4}'`
    endtime=$(($starttime + $duration))
    combo=`echo $infile | awk 'gsub("-"," ") {print $1}'`
    
    if [ $combo == H1H2L1V1 ] ; then
      triggers="--h1-triggers --h2-triggers --l1-triggers --v1-triggers"    
    elif [ $combo == H1H2L1 ] ; then
      triggers="--h1-triggers --h2-triggers --l1-triggers"
    elif [ $combo == H1H2V1 ] ; then
      triggers="--h1-triggers --h2-triggers --v1-triggers"
    elif [ $combo == H1L1V1 ] ; then
      triggers="--h1-triggers --l1-triggers --v1-triggers"
    elif [ $combo == H2L1V1 ] ; then
      triggers="--h2-triggers --l1-triggers --v1-triggers"
    elif [ $combo == H1H2 ] ; then
      triggers="--h1-triggers --h2-triggers"
    elif [ $combo == H1L1 ] ; then
      triggers="--h1-triggers --l1-triggers"
    elif [ $combo == H2L1 ] ; then
      triggers="--h2-triggers --l1-triggers"
    elif [ $combo == H1V1 ] ; then
      triggers="--h1-triggers --v1-triggers"
    elif [ $combo == H2V1 ] ; then
      triggers="--h2-triggers --v1-triggers"
    elif [ $combo == L1V1 ] ; then
      triggers="--l1-triggers --v1-triggers"
      
    fi

    echo "JOB $job_name septime_zero_lag.septime.sub"
    echo "RETRY $job_name 3"
    echo "VARS $job_name macroinfile=\"$file\" macrotriggers=\"$triggers\" macrogpsstarttime=\"$starttime\" macrogpsendtime=\"$endtime\""
    echo "CATEGORY $job_name septime"
    echo "## JOB $job_name requires input file $infile"
  done
  echo "MAXJOBS septime 20"
fi > septime_zero_lag.dag

if [ 1 ]; then
  echo "universe = vanilla"
  echo "executable = ${septime_path}"
  echo "arguments = --thinca \$(macroinfile) \$(macrotriggers) --veto-file vetoes_${cat}.xml.gz --output-dir septime_files/${cat} --gps-start-time \$(macrogpsstarttime) --gps-end-time \$(macrogpsendtime)"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/septime-\$(cluster)-\$(process).err"
  echo "output = logs/septime-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi > septime_zero_lag.septime.sub
echo -e "\n...done."

#generate injection dag
echo "Generating septime_injection.dag and .sub files... "
num_thincas=`grep THINCA_SECOND.*INJ_${cat} ${hipe_cache} | awk '{print $5}' | sed s+file://localhost++g | wc -l`
thinca_idx=1

if [ 1 ]; then
  for injstring in BBHLININJ BBHLOGINJ BBHSPINLININJ BBHSPINLOGINJ BNSLININJ BNSLOGINJ BNSSPINLININJ BNSSPINLOGINJ NSBHLININJ NSBHLOGINJ NSBHSPINLININJ NSBHSPINLOGINJ BNS001INJ; do
    for file in `grep THINCA_SECOND.*${injstring}_${cat} ${hipe_cache} | awk '{print $5}' | sed s+file://localhost++g`; do
      echo -ne "processing ${thinca_idx} / ${num_thincas}\r" >&2
      thinca_idx=$(( ${thinca_idx} + 1 ))
      infile=`basename $file`
      job_name=`echo $infile | awk 'gsub("THINCA","SEPTIME")'`
      starttime=`echo $infile | awk 'gsub("-"," ") {print $3}'`
      duration=`echo $infile | sed 's/\./ /g' | awk 'gsub("-"," ") {print $4}'`
      endtime=$(($starttime + $duration))
      combo=`echo $infile | awk 'gsub("-"," ") {print $1}'`
      if [ $combo == H1H2L1V1 ] ; then
        triggers="--h1-triggers --h2-triggers --l1-triggers --v1-triggers"
      elif [ $combo == H1H2L1 ] ; then
        triggers="--h1-triggers --h2-triggers --l1-triggers"
      elif [ $combo == H1H2V1 ] ; then
        triggers="--h1-triggers --h2-triggers --v1-triggers"
      elif [ $combo == H1L1V1 ] ; then
        triggers="--h1-triggers --l1-triggers --v1-triggers"
      elif [ $combo == H2L1V1 ] ; then
        triggers="--h2-triggers --l1-triggers --v1-triggers"
      elif [ $combo == H1H2 ] ; then
        triggers="--h1-triggers --h2-triggers"
      elif [ $combo == H1L1 ] ; then
        triggers="--h1-triggers --l1-triggers"
      elif [ $combo == H2L1 ] ; then
        triggers="--h2-triggers --l1-triggers"
      elif [ $combo == H1V1 ] ; then
        triggers="--h1-triggers --v1-triggers"
      elif [ $combo == H2V1 ] ; then
        triggers="--h2-triggers --v1-triggers"
      elif [ $combo == L1V1 ] ; then
        triggers="--l1-triggers --v1-triggers"
      fi

      echo "JOB $job_name septime_${injstring}.septime.sub"
      echo "RETRY $job_name 3"
      echo "VARS $job_name macroinfile=\"$file\" macrotriggers=\"$triggers\" macrogpsstarttime=\"$starttime\" macrogpsendtime=\"$endtime\""
      echo "CATEGORY $job_name septime"
      echo "## JOB $job_name requires input file $infile"
    done
  done
  echo "MAXJOBS septime 20"
fi > septime_injection.dag

for injstring in BBHLININJ BBHLOGINJ BBHSPINLININJ BBHSPINLOGINJ BNSLININJ BNSLOGINJ BNSSPINLININJ BNSSPINLOGINJ NSBHLININJ NSBHLOGINJ NSBHSPINLININJ NSBHSPINLOGINJ BNS001INJ; do 
  if [ 1 ]; then
    echo "universe = vanilla"
    echo "executable = ${septime_path}"
    echo "arguments = --thinca \$(macroinfile) \$(macrotriggers) --veto-file vetoes_${cat}.xml.gz --output-dir septime_files/${injstring}_${cat} --gps-start-time \$(macrogpsstarttime) --gps-end-time \$(macrogpsendtime)"
    echo "getenv = True"
    echo "log = " `mktemp -p ${log_path}`
    echo "error = logs/septime-\$(cluster)-\$(process).err"
    echo "output = logs/septime-\$(cluster)-\$(process).out"
    echo "notification = never"
    echo "priority=${condor_priority}"
    echo "queue 1"
  fi > septime_${injstring}.septime.sub
done
echo -e "\n...done."

#generate time-slide dag
echo "Generating septime_slide.dag and .sub files..."
num_thincas=`grep THINCA_SLIDE_SECOND.*FULL_DATA_CAT_3 ${hipe_cache} | awk '{print $5}' | sed s+file://localhost++g | wc -l` 
thinca_idx=1


if [ 1 ]; then
  for file in `grep THINCA_SLIDE_SECOND.*FULL_DATA_CAT_3 ${hipe_cache} | awk '{print $5}' | sed s+file://localhost++g`; do
    echo -ne "processing ${thinca_idx} / ${num_thincas}\r" >&2
    thinca_idx=$(( ${thinca_idx} + 1 ))
    infile=`basename $file`
    job_name=`echo $infile | awk 'gsub("THINCA","SEPTIME")'`
    starttime=`echo $infile | awk 'gsub("-"," ") {print $3}'`
    duration=`echo $infile | sed 's/\./ /g' | awk 'gsub("-"," ") {print $4}'`
    endtime=$(($starttime + $duration))
    combo=`echo $infile | awk 'gsub("-"," ") {print $1}'`
    
    if [ $combo == H1H2L1V1 ] ; then
      triggers="--h1-triggers --h2-triggers --l1-triggers --v1-triggers"
    elif [ $combo == H1H2L1 ] ; then
      triggers="--h1-triggers --h2-triggers --l1-triggers"
    elif [ $combo == H1H2V1 ] ; then
      triggers="--h1-triggers --h2-triggers --v1-triggers"
    elif [ $combo == H1L1V1 ] ; then
      triggers="--h1-triggers --l1-triggers --v1-triggers"
    elif [ $combo == H2L1V1 ] ; then
      triggers="--h2-triggers --l1-triggers --v1-triggers"
    elif [ $combo == H1H2 ] ; then
      triggers="--h1-triggers --h2-triggers"
    elif [ $combo == H1L1 ] ; then
      triggers="--h1-triggers --l1-triggers"
    elif [ $combo == H2L1 ] ; then
      triggers="--h2-triggers --l1-triggers"
    elif [ $combo == H1V1 ] ; then
      triggers="--h1-triggers --v1-triggers"
    elif [ $combo == H2V1 ] ; then
      triggers="--h2-triggers --v1-triggers"
    elif [ $combo == L1V1 ] ; then
      triggers="--l1-triggers --v1-triggers"
    

    fi

    echo "JOB $job_name septime_slide.septime.sub"
    echo "RETRY $job_name 3"
    echo "VARS $job_name macroinfile=\"$file\" macrotriggers=\"$triggers\" macrogpsstarttime=\"$starttime\" macrogpsendtime=\"$endtime\""
    echo "CATEGORY $job_name septime_slide"
    echo "## JOB $job_name requires input file $infile"
  done
  echo "MAXJOBS septime_slide 40"
fi> septime_slide.dag

# Slides are hard coded but looks correct. Make sure that below the time shifts and number of time slides are correct.

if [ 1 ]; then
  echo "universe = vanilla"
  echo "executable = ${septime_path}"
  echo "arguments = --thinca \$(macroinfile) \$(macrotriggers) --veto-file vetoes_${cat}.xml.gz --output-dir septime_files/${cat} --gps-start-time \$(macrogpsstarttime) --gps-end-time \$(macrogpsendtime) --h1-slide 0 --h2-slide 0 --l1-slide 5 --v1-slide 15 --num-slides 50 --write-antime-file"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/septime_slide-\$(cluster)-\$(process).err"
  echo "output = logs/septime_slide-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi > septime_slide.septime.sub
echo -e "\n...done."

#setup directory structure
if [ ! -d septime_files ] ; then
  mkdir septime_files
fi
if [ ! -d septime_files/${cat} ] ; then
  mkdir septime_files/${cat}
fi
for string in BBHLININJ BBHLOGINJ BBHSPINLININJ BBHSPINLOGINJ BNSLININJ BNSLOGINJ BNSSPINLININJ BNSSPINLOGINJ NSBHLININJ NSBHLOGINJ NSBHSPINLININJ NSBHSPINLOGINJ BNS001INJ; do
  if [ ! -d septime_files/${string}_${cat} ] ; then
    mkdir septime_files/${string}_${cat}
  fi
done
if [ ! -d logs ] ; then
  mkdir logs
fi

echo "******************************************************"
echo "  Now run: condor_submit_dag septime_zero_lag.dag"
echo "      and: condor_submit_dag septime_slide.dag"
echo "      and: condor_submit_dag septime_injection.dag"
echo "  These dags can be run simutaneously."
echo "******************************************************"

