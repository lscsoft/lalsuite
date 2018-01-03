#!/bin/bash 

################################################################################
# get needed options from ini file

cat=`cat write_ifar_scripts_lv.ini | grep 'cat' | awk '{print $3}'`

coire_path=`cat write_ifar_scripts_lv.ini | grep 'coire_path' | awk '{print $3}'`

log_path=`cat write_ifar_scripts_lv.ini | grep 'log_path' | awk '{print $3}'`
condor_priority=`cat write_ifar_scripts_lv.ini | grep 'condor_priority' | awk '{print $3}'`

#Print options out to screen for verification
echo "Options used are:"
echo "  cat = ${cat}"
echo "  coire_path = ${coire_path}"
echo "  log_path = ${log_path}"
echo "  condor_priority = ${condor_priority}"
echo
################################################################################

#get septime zero-lag files
/bin/echo -n "Generating septime file list..."
num_septimes=0
pushd septime_files/${cat} > /dev/null
for file in *SEPTIME_H*xml.gz; do
  echo ${file}
  num_septimes=$(( ${num_septimes} + 1 ))
done > ../septime_${cat}.cache
#popd > /dev/null
#echo "done"



for file in *SEPTIME_L*xml.gz; do
  echo ${file}
  num_septimes=$(( ${num_septimes} + 1 ))
done >> ../septime_${cat}.cache
popd > /dev/null
echo " done."


#generate zero-lag coire dag
echo "Generating zero-lag first_coire.dag and .sub files..."
septime_idx=1
if [ 1 ]; then
  for infile in `cat septime_files/septime_${cat}.cache`; do
    echo -ne "processing ${septime_idx} / ${num_septimes}\r" >&2
    septime_idx=$(( ${septime_idx} + 1))
    outfile=`echo $infile | sed s/SEPTIME/COIRE_${cat}/g`
    echo "JOB $outfile first_coire.coire.sub"
    echo "RETRY $outfile 1"
    echo "VARS $outfile macroinfile=\"$infile\" macrooutfile=\"$outfile\""
    echo "CATEGORY $outfile coire"
    echo "## JOB $outfile requires input file $infile"
  done 
  echo "MAXJOBS coire 20"
fi > first_coire.dag

if [ 1 ]; then
  echo "universe = standard"
  echo "executable = ${coire_path}"
  echo "arguments = --glob septime_files/${cat}/\$(macroinfile) --output first_coire_files/${cat}/\$(macrooutfile) --data-type all_data --coinc-stat effective_snrsq --cluster-time 10000"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/coire-\$(cluster)-\$(process).err"
  echo "output = logs/coire-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > first_coire.coire.sub
echo -e "\n...done."

#get septime injection files
/bin/echo -n "Generating septime file list..."
num_septimes=0
for injstring in BBHLININJ BBHLOGINJ BBHSPINLININJ BBHSPINLOGINJ BNSLININJ BNSLOGINJ BNSSPINLININJ BNSSPINLOGINJ NSBHLININJ NSBHLOGINJ NSBHSPINLININJ NSBHSPINLOGINJ BNS001INJ; do
  pushd septime_files/${injstring}_${cat} > /dev/null
  for file in *SEPTIME_H*xml.gz; do
    echo ${file}
    num_septimes=$(( ${num_septimes} + 1 ))
  done > ../septime_${injstring}_${cat}.cache

  for file in *SEPTIME_L*xml.gz; do
    echo ${file}
    num_septimes=$(( ${num_septimes} + 1 ))
  done >> ../septime_${injstring}_${cat}.cache

  popd > /dev/null
done
echo " done."

#generate injection coire dag
echo "Generating injection first_coire.dag and .sub files..."
septime_idx=1
for injstring in BBHLININJ BBHLOGINJ BBHSPINLININJ BBHSPINLOGINJ BNSLININJ BNSLOGINJ BNSSPINLININJ BNSSPINLOGINJ NSBHLININJ NSBHLOGINJ NSBHSPINLININJ NSBHSPINLOGINJ BNS001INJ; do
  #get HL-INJ file
  for file in `ls ../*inj/HL*${injstring}*`; do
    hlinjfile=$file
  done
  if [ 1 ]; then
    for infile in `cat septime_files/septime_${injstring}_${cat}.cache`; do
      echo -ne "processing ${septime_idx} / ${num_septimes}\r" >&2
      septime_idx=$(( ${septime_idx} + 1))
      outfile=`echo $infile | sed s/SEPTIME/COIRE_${injstring}_${cat}/g`
      echo "JOB $outfile first_coire_${injstring}.coire.sub"
      echo "RETRY $outfile 1"
      echo "VARS $outfile macroinfile=\"$infile\" macrooutfile=\"$outfile\" macroinjectionfile=\"$hlinjfile\" "
      echo "CATEGORY $outfile coire"
      echo "## JOB $outfile requires input file $infile"
    done
  fi
  echo "MAXJOBS coire 200"
done > first_coire_injection.dag

for injstring in BBHLININJ BBHLOGINJ BBHSPINLININJ BBHSPINLOGINJ BNSLININJ BNSLOGINJ BNSSPINLININJ BNSSPINLOGINJ NSBHLININJ NSBHLOGINJ NSBHSPINLININJ NSBHSPINLOGINJ BNS001INJ; do
  if [ 1 ]; then
    echo "universe = standard"
    echo "executable = ${coire_path}"
    echo "arguments = --glob septime_files/${injstring}_${cat}/\$(macroinfile) --output first_coire_files/${injstring}/\$(macrooutfile) --data-type all_data --coinc-stat effective_snrsq --cluster-time 10000 --injection-file \$(macroinjectionfile) --injection-window 100"
    echo "log = " `mktemp -p ${log_path}`
    echo "error = logs/coire-\$(cluster)-\$(process).err"
    echo "output = logs/coire-\$(cluster)-\$(process).out"
    echo "notification = never"
    echo "priority = ${condor_priority}"
    echo "queue 1"
  fi > first_coire_${injstring}.coire.sub
done
echo -e "\n...done."

#get septime_slide files
/bin/echo -n "Generating septime_slide file list..."
num_septimes=0
pushd septime_files/${cat}/ > /dev/null
for file in *SEPTIME_SLIDE_H*xml.gz; do
  echo ${file}
  num_septimes=$(( ${num_septimes} + 1 ))
done > ../septime_slide_${cat}.cache

for file in *SEPTIME_SLIDE_L*xml.gz; do
  echo ${file}
  num_septimes=$(( ${num_septimes} + 1 ))
done >> ../septime_slide_${cat}.cache
popd > /dev/null
echo " done."

#generate time slide coire dag
echo "Generating first_coire_slide.dag and .sub files..."
septime_idx=1
if [ 1 ]; then
  for infile in `cat septime_files/septime_slide_${cat}.cache`; do
    echo -ne "processing ${septime_idx} / ${num_septimes}\r" >&2
    septime_idx=$(( ${septime_idx} + 1))
    outfile=`echo $infile | sed s/SEPTIME_SLIDE/COIRE_SLIDE_${cat}/g`
    echo "JOB $outfile first_coire_slide.coire.sub"
    echo "RETRY $outfile 1"
    echo "VARS $outfile macroinfile=\"$infile\" macrooutfile=\"$outfile\""
    echo "CATEGORY $outfile first_coire_slide"
    echo "## JOB $outfile requires input file $infile"
  done 
  echo "MAXJOBS coire_slide 40"
fi > first_coire_slide.dag

if [ 1 ]; then
  echo "universe = vanilla"
  echo "executable = ${coire_path}"
  echo "arguments = --glob septime_files/${cat}/\$(macroinfile) --output first_coire_files/${cat}/\$(macrooutfile) --data-type all_data --coinc-stat effective_snrsq --cluster-time 10000 --num-slides 50"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/coire_slide-\$(cluster)-\$(process).err"
  echo "output = logs/coire_slide-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > first_coire_slide.coire.sub
echo -e "\n...done."

#setup directory structure
if [ ! -d first_coire_files ]; then
  mkdir first_coire_files
fi

for string in ${cat} BBHLININJ BBHLOGINJ BBHSPINLININJ BBHSPINLOGINJ BNSLININJ BNSLOGINJ BNSSPINLININJ BNSSPINLOGINJ NSBHLININJ NSBHLOGINJ NSBHSPINLININJ NSBHSPINLOGINJ BNS001INJ; do
  if [ ! -d first_coire_files/${string} ]; then
    mkdir first_coire_files/${string}
  fi
done

echo "******************************************************"
echo "  Now run: condor_submit_dag first_coire.dag"
echo "      and: condor_submit_dag first_coire_slide.dag"
echo "      and: condor_submit_dag first_coire_injection.dag"
echo "  These dags can be run simutaneously."
echo "******************************************************"

