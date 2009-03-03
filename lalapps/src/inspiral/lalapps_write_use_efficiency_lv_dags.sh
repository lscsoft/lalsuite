#!/bin/bash

################################################################################
# get needed options from ini file

month_gps_time=`cat write_ifar_scripts_lv.ini | grep 'month_gps_time' | awk '{print $3}'`
month_duration=`cat write_ifar_scripts_lv.ini | grep 'month_duration' | awk '{print $3}'`
cat=`cat write_ifar_scripts_lv.ini | grep 'cat' | awk '{print $3}'`

septime_path=`cat write_ifar_scripts_lv.ini | grep 'septime_path' | awk '{print $3}'`
findEfficiencyFactors_path=`cat write_ifar_scripts_lv.ini | grep 'findEfficiencyFactors_path' | awk '{print $3}'`
addEfficiencyFactors_path=`cat write_ifar_scripts_lv.ini | grep 'addEfficiencyFactors_path' | awk '{print $3}'`
findLoudestEvents_path=`cat write_ifar_scripts_lv.ini | grep 'findLoudestEvents_path' | awk '{print $3}'`
search_summary_path=`cat write_ifar_scripts_lv.ini | grep 'search_summary_path' | awk '{print $3}'`

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

# MAKE THE findEfficiencyFactors SUB FILE
################################################################################

if [ 1 ]; then
  echo "universe = local"
  echo "executable = ${findEfficiencyFactors_path}"
  echo "arguments = --input-cache-file ${hipe_cache} --found-injection-pattern \$(macrofoundinj) --missed-injection-pattern \$(macromissedinj) --veto-file vetoes_${cat}.xml.gz --h1-triggers --h2-triggers --l1-triggers --v1-triggers --statistic effective_snr --output-file \$(macrooutput) --verbose "
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/findEfficiencyFactors-\$(cluster)-\$(process).err"
  echo "output = logs/findEfficiencyFactors-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi >  findEfficiencyFactors.sub

# MAKE THE addEfficiencyFactors SUB FILE 
###############################################################################

if [ 1 ]; then
  echo "universe = local"
  echo "executable = ${addEfficiencyFactors_path}"
  echo "arguments = --config-file addEfficiencyFactors.ini --data-type all-data --output-file addefficiency.out"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/addEfficiencyFactors-\$(cluster)-\$(process).err"
  echo "output = logs/addEfficiencyFactors-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi >  addEfficiencyFactors.sub

# MAKE THE findLoudestEvents SUB FILE
##############################################################################

if [ 1 ]; then
  echo "universe = local"
  echo "executable = ${findLoudestEvents_path}"
  echo "arguments = --slides-glob=corse_all_data_files/CAT_3/all_data/IFO_TIME-CORSE_SLIDE_ALL_DATA_mchirp_bin*.xml.gz --h1-triggers --h2-triggers --l1-triggers --v1-triggers --statistic=effective_snr --num-slides=50 --verbose --output-file loudest_stats.stat"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/findLoudestEvents-\$(cluster)-\$(process).err"
  echo "output = logs/findLoudestEvents-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi >  findLoudestEvents.sub

# MAKE THE search_summary SUB FILE
##############################################################################

if [ 1 ]; then
  echo "universe = local"
  echo "executable = ${search_summary_path}"
  echo "arguments = --slides-glob "CAT_3/all_data/*SLIDE*" --zero-lag-glob \$(macrozerolag) --statistic lvS5stat --num-events 20 --save-background-stats --output-background-file \$(macropickle) --num-slides 50 --output-path \$(macrooutput) --enable-output"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/search_summary-\$(cluster)-\$(process).err"
  echo "output = logs/search_summary-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi >  search_summary.sub


# MAKE THE addEfficiencyFactors.ini FILE
###############################################################################

if [ 1 ]; then
  echo ";This is the configuration file for the likelihood DAG generation"
  echo ""
  echo "[playground-only]"
  echo "zero-lag-glob = corse_all_data_files/CAT_3/playground_only/*-CORSE_PLAYGROUND_ONLY_*.xml.gz"
  echo "slides-glob = corse_all_data_files/CAT_3/playground_only/*-CORSE_SLIDE_PLAYGROUND_ONLY_*.xml.gz"
  echo "no-background-candidates = no_background_candidates_playground_only.xml"
  echo ""
  echo "[all-data]"
  echo "zero-lag-glob = corse_all_data_files/CAT_3/all_data/*-CORSE_ALL_DATA*.xml.gz"
  echo "slides-glob = corse_all_data_files/CAT_3/all_data/*-CORSE_SLIDE_ALL_DATA*.xml.gz"
  echo "no-background-candidates = no_background_candidates_all_data.xml"
  echo ""
  echo "[injections]"
  echo "bns001inj  = "
  echo "nsbh001inj = "
  echo "bbh001inj  = "
  echo "spin001inj = "
  echo "fullrange001inj = "
  echo ""
  echo "[bns001inj]"
  echo "inj-glob = corse_all_data_files/BNS001INJ/*-CORSE_FOUND_BNS001INJ_*.xml.gz"
  echo ""
  echo "[nsbh001inj]"
  echo "inj-glob = corse_all_data_files/NSBH001INJ/*-CORSE_FOUND_NSBH001INJ_*.xml.gz"
  echo ""
  echo "[bbh001inj]"
  echo "inj-glob = corse_all_data_files/BBH001INJ/*-CORSE_FOUND_BBH001INJ_*.xml.gz"
  echo ""
  echo "[spin001inj]"
  echo "inj-glob = corse_all_data_files/SPIN001INJ/*-CORSE_FOUND_SPIN001INJ_*.xml.gz"
  echo ""
  echo "[fullrange001inj]"
  echo "inj-glob = corse_all_data_files/FULLRANGE001INJ/*-CORSE_FOUND_FULLRANGE001INJ_*.xml.gz"
  echo ""
  echo "[addEfficiencyFactors]"
  echo "sep-time-files-dir = septime_files"
  echo "EFactors-file = bns001inj_efficiency_factors.stat"
  echo "loudest-stats-file = loudest_stats.stat"
  echo "; choices are constant, exponent, two-sigma"
  echo "extrapolation-method = two-sigma"

fi > addEfficiencyFactors.ini


# MAKE THE DAG
###############################################################################

if [ 1 ]; then
  echo "JOB 1001 findEfficiencyFactors.sub"
  echo "VARS 1001 macrofoundinj=\"COIRE_INJECTIONS_1234_BNS001INJ_FOUND_SECOND_*_BNS001INJ_CAT_3_VETO\" macromissedinj=\"COIRE_INJECTIONS_1234_BNS001INJ_MISSED_SECOND_*_BNS001INJ_CAT_3_VETO\" macrooutput=\"bns001inj_efficiency_factors.stat\" "

  echo ""
  echo "JOB 1002 findLoudestEvents.sub"
 
  echo ""
  echo "JOB 1003 addEfficiencyFactors.sub"

  echo ""
  echo "JOB 1004 search_summary.sub"
  echo "VARS 1004 macrozerolag=\"CAT_3/playground_only/*LIKELIHOOD_PLAY*\" macropickle=\"playground_loudestbg.pickle\" macrooutput=\"search_summary_playground_only/\""

  echo ""
  echo "JOB 1005 search_summary.sub"
  echo "VARS 1005 macrozerolag=\"CAT_3/all_data/*LIKELIHOOD_ALL_DATA*\" macropickle=\"all_data_loudestbg.pickle\" macrooutput=\"search_summary_all_data/\""

  echo ""
  echo "PARENT 1001 CHILD 1002"
  echo "PARENT 1002 CHILD 1003"
  echo "PARENT 1003 CHILD 1004"
  echo "PARENT 1004 CHILD 1005"

fi > use_efficiency.dag
