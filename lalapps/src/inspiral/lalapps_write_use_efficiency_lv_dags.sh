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
sumtimes_path=`cat write_ifar_scripts_lv.ini | grep 'sumtimes_path' | awk '{print $3}'`
calculatefar_path=`cat write_ifar_scripts_lv.ini | grep 'calculatefar_path' | awk '{print $3}'`

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
  echo "arguments = --slides-glob \$(macroslideglob) --zero-lag-glob \$(macrozerolag) --statistic lvS5stat --num-events 20 --save-background-stats --output-background-file \$(macropickle) --num-slides 50 --output-path \$(macrooutput) --enable-output --make_zoomed_histogram --no-bg-zero-lag-glob \$(macronobg) --verbose --ignore-IFO-times H1H2_H1H2,H1H2_H1H2L1,H1H2_H1H2V1,H1H2_H1H2L1V1,H2L1_H1H2L1,H2L1_H1H2L1V1,H2V1_H1H2V1,H2V1_H1H2L1V1,H2L1V1_H1H2L1V1"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/search_summary-\$(cluster)-\$(process).err"
  echo "output = logs/search_summary-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi >  search_summary.sub

# MAKE THE sumtimes SUB FILE
##############################################################################

if [ 1 ]; then
  echo "universe = local"
  echo "executable = ${sumtimes_path}"
  echo "arguments = --input-glob \$(macroinputglob) --output-file combinedtimes_exclude_play.txt"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/sumtimes-\$(cluster)-\$(process).err"
  echo "output = logs/sumtimes-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi >  sumtimes.sub

# MAKE THE CalculateFAR_zerolag SUB FILE
##############################################################################

if [ 1 ]; then
  echo "universe = local"
  echo "executable = ${calculatefar_path}"
  echo "arguments = --slides-glob \$(macroslideglob) --events-glob \$(macroeventglob) --statistic \$(macrostat) --save-background-stats --output-background-file \$(macrobgfile) --verbose --output-path \$(macrooutputpath) --combine-output \$(macrocombineoutput) --num-slides 50 --septime-files \$(macroseptimefiles)"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/calculatefar-\$(cluster)-\$(process).err"
  echo "output = logs/calculatefar-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi >  calculatefar_zerolag.sub

# MAKE THE CalculateFAR (slides)  SUB FILE
##############################################################################
if [ 1 ]; then  echo "universe = local"
  echo "executable = ${calculatefar_path}"
  echo "arguments = --events-glob \$(macroeventglob) --statistic \$(macrostat) --skip-timeslides --background-stats-file \$(macrobgfile) --
verbose --output-path \$(macrooutputpath) --combine-output \$(macrocombineoutput) --num-slides 50 --septime-files \$(macroseptimefiles)"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/calculatefar-\$(cluster)-\$(process).err"
  echo "output = logs/calculdatefar-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi >  calculatefar_slides.sub

# MAKE THE CalculateFAR (inj)  SUB FILE
##############################################################################

if [ 1 ]; then
  echo "universe = local"
  echo "executable = ${calculatefar_path}"
  echo "arguments = --events-glob \$(macroeventglob) --statistic \$(macrostat) --skip-timeslides --background-stats-file \$(macrobgfile) --verbose --output-path \$(macrooutputpath) --combine-output \$(macrocombineoutput) --num-slides 50 --septime-files \$(macroseptimefiles) --ensure-search-summary-table-uniqueness"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/calculatefar-\$(cluster)-\$(process).err"
  echo "output = logs/calculdatefar-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority=${condor_priority}"
  echo "queue 1"
fi >  calculatefar_inj.sub

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
  echo "[exclude-playground-data]"
  echo "zero-lag-glob = corse_all_data_files/CAT_3/exclude_play/*-CORSE_EXCLUDE*.xml.gz"
  echo "slides-glob = corse_all_data_files/CAT_3/all_data/*-CORSE_SLIDE_ALL_DATA*.xml.gz"
  echo "no-background-candidates = no_background_candidates_exclude_play.xml"
  echo ""
  echo "[injections]"
  echo "bbhlininj = "
  echo "bbhloginj = "
  echo "bbhspinlininj = "
  echo "bbhspinloginj = "
  echo "bnslininj = "
  echo "bnsloginj = "
  echo "bnsspinlininj = "
  echo "bnsspinloginj = "
  echo "nsbhlininj = "
  echo "nsbhloginj = "
  echo "nsbhspinlininj = "
  echo "nsbhspinloginj = "
  echo "bns001inj = "
  echo ""
  echo "[bbhlininj]"
  echo "inj-glob = corse_all_data_files/BBHLININJ/*-CORSE_FOUND_BBHLININJ_*.xml.gz"
  echo ""
  echo "[bbhloginj]"
  echo "inj-glob = corse_all_data_files/BBHLOGINJ/*-CORSE_FOUND_BBHLOGINJ_*.xml.gz"
  echo ""
  echo "[bbhspinlininj]"
  echo "inj-glob = corse_all_data_files/BBHSPINLININJ/*-CORSE_FOUND_BBHSPINLININJ_*.xml.gz"
  echo ""
  echo "[bbhspinloginj]"
  echo "inj-glob = corse_all_data_files/BBHSPINLOGINJ/*-CORSE_FOUND_BBHSPINLOGINJ_*.xml.gz"
  echo ""
  echo "[bnslininj]"
  echo "inj-glob = corse_all_data_files/BNSLININJ/*-CORSE_FOUND_BNSLININJ_*.xml.gz"
  echo ""
  echo "[bnsloginj]"
  echo "inj-glob = corse_all_data_files/BNSLOGINJ/*-CORSE_FOUND_BNSLOGINJ_*.xml.gz"
  echo ""
  echo "[bnsspinlininj]"
  echo "inj-glob = corse_all_data_files/BNSSPINLININJ/*-CORSE_FOUND_BNSSPINLININJ_*.xml.gz"
  echo ""
  echo "[bnsspinloginj]"
  echo "inj-glob = corse_all_data_files/BNSSPINLOGINJ/*-CORSE_FOUND_BNSSPINLOGINJ_*.xml.gz"
  echo ""
  echo "[nsbhlininj]"
  echo "inj-glob = corse_all_data_files/NSBHLININJ/*-CORSE_FOUND_NSBHLININJ_*.xml.gz"
  echo ""
  echo "[nsbhloginj]"
  echo "inj-glob = corse_all_data_files/NSBHLOGINJ/*-CORSE_FOUND_NSBHLOGINJ_*.xml.gz"
  echo ""
  echo "[nsbhspinlininj]"
  echo "inj-glob = corse_all_data_files/NSBHSPINLININJ/*-CORSE_FOUND_NSBHSPINLININJ_*.xml.gz"
  echo ""
  echo "[nsbhspinloginj]"
  echo "inj-glob = corse_all_data_files/NSBHSPINLOGINJ/*-CORSE_FOUND_NSBHSPINLOGINJ_*.xml.gz"
  echo ""
  echo "[bns001inj]"
  echo "inj-glob = corse_all_data_files/BNS001INJ/*-CORSE_FOUND_BNS001INJ_*.xml.gz"
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
  echo "VARS 1004 macroslideglob=\"CAT_3/playground_only/*SLIDE*\"  macrozerolag=\"CAT_3/playground_only/*LIKELIHOOD_PLAY*\" macropickle=\"playground_loudestbg.pickle\" macrooutput=\"search_summary_playground_only/\" macronobg=\"CAT_3/playground_only/no_background*\""

  echo ""
  echo "JOB 1005 search_summary.sub"
  echo "VARS 1005 macroslideglob=\"CAT_3/all_data/*SLIDE*\"  macrozerolag=\"CAT_3/all_data/*LIKELIHOOD_ALL_DATA*\" macropickle=\"all_data_loudestbg.pickle\" macrooutput=\"search_summary_all_data/\" macronobg=\"CAT_3/all_data/no_background*\""

  echo ""
  echo "JOB 1006 search_summary.sub"
  echo "VARS 1006 macroslideglob=\"CAT_3/exclude_play/*SLIDE*\"  macrozerolag=\"CAT_3/exclude_play/*LIKELIHOOD_EXCLUDE_PLAY*\" macropickle=\"exclude_play_loudestbg.pickle\" macrooutput=\"search_summary_exclude_play/\" macronobg=\"CAT_3/exclude_play/no_background*\""

  echo ""
  echo "JOB 1007 sumtimes.sub"
  echo "VARS 1007 macroinputglob=\"corse_all_data_files/CAT_3/exclude_play/*CORSE_EXCLUDE*txt\""

  echo ""
  echo "JOB 1008 calculatefar_zerolag.sub"
  echo "VARS 1008 macroslideglob=\"CAT_3/exclude_play/*SLIDE*\" macroeventglob=\"search_summary_exclude_play/loudest_zerolag.xml\" macrostat=\"lvS5stat\" macrobgfile=\"exclude_play_bg.stat\" macrooutputpath=\"uncombined_ifar_files/exclude_play/\" macrocombineoutput=\"combined_ifar_files/exclude_play/combinedifar_zerolag.xml\" macroseptimefiles=\"septime_files/*CAT*txt\" "

  echo ""
  echo "JOB 1009 calculatefar_slides.sub"
  echo "VARS 1009 macroeventglob=\"search_summary_exclude_play/loudest_timeslides.xml\" macrostat=\"lvS5stat\" macrobgfile=\"exclude_play_bg.stat\" macrooutputpath=\"uncombined_ifar_files/exclude_play/\" macrocombineoutput=\"combined_ifar_files/exclude_play/combinedifar_slides.xml\" macroseptimefiles=\"septime_files/*CAT*txt\" "


  injcounter=1010
  for injstring in bbhlininj bbhloginj bbhspinlininj bbhspinloginj bnslininj bnsloginj bnsspinlininj bnsspinloginj nsbhlininj nsbhloginj nsbhspinlininj nsbhspinloginj bns001inj; do
    echo ""
    echo "JOB ${injcounter} calculatefar_inj.sub"
    echo "VARS ${injcounter} macroeventglob=\"${injstring}/*LIKELIHOOD*FOUND*xml\" macrostat=\"lvS5stat\" macrobgfile=\"exclude_play_bg.stat\" macrooutputpath=\"uncombined_ifar_files/${injstring}/\" macrocombineoutput=\"combined_ifar_files/${injstring}/combinedifar_${injstring}.xml\" macroseptimefiles=\"septime_files/*CAT*txt\" "
    injcounter+=1
  done

  echo ""
  echo "PARENT 1001 CHILD 1002"
  echo "PARENT 1002 CHILD 1003"
  echo "PARENT 1003 CHILD 1004"
  echo "PARENT 1004 CHILD 1005"
  echo "PARENT 1005 CHILD 1006"
  echo "PARENT 1006 CHILD 1007"
  echo "PARENT 1007 CHILD 1008"
  echo "PARENT 1008 CHILD 1009"
  echo "PARENT 1009 CHILD 1010"
  echo "PARENT 1010 CHILD 1011"
  echo "PARENT 1011 CHILD 1012"
  echo "PARENT 1012 CHILD 1013"
  echo "PARENT 1013 CHILD 1014"
  echo "PARENT 1014 CHILD 1015"
  echo "PARENT 1015 CHILD 1016"
  echo "PARENT 1016 CHILD 1017"
  echo "PARENT 1017 CHILD 1018"
  echo "PARENT 1018 CHILD 1019"
  echo "PARENT 1019 CHILD 1020"
  echo "PARENT 1020 CHILD 1021"
  echo "PARENT 1021 CHILD 1022"
  echo "PARENT 1022 CHILD 1023"

fi > use_efficiency.dag
