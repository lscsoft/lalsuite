#!/bin/bash 

################################################################################
month_gps_time=`cat write_ifar_scripts_lv.ini | grep 'month_gps_time' | awk '{print $3}'`
month_duration=`cat write_ifar_scripts_lv.ini | grep 'month_duration' | awk '{print $3}'`
cat=`cat write_ifar_scripts_lv.ini | grep 'cat' | awk '{print $3}'`

coire_path=`cat write_ifar_scripts_lv.ini | grep 'coire_path' | awk '{print $3}'`
corse_path=`cat write_ifar_scripts_lv.ini | grep 'corse_path' | awk '{print $3}'`

log_path=`cat write_ifar_scripts_lv.ini | grep 'log_path' | awk '{print $3}'`
condor_priority=`cat write_ifar_scripts_lv.ini | grep 'condor_priority' | awk '{print $3}'`

#Print options out to screen for verification
echo "Options used are:"
echo "  month_gps_time = ${month_gps_time}"
echo "  month_duration = ${month_duration}"
echo "  cat = ${cat}"
echo "  coire_path = ${coire_path}"
echo "  corse_path = ${corse_path}"
echo "  log_path = ${log_path}"
echo "  condor_priority = ${condor_priority}"
echo
################################################################################

#get septime txt files organized by combo and run add_septime to add up times
echo "Getting SEPTIME_SLIDE*.txt file names and writing to cache files..."
pushd septime_files/${cat}/ > /dev/null
for combo in H1H2L1V1 H1H2L1 H1H2V1 H1L1V1 H2L1V1 H1H2 H1L1 H2L1 H1V1 H2V1 L1V1; do
  echo -n "  Getting ${cat}/${combo}-SEPTIME_SLIDE*.txt files..."
  for file in ${combo}-SEPTIME_SLIDE_*.txt; do
    echo "septime_files/${cat}/${file}"
  done > ../septime_${cat}_${combo}_times.cache
  echo " done."
done
popd > /dev/null
echo "...done."

#generate corse dag
/bin/echo -n "Generating corse_all_data.dag and .sub files..."
if [ 1 ]; then
  #write add_septime job
  for combo in H1H2L1V1 H1H2L1 H1H2V1 H1L1V1 H2L1V1 H1L1 H2L1 H1V1 H2V1 L1V1; do
    infile=septime_${cat}_${combo}_times.cache
    outfile=${combo}_V3_${cat}.txt
    echo "JOB $outfile corse_all_data.add_septime.sub"
    echo "RETRY $outfile 0"
    echo "VARS $outfile macrooutfile=\"$outfile\" macroinfile=\"$infile\""
    echo "CATEGORY $outfile add_septime"
    echo
  done
  # setting up corse jobs for zero-lag 
  for data in all_data playground_only exclude_play; do
    for mass in mchirp_2_8 mchirp_8_17 mchirp_17_35; do
    #write corse jobs for double-in_double files
      for combo in H1L1 H2L1 H1V1 H2V1 L1V1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_${cat}_${combo}-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_${combo}-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/${combo}_${combo}-CORSE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/${combo}_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="${combo}_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done
    #write corse jobs for double-in_triple and triple-in_triple files
      for combo in H1L1 H2L1 H1H2L1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_${cat}_H1H2L1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2L1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/H1H2L1_${combo}-CORSE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1H2L1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="H1H2L1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done
	  
    # write corse jobs for H1L1V1 time
      for combo in H1L1V1 H1V1 H1L1 L1V1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_${cat}_H1L1V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1L1V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/H1L1V1_${combo}-CORSE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1L1V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="H1L1V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done
	
	# write corse jobs for H2L1V1 time
      for combo in H2L1V1 L1V1 H2V1 H2L1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_${cat}_H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/H2L1V1_${combo}-CORSE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H2L1V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="H2L1V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done

    # write corse jobs for H1H2V1 time
      for combo in H1H2V1 H1V1 H2V1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_${cat}_H1H2V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/H1H2V1_${combo}-CORSE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1H2V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="H1H2V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done

	# write corse jobs for H1H2L1V1 time
      for combo in H1H2L1V1 H1H2L1 H1L1V1 H2L1V1 H1H2V1 H1L1 H2L1 H1V1 H2V1 L1V1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_${cat}_H1H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/H1H2L1V1_${combo}-CORSE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1H2L1V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="H1H2L1V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done
    done
  done
  
  # set up corse jobs for time-slides
  
  for data in all_data playground_only exclude_play; do
    for mass in mchirp_2_8 mchirp_8_17 mchirp_17_35; do
    #write corse jobs for double-in_double files
      for combo in H1L1 H2L1 H1V1 H2V1 L1V1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_${combo}-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_${combo}-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/${combo}_${combo}-CORSE_SLIDE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/${combo}_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="${combo}_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done
    #write corse jobs for double-in_triple and triple-in_triple files
      for combo in H1L1 H2L1 H1H2L1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2L1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2L1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/H1H2L1_${combo}-CORSE_SLIDE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1H2L1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="H1H2L1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done
	  
    # write corse jobs for H1L1V1 time
      for combo in H1L1V1 H1V1 H1L1 L1V1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1L1V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1L1V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/H1L1V1_${combo}-CORSE_SLIDE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1L1V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="H1L1V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done
	
	# write corse jobs for H2L1V1 time
      for combo in H2L1V1 L1V1 H2V1 H2L1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/H2L1V1_${combo}-CORSE_SLIDE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H2L1V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="H2L1V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done

    # write corse jobs for H1H2V1 time
      for combo in H1H2V1 H1V1 H2V1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/H1H2V1_${combo}-CORSE_SLIDE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1H2V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="H1H2V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done

	# write corse jobs for H1H2L1V1 time
      for combo in H1H2L1V1 H1H2L1 H1L1V1 H2L1V1 H1H2V1 H1L1 H2L1 H1V1 H2V1 L1V1; do
        zerofile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${cat}/${data}/H1H2L1V1_${combo}-CORSE_SLIDE_`echo ${data} | tr '[a-z]' '[A-Z]'`_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1H2L1V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"$data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\""
        echo "CATEGORY $outfile corse"
        parent_file="H1H2L1V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done
    done
  done

  ########### SET UP INJECTIONS #######

  for injstring in BBHLININJ BBHLOGINJ BBHSPINLININJ BBHSPINLOGINJ BNSLININJ BNSLOGINJ BNSSPINLININJ BNSSPINLOGINJ NSBHLININJ NSBHLOGINJ NSBHSPINLININJ NSBHSPINLOGINJ BNS001INJ; do

    #get HL-INJ file
    for file in `ls ../*inj/HL*${injstring}*`; do
      hlinjfile=$file
    done

    # setting up corse jobs for injections 
    for mass in mchirp_2_8 mchirp_8_17 mchirp_17_35; do
    #write corse jobs for double-in_double files
      for combo in H1L1 H2L1 H1V1 H2V1 L1V1; do
        zerofile="second_coire_files/${injstring}/${mass}/${combo}-SECOND_COIRE_${injstring}_${cat}_${combo}-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_${combo}-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${injstring}/${combo}_${combo}-CORSE_FOUND_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        missedfile="corse_all_data_files/${injstring}/${combo}_${combo}-CORSE_MISSED_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/${combo}_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.injection.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"all_data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\" macroinjectionfile=\"$hlinjfile\" macromissedinjectionfile=\"$missedfile\""
        echo "CATEGORY $outfile corse"
        parent_file="${combo}_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done

      for combo in H1H2L1 H1L1 H2L1; do
        zerofile="second_coire_files/${injstring}/${mass}/${combo}-SECOND_COIRE_${injstring}_${cat}_H1H2L1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2L1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${injstring}/H1H2L1_${combo}-CORSE_FOUND_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        missedfile="corse_all_data_files/${injstring}/H1H2L1_${combo}-CORSE_MISSED_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1H2L1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.injection.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"all_data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\" macroinjectionfile=\"$hlinjfile\" macromissedinjectionfile=\"$missedfile\""
        echo "CATEGORY $outfile corse"
        parent_file="H1H2L1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done

      for combo in H1H2V1 H1V1 H2V1; do
        zerofile="second_coire_files/${injstring}/${mass}/${combo}-SECOND_COIRE_${injstring}_${cat}_H1H2V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${injstring}/H1H2V1_${combo}-CORSE_FOUND_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        missedfile="corse_all_data_files/${injstring}/H1H2V1_${combo}-CORSE_MISSED_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1H2V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.injection.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"all_data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\" macroinjectionfile=\"$hlinjfile\" macromissedinjectionfile=\"$missedfile\""
        echo "CATEGORY $outfile corse"
        parent_file="H1H2V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done

      for combo in H1L1V1 H1L1 H1V1 L1V1; do
        zerofile="second_coire_files/${injstring}/${mass}/${combo}-SECOND_COIRE_${injstring}_${cat}_H1L1V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1L1V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${injstring}/H1L1V1_${combo}-CORSE_FOUND_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        missedfile="corse_all_data_files/${injstring}/H1L1V1_${combo}-CORSE_MISSED_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1L1V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.injection.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"all_data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\" macroinjectionfile=\"$hlinjfile\" macromissedinjectionfile=\"$missedfile\""
        echo "CATEGORY $outfile corse"
        parent_file="H1L1V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done

      for combo in H2L1V1 H2L1 H2V1 L1V1; do
        zerofile="second_coire_files/${injstring}/${mass}/${combo}-SECOND_COIRE_${injstring}_${cat}_H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${injstring}/H2L1V1_${combo}-CORSE_FOUND_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        missedfile="corse_all_data_files/${injstring}/H2L1V1_${combo}-CORSE_MISSED_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H2L1V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.injection.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"all_data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\" macroinjectionfile=\"$hlinjfile\" macromissedinjectionfile=\"$missedfile\""
        echo "CATEGORY $outfile corse"
        parent_file="H2L1V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done

      for combo in H1H2L1V1 H1H2L1 H1H2V1 H1L1V1 H2L1V1 H1L1 H1V1 H2L1 H2V1 L1V1; do
        zerofile="second_coire_files/${injstring}/${mass}/${combo}-SECOND_COIRE_${injstring}_${cat}_H1H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        slidefile="second_coire_files/${cat}/${mass}/${combo}-SECOND_COIRE_SLIDE_${cat}_H1H2L1V1-${month_gps_time}-${month_duration}.xml.gz"
        outfile="corse_all_data_files/${injstring}/H1H2L1V1_${combo}-CORSE_FOUND_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        missedfile="corse_all_data_files/${injstring}/H1H2L1V1_${combo}-CORSE_MISSED_${injstring}_${mass}_${cat}-${month_gps_time}-${month_duration}.xml.gz"
        summfile=`echo ${outfile} | sed s/.xml.gz/.txt/g`
        timeanalyzedfile="septime_files/H1H2L1V1_V3_${cat}.txt"
        echo "JOB $outfile corse_all_data.injection.corse.sub"
        echo "RETRY $outfile 1"
        echo "VARS $outfile macrozerofile=\"$zerofile\" macroslidefile=\"$slidefile\" macrooutfile=\"$outfile\" macrotimeanalyzedfile=\"$timeanalyzedfile\" macrodatatype=\"all_data\" macrosummaryfile=\"$summfile\" macromasstag=\"$mass\" macroinjectionfile=\"$hlinjfile\" macromissedinjectionfile=\"$missedfile\""
        echo "CATEGORY $outfile corse"
        parent_file="H1H2L1V1_V3_${cat}.txt"
        echo "PARENT $parent_file CHILD $outfile"
        echo
      done

    done
  done
  
  echo "MAXJOBS corse 20"
  echo "MAXJOBS add_septime 3"
fi > corse_all_data.dag

#write add_septime.sub file
if [ 1 ] ; then
  echo "universe = local"
  echo "executable = add_septime.py"
  # add --set-zero-lag-time-to 31536000 below to normilze IFARs to 1 year
  echo "arguments = --input-file septime_files/\$(macroinfile) --output-file septime_files/\$(macrooutfile) --num-slides 50 --set-zero-lag-time-to 31536000"
  echo "getenv = True"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/add_septime-\$(cluster)-\$(process).err"
  echo "output = logs/add_septime-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > corse_all_data.add_septime.sub
 
#write corse.sub files
if [ 1 ] ; then
  echo "universe = local"
  echo "executable = ${corse_path}"
  echo "arguments = --glob-zero \$(macrozerofile) --glob-slide \$(macroslidefile) --output \$(macrooutfile) --data-type \$(macrodatatype) --coinc-stat effective_snrsq --num-slides 50 --time-analyzed-file  \$(macrotimeanalyzedfile) --summary-file \$(macrosummaryfile) --mass-tag \$(macromasstag)"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/corse-\$(cluster)-\$(process).err"
  echo "output = logs/corse-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > corse_all_data.corse.sub
echo " done."

#write injection.corse.sub files
if [ 1 ] ; then
  echo "universe = local"
  echo "executable = ${corse_path}"
  echo "arguments = --glob-zero \$(macrozerofile) --glob-slide \$(macroslidefile) --output \$(macrooutfile) --data-type \$(macrodatatype) --coinc-stat effective_snrsq --num-slides 50 --time-analyzed-file  \$(macrotimeanalyzedfile) --summary-file \$(macrosummaryfile) --mass-tag \$(macromasstag) --injection-window 100 --missed-injections \$(macromissedinjectionfile) --injection-file \$(macroinjectionfile)"
  echo "log = " `mktemp -p ${log_path}`
  echo "error = logs/corse-\$(cluster)-\$(process).err"
  echo "output = logs/corse-\$(cluster)-\$(process).out"
  echo "notification = never"
  echo "priority = ${condor_priority}"
  echo "queue 1"
fi > corse_all_data.injection.corse.sub
echo " done."

#setup directory structure
/bin/echo -n "Creating corse_all_data_files directory structure..."
if [ ! -d corse_all_data_files ] ; then
  mkdir corse_all_data_files
fi
if [ ! -d corse_all_data_files/${cat} ] ; then
  mkdir corse_all_data_files/${cat}
fi
for data in all_data playground_only exclude_play; do
  if [ ! -d corse_all_data_files/${cat}/${data} ] ; then
    mkdir corse_all_data_files/${cat}/${data}
  fi
done

for injstring in BBHLININJ BBHLOGINJ BBHSPINLININJ BBHSPINLOGINJ BNSLININJ BNSLOGINJ BNSSPINLININJ BNSSPINLOGINJ NSBHLININJ NSBHLOGINJ NSBHSPINLININJ NSBHSPINLOGINJ BNS001INJ; do
  if [ ! -d corse_all_data_files/${injstring} ] ; then
    mkdir corse_all_data_files/${injstring}
  fi
done

echo " done."
echo "******************************************************"
echo "  Now run: condor_submit_dag corse_all_data.dag"
echo "******************************************************"

