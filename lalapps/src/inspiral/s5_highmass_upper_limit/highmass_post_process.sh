#Executables (put the paths to tagged versions when appropriate)
SQLITE=/usr/bin/sqlite3
LIGOLW_SQLITE=/archive/home/channa/opt/bin/ligolw_sqlite
LIGOLW_INSPINJFIND=./ligolw_inspinjfind
LALAPPS_NEWCORSE=/usr1/channa/pylal/bin/lalapps_newcorse
LIGOLW_SEGMENTS=/archive/home/channa/opt/bin/ligolw_segments
LIGOLW_THINCA_TO_COINC=/archive/home/channa/opt/bin/ligolw_thinca_to_coinc


#Necessary input files from ihope
INJCACHE=/archive/home/channa/analysis/final_highmass_analysis_injections/871147814-875232014/highmass_ihope.cache
FULLDATACACHE=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/highmass_ihope.cache
H1VETOSEGMENTS=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/segments/H1-CATEGORY_3_VETO_SEGS-871147814-4084200.txt
H2VETOSEGMENTS=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/segments/H2-CATEGORY_3_VETO_SEGS-871147814-4084200.txt
L1VETOSEGMENTS=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/segments/L1-CATEGORY_3_VETO_SEGS-871147814-4084200.txt

CATS="CAT_2 CAT_3"

#for cat in $CATS
#       do
#       echo "Using " ${LIGOLW_SEGMENTS} " to make vetoes"
#       if ${LIGOLW_SEGMENTS} --insert-from-segwizard=H1=${H1VETOSEGMENTS} --insert-from-segwizard=H2=${H2VETOSEGMENTS} --insert-from-segwizard=L1=${L1VETOSEGMENTS} --name vetoes -o vetoes_${cat}.xml.gz; then echo "Done."; else exit; fi
#       done


#Run thinca_to_coinc on zero lag and time slides
for type in 'FULL_DATA'
        do
        for cat in $CATS
                do
                grep ${type}'.*'${cat} ${FULLDATACACHE} > ${type}${cat}.cache
                echo "Using " ${LIGOLW_THINCA_TO_COINC} " on " ${type} ${cat} " ..."
                condor_run process_full.sh ${LIGOLW_THINCA_TO_COINC} ${type} ${LIGOLW_SQLITE} ${SQLITE} ${cat} &
                done
        done

#Run thinca_to_coinc on injections
for f in $(grep HL- ${INJCACHE} | awk '{print $2}')
        do
        INJ_DESC=$(echo ${f} | sed -e 's@_@ @g' | awk '{print $3"_"$4"_"$5}')
        INJ_FILE=$(grep HL-${f} ${INJCACHE} |awk '{print $5}')
        for cat in $CATS
                do
                grep ${INJ_DESC}'.*'${cat} ${INJCACHE} > ${INJ_DESC}${cat}.cache
                echo "Running process_inj.sh on " ${INJ_DESC} ${cat} " ..."
                condor_run process_inj.sh ${LIGOLW_THINCA_TO_COINC} ${INJ_DESC} ${LIGOLW_SQLITE} ${INJ_FILE} ${SQLITE} ${LIGOLW_INSPINJFIND} ${cat} &
                done
        done

wait

exit

#Stuff here down still experimental

#Clean up after ourselves
find -name "S5_HM*.xml.gz" -delete
find -name "S5_HM_INJ*.xml.gz" -delete

#Run new corse *Actual total mass range should be ~ [25-100], 
#here boundaries are 0 and inf just to catch things that aren't 
#quite in the boundaries
for cat in $CATS
	do
	echo "Running new corse on " ${cat} 
	condor_run process_ifar.sh ${LALAPPS_NEWCORSE} ${cat} &
	done

wait


#Run upper limit


