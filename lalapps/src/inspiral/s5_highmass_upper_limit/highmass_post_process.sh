#Executables (put the paths to tagged versions when appropriate)
SQLITE=/usr/bin/sqlite3
LIGOLW_SQLITE=/archive/home/channa/opt/bin/ligolw_sqlite
LIGOLW_INSPINJFIND=./ligolw_inspinjfind
LALAPPS_NEWCORSE=/usr1/channa/pylal/bin/lalapps_newcorse
LIGOLW_SEGMENTS=/archive/home/channa/opt/bin/ligolw_segments
LIGOLW_THINCA_TO_COINC=./ligolw_thinca_to_coinc

#Necessary input files from ihope
INJCACHE=/archive/home/channa/analysis/final_highmass_analysis_injections/871147814-875232014/highmass_ihope.cache
FULLDATACACHE=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/highmass_ihope.cache
H1VETOSEGMENTS=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/segments/H1-CATEGORY_3_VETO_SEGS-871147814-4084200.txt
H2VETOSEGMENTS=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/segments/H2-CATEGORY_3_VETO_SEGS-871147814-4084200.txt
L1VETOSEGMENTS=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/segments/L1-CATEGORY_3_VETO_SEGS-871147814-4084200.txt



#Make veto XML files
echo "Using " ${LIGOLW_SEGMENTS} " to make vetoes_CAT_3.xml.gz..."
if ${LIGOLW_SEGMENTS} --insert-from-segwizard=H1=${H1VETOSEGMENTS} --insert-from-segwizard=H2=${H2VETOSEGMENTS} --insert-from-segwizard=L1=${L1VETOSEGMENTS} --name vetoes -o vetoes_CAT_3.xml.gz; then echo "Done."; else exit; fi

#Run thinca_to_coinc on zero lag and time slides
echo "Using " ${LIGOLW_THINCA_TO_COINC} " on full data"
if ${LIGOLW_THINCA_TO_COINC} --cache-description=FULL_DATA --ihope-cache=${FULLDATACACHE} --veto-segments=vetoes_CAT_3.xml.gz --veto-segments-name=vetoes --output=FULL_DATA.xml.gz --verbose; then echo "Done."; else exit; fi

echo "Using " ${LIGOLW_SQLITE} " to insert into database"
if ${LIGOLW_SQLITE} -d FULL_DATA.sqlite -t /tmp -v FULL_DATA.xml.gz; then echo "Done."; else exit; fi

#clean up after ourselves
echo "cleaning up xml files..."
if rm *FULL_DATA.xml.gz; then echo "Done."; else exit; fi

#Run thinca_to_coinc on injections
for f in $(grep HL- ${INJCACHE} | awk '{print $2}')
	do 
	INJ_DESC=$(echo ${f} | sed -e 's@_@ @g' | awk '{print $3"_"$4"_"$5}')
	INJ_FILE=$(grep HL-${f} ${INJCACHE} |awk '{print $5}')

	echo "Using " ${LIGOLW_THINCA_TO_COINC} " on " ${INJ_DESC} " ..."
	if ${LIGOLW_THINCA_TO_COINC} --cache-description=${INJ_DESC} --ihope-cache=${INJCACHE} --veto-segments=vetoes_CAT_3.xml.gz --veto-segments-name=vetoes --output=${INJ_DESC}.xml.gz --verbose; then echo "Done."; else exit; fi

        echo "Using " ${{LIGOLW_SQLITE} " on " ${INJ_DESC} ".xml.gz ..."
	if ${LIGOLW_SQLITE} -d ${INJ_DESC}.sqlite -t /tmp -v ${INJ_DESC}.xml.gz; then echo "Done."; else exit; fi
	#clean up after ourselves

        echo "Using " ${{LIGOLW_SQLITE} " to insert sims from " ${INJ_FILE} " ..."
        if ${LIGOLW_SQLITE} -d ${INJ_DESC}.sqlite -t /tmp -v ${INJ_FILE}; then echo "Done."; else exit; fi
        #clean up after ourselves

	echo "cleaning up xml files..."
	if rm *${INJ_DESC}.xml.gz; then echo "Done."; else exit; fi
	done

#Prepare and cluster the sqlite files
for f in *.sqlite
	do
	echo "Simplifying " ${f}
	if ${SQLITE} ${f} < simplify.sql; then echo "Done."; else exit; fi
        echo "Removing H1H2 in " ${f}
	if ${SQLITE} ${f} < remove_h1h2.sql; then echo "Done."; else exit; fi
        echo "Clustering " ${f}
	if ${SQLITE} ${f} < cluster.sql; then echo "Done."; else exit; fi
	done

#Put the injections back to XML so that the injection
#finder will work
for f in *INJ*.sqlite
	do
	echo "transforming " ${f} " to XML..."
	if ${LIGOLW_SQLITE} -d ${f} -v -x ${f}.xml.gz; then echo "Done."; else exit; fi
	done

#Find the injections
echo "Finding injections..."
if ${LIGOLW_INSPINJFIND} -v *INJ*.sqlite.xml.gz; then echo "Done."; else exit; fi

#Put the injections back to DB 
#WARNING REPLACES ORIGINAL DB FOR INJECTIONS
for f in *INJ*.sqlite
        do
        echo "transforming " ${f} " to sqlite..."
        if ${LIGOLW_SQLITE} -d ${f} -v -r ${f}.xml.gz; then echo "Done."; else exit; fi
        done

#Run new corse *Actual total mass range should be ~ [25-100], 
#here boundaries are 0 and inf just to catch things that aren't 
#quite in the boundaries
echo "Running new corse..."
if ${LALAPPS_NEWCORSE} -i H1,H2,L1 -b 0,50,85,inf -p thinca --veto-segments vetoes_CAT_3.xml.gz --veto-segments-name=vetoes -v *.sqlite; then echo "Done."; else exit; fi

#Run plot ifar

#Run upper limit


