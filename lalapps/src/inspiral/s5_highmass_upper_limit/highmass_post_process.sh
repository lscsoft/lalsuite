#Executables (put the paths to tagged versions when appropriate)
SQLITE=/usr/bin/sqlite3
LIGOLW_SQLITE=/archive/home/channa/opt/bin/ligolw_sqlite
LIGOLW_INSPINJFIND=./ligolw_inspinjfind
LALAPPS_NEWCORSE=/usr1/channa/pylal/bin/lalapps_newcorse
LIGOLW_SEGMENTS=/archive/home/channa/opt/bin/ligolw_segments

#Necessary input files from ihope
INJCACHE=/archive/home/channa/analysis/final_highmass_analysis_injections/871147814-875232014/highmass_ihope.cache
FULLDATACACHE=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/highmass_ihope.cache
H1VETOSEGMENTS=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/segments/H1-CATEGORY_3_VETO_SEGS-871147814-4084200.txt
H2VETOSEGMENTS=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/segments/H2-CATEGORY_3_VETO_SEGS-871147814-4084200.txt
L1VETOSEGMENTS=/archive/home/channa/analysis/final_highmass_analysis/871147814-875232014/segments/L1-CATEGORY_3_VETO_SEGS-871147814-4084200.txt



#Make veto XML files
echo "Using " ${LIGOLW_SEGMENTS} " to make vetoes_CAT_3.xml.gz..."
${LIGOLW_SEGMENTS} --insert-from-segwizard=H1=${H1VETOSEGMENTS} --insert-from-segwizard=H2=${H2VETOSEGMENTS} --insert-from-segwizard=L1=${L1VETOSEGMENTS} --name vetoes -o vetoes_CAT_3.xml.gz
echo "Done."

#Run thinca_to_coinc 

#Prepare and cluster the sqlite files
for f in *.sqlite
	do
	echo "Simplifying " ${f}
	${SQLITE} ${f} < simplify.sql
        echo "Removing H1H2 in " ${f}
	${SQLITE} ${f} < remove_h1h2.sql
        echo "Clustering " ${f}
	${SQLITE} ${f} < cluster.sql 
	done

#Put the injections back to XML so that the injection
#finder will work
for f in *INJ*.sqlite
	do
	${LIGOLW_SQLITE} -d ${f} -v -x ${f}.xml.gz
	done

#Find the injections
${LIGOLW_INSPINJFIND} -v *INJ*.sqlite.xml.gz

#Put the injections back to DB 
#WARNING REPLACES ORIGINAL DB FOR INJECTIONS
for f in *INJ*.sqlite
        do
        ${LIGOLW_SQLITE} -d ${f} -v -r ${f}.xml.gz
        done

#Run new corse *Actual total mass range should be ~ [25-100], 
#here boundaries are 0 and inf just to catch things that aren't 
#quite in the boundaries
${LALAPPS_NEWCORSE} -i H1,H2,L1 -b 0,50,85,inf -p thinca --veto-segments vetoes_CAT_3.xml.gz --veto-segments-name=vetoes -v *.sqlite

#Run plot ifar

#Run upper limit


