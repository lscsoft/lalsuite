#!/bin/sh

echo "Check job outputs and search for errors..."

here=`pwd`

for d in injections* noninjections; do

    if [ ! -d ${d} ]; then continue; fi
    
    cd ${d}/logs

    # check if error files are empty
    for errfile in datafind* lalapps_binj*.out lalapps_StringSearch*.err ligolw_add* ligolw_burca* ligolw_sqlite*; do

	if [ -s ${errfile} ]; then
	    echo "*****************************************************"
	    echo "***** ${d}/logs/${errfile} :"
	    cat ${errfile}
	fi
    done


    # check if the string job went to the end
    for outfile in lalapps_StringSearch*.out; do
    

	if ! grep -q "StringJob is done" ${outfile}; then
	    echo "*****************************************************"
	    echo "***** ${d}/logs/${outfile} :"
	    echo "      JOB aborted or not finished !"
	fi
    done

    cd ${here}

done

echo "Error scan is finished."
exit
