#!/bin/sh

echo "Check job outputs and search for errors..."

here=`pwd`

for d in injections* training* noninjections; do

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

    # check the number of files consistency
    n_add=`ls ligolw_add*.out | wc | awk '{print $1}'`
    n_burca=`ls ligolw_burca*.out | wc | awk '{print $1}'`
    n_sqlite=`ls ligolw_sqlite*.out | wc | awk '{print $1}'`

    if [ ${n_add} -ne ${n_burca} ]; then 
	echo "*****************************************************"
	echo "***** ${d}/logs/ : missing ligolw_add or ligolw_burca jobs"
    fi

    if [ ${n_add} -ne ${n_sqlite} ]; then 
	echo "*****************************************************"
	echo "***** ${d}/logs/ : missing ligolw_add or ligolw_sqlite jobs"
    fi



    cd ${here}

done

echo "Error scan is finished."
exit
