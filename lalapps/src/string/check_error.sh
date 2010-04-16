#!/bin/sh

echo "Check job outputs and search for errors..."

here=`pwd`

for d in injections* noninjections; do
    
    cd ${d}/logs

    for errfile in datafind*.err lalapps_binj*.out lalapps_StringSearch*.err burca*.err; do

	if [ ! -s ${errfile} ]; then
	    echo "*****************************************************"
	    echo "***** ${d}/logs/${errfile} :"
	    cat ${errfile}
	fi
    done

    cd ${here}

done

echo "Error scan is finished."
exit
