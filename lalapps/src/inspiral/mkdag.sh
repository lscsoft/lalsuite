#!/bin/bash

function error() {
	echo "${0}: error: ${*}" >&2
	exit 1
}

unset GPS_START GPS_END TRIG_START TRIG_END DAGFILE
while getopts "a:b:e:f:s:t:" OPT ; do
	case "${OPT}" in
	# GPS start time
	s)
		GPS_START="${OPTARG}" ;;
	
	# GPS end time
	e)
		GPS_END="${OPTARG}" ;;
	
	# GPS trigger start time
	a)
		TRIG_START="${OPTARG}" ;;
	
	# GPS trigger end time
	b)
		TRIG_END="${OPTARG}" ;;
	
	# DAG filename
	f)
		DAGFILE="${OPTARG}" ;;
	
	# OnASys run-time configuration top-dir
	t)
		DATADIR="${OPTARG}" ;;
	esac
done

if [ $((${GPS_END} - ${GPS_START})) -lt 2048 ] ; then
cat > ${DAGFILE} << EOF
job a true.sub
EOF
cat > true.sub << EOF
universe = vanilla
executable = /bin/true
log = true.log
queue
EOF
else
cp ${DATADIR}/online.ini .
echo 0 ${GPS_START} ${GPS_END} $((${GPS_END} - ${GPS_START})) > segment.txt
lalapps_inspiral_online_pipe --config-file online.ini --log-path /usr1/dbrown/E12/L1 || exit 1
ln -sf online.dag ${DAGFILE}
fi
