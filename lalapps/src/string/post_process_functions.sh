#
# bash functions for use in the post-processing scripts.
#

function condor_run_burca() {
	# runs burca on the cluster, 1 job per argument

	# start a .sub file
	cat <<_EOF >ligolw_burca.sub
universe = vanilla
executable = \$ENV(HOME)/scratch/bin/ligolw_burca
arguments = --verbose --program StringSearch --window H1,H2=0.002 --window H1,L1=0.012 --window H2,L1=0.012 --string-compare 3.0,0.5 \$(macroinput)
getenv = True
log = logs/ligolw_burca.log
error = logs/ligolw_burca-\$(cluster)-\$(process).err
output = logs/ligolw_burca-\$(cluster)-\$(process).out
notification = never
_EOF

	# add a "queue" command for each argument
	for f in "$@" ; do
		echo -e "\nmacroinput = $f\nqueue 1" >>ligolw_burca.sub

	# submit the jobs
	condor_submit ligolw_burca.sub

	# wait until the word "termination" appears in the log as many
	# times as there were jobs
	while ! [ -f logs/ligolw_burca.log ] || [ $(grep termination logs/ligolw_burca.log | wc -l) != $# ] ; do
		sleep 10
	done

	# clean up
	rm -f ligolw_burca.sub
}
