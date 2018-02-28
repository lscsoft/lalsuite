#!/usr/bin/env python
#
#       coincToLALInference.py
#
#       Copyright 2012
#       Vivien Raymond <vivien.raymond@ligo.org>
#       Chris Pankow <chris.pankow@ligo.org>
#       John Veitch <john.veitch@ligo.org>
#
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

usage = """
Python script to create all necessary files to run lalinference_mcmc in the 
working directory. 

       cd <submit_directory>
       ./coincToLALInference.py [options] [--mcmc|--nest] coinc.xml

It will use ligo_data_find to generate cache files, and print an error message
if not enough contiguous data is available.

Tested only on the CIT cluster.
"""

import sys, os
import math
import commands

from glue.ligolw import utils
from glue.ligolw import lsctables

from optparse import OptionParser

VALID_SETS = {
	"S5": 
		{'H1': 'H1_RDS_C03_L2',
		'H2': 'H2_RDS_C03_L2', 
		'L': 'L1_RDS_C03_L2', 
		'V': 'HrecV2_16384'},
	"S6VSR2":  
		{'H': 'H1_LDAS_C02_L2', 
		'L': 'L1_LDAS_C02_L2', 
		'V': 'HrecV3'},
	"S6VSR3": 
		{'H': 'H1_LDAS_C02_L2', 
		'L': 'L1_LDAS_C02_L2', 
		'V': 'HrecV2'},
	"ER1": 
		{'H': 'H1H2_DMT_ERHOFT', 
		'L': 'L1_DMT_ERHOFT', 
		'V': 'ER1_hoft'}
}

ligo_data_find_cmd = 'ligo_data_find -o %(ifo)s -t %(type)s -s %(gps_start)d -e %(gps_end)d -u file -l'


opt = OptionParser()
opt.usage = usage
opt.add_option("-T", "--default-frame-set", action="store", default='ER1', help="Use a default set of frame types for the data, e.g. S5, S6VSR2, ER1, etc. Allowed values are %s. Default is 'ER1'." % VALID_SETS.keys())
opt.add_option("-t", "--frame-data-type", action="store", help="Associate a data type with an instrument. E.g. --frame-data-type H:H1H2_DMT_ERHOFT. Use a comma separated list to indicate more than one. Specifying this after '-T' will override that instrument in the default set of values.")
opt.add_option("-p", "--psd-segments", action="store", type=int, default=32, help="Number of PSD segments to use. Default is 32.")
opt.add_option("-a", "--approximant", action="store", default="TaylorF2", help="Approximant to use. Default is TaylorF2.")
opt.add_option("-A", "--addtl-lalinf-args", action="store", help="Additional args to be sent to the lalinference_mcmc command line.")
opt.add_option("-v", "--verbose", action="store_true", help="Be verbose.")
opt.add_option("-i", "--nest-ini-file", action="store", help="Initialization file for lalapps_nest_pipe.")
#opt.add_option("-g", "--gracedb-id", action="store", help="Use information from indicated id.")
opt.add_option("--mcmc", action="store_true", help="Use lalinference_mcmc.")
opt.add_option("--nest", action="store_true", help="Use lalapps_nest_pipe.")
opt.add_option("-s", "--sub-filename", action="store", help="Make subfile with this name.")
# TODO: Integrate graceDB info

opts, args = opt.parse_args()
verb = opts.verbose

if not (opts.mcmc or opts.nest):
	sys.exit("Indicate one of --mcmc or --nest.")
elif opts.mcmc and opts.nest:
	sys.exit("Indicate one of --mcmc or --nest.")

"""
if len(args) == 0 and opts.gracedb_id is None:
	sys.exit("No coinc tables to process and no GDB entry to examine, exiting.")

gdbid = ""
if opts.gracedb_id is not None:
	# FIXME: Either enhance or remove this based on how GDB interaction 
	# crystallizes
	raise NotImplementedError
	gdbid = opts.gracedb_id
	if gdbid[0] != "G":
		if verb: print "Normalizing GDB ID..."
		gdbid = "G" + gdbid
	print "Using graceDB ID %s, ignoring other arguments." % gdbid
	# TODO: Retrieve coinc table and set argus to that
"""

if opts.nest and not opts.nest_ini_file:
	sys.exit("--nest option requires an initilization file")

types = opts.default_frame_set
if not types in VALID_SETS.keys():
	sys.exit("Invalid set type. Valid set choices are %s" % VALID_SETS.keys())
else:
	if verb: print "Using frame data type " + types
	types = VALID_SETS[types]

if opts.frame_data_type is not None:
	for dt in opts.frame_data_type.split(","):
		inst, type = dt.split(":")
		types[inst] = type

if verb: print "Full frame types requested: " + str(types)

psdsegments = opts.psd_segments

coinc_events = []

for arg in args:
	xmldoc = utils.load_filename(arg)
	coinctable = lsctables.CoincInspiralTable.get_table(xmldoc)
	sngltable = lsctables.SnglInspiralTable.get_table(xmldoc)
	coincmap = lsctables.CoincMapTable.get_table(xmldoc)

	coinc_events += [event.coinc_event_id for event in coinctable]
	if verb: print "Found %d coinc events in table." % len(coinc_events)

	for cevent in coinctable:
		cid = cevent.coinc_event_id
		end_time = cevent.end_time
		end_time_ns = cevent.end_time_ns
		if verb: print "Checking coinc event id %d, end time %10.2f" % (cid, end_time + end_time_ns)
		snr = cevent.snr
		# ifos = [event.ifos for event in insptable]

		sngl_insp_id = [map.event_id for map in coincmap if map.coinc_event_id == cid]
		sngl_insp = []
		for iid in sngl_insp_id:
			sngl_insp += [event for event in sngltable if event.event_id == iid]

		# Determine max template duration for use in constructing cache
		template_duration = max([event.template_duration for event in sngl_insp if event.event_id in sngl_insp_id])

		cachefilenames = []
		for event in sngl_insp:
			#template_duration = event.template_duration
			#channel += [event.channel for event in sngltable]
			if verb: print "Checking data for sngl_inspiral event %d" % event.event_id

			# Try to determine the frame type needed for this ifo
			dtype = None
			try: dtype = types[event.ifo]
			except KeyError: dtype = types[event.ifo[0]]

			e_info = {'ifo': event.ifo[0], 'type': dtype, 
				'gps_start': math.floor(end_time-template_duration*1.2*(psdsegments+1)) ,
				'gps_end': end_time
			}
			ligo_data_find_seg = (ligo_data_find_cmd % e_info ) + " --show-times"
			if verb: print ligo_data_find_seg
			status, output = commands.getstatusoutput(ligo_data_find_seg)
			# FIXME: Unhardcode
			if output.strip() == "No segments found!":
				sys.exit("Could not find any segments for ifo and type requested.")
			if output.count('\n') > 0:
				exit_str='Data for ' + event.ifo + ' is not contiguous. Aborting.\n'+output
				sys.exit(exit_str)
			if verb: print "Making cache for " + event.ifo
			status, output = commands.getstatusoutput(ligo_data_find_cmd % e_info)
			cachefilenames.append( '%s/%s-%d-%d.cache' % (os.getcwd(), event.ifo, end_time, cid) )
			cache_file = open( cachefilenames[-1],'w')
			cache_file.write(output)

		srate = [2**x for x in range(15)]
		freq=16384
		for f in reversed(srate):
			if f > max([event.f_final for event in sngl_insp]):
				freq = f

		if opts.mcmc:
			if verb: print "Constructing lalinference_mcmc command."
			command = '--trigtime ' + str(end_time) + '.' + str(end_time_ns)
			command += ' --trigSNR ' + str(snr)
			command += ' --seglen ' + str(math.floor(template_duration*1.2))
			command += ' --psdlength ' + str(math.floor(template_duration*1.2)*psdsegments)
			command += ' --psdstart ' + str(math.floor(end_time-template_duration*1.2*(psdsegments+1)))
			command += ' --ifo ' + "[%s]" % ",".join([e.ifo for e in sngl_insp])
			command += ' --cache ' + "[%s]" % ",".join(cachefilenames)
			command += ' --channel ' + "[%s]" % ",".join([e.ifo + ":" + e.channel for event in sngl_insp])
			# FIXME: Unhardcode... make option?
			flow = [ "10" for i in range(len(sngl_insp)) ]
			command += ' --flow ' + "[%s]" % ",".join(flow)
			command += ' --srate ' + str(freq)
			command += ' --approximant ' + opts.approximant
			command += ' --Dmax 1000 '
			#command += ' --Niter 5000000'
			if opts.addtl_lalinf_args is not None:
				command += opts.addtl_lalinf_args

			if verb: print command

			if verb: print "Creating submit file."

			if opts.sub_filename is None:
				submit_file = open('%s/lalinference_mcmc_%d.sub' % (os.getcwd(), cid),'w')
			else: 
				submit_file = open(opts.sub_filename,'w')

			submit_str ="""
universe = parallel
environment=CONDOR_MPI_PATH=/usr/lib64/openmpi
getenv = True
executable = /archive/home/vivien/condor_mpirun
arguments = --verbose --stdout cluster$(CLUSTER).proc$(PROCESS).mpiout --stderr cluster$(CLUSTER).proc$(PROCESS).mpierr /archive/home/vivien/master/bin/lalinference_mcmc -- %s
machine_count = 8
log = cluster$(CLUSTER).proc$(PROCESS).log
output = cluster$(CLUSTER).proc$(PROCESS).subproc$(NODE).out
error = cluster$(CLUSTER).proc$(PROCESS).subproc$(NODE).err
notification = Always
on_exit_remove = (ExitBySignal == True) || (ExitCode != 143)
rank = (40 - (2.0 * TotalCondorLoadAvg))
queue
""" % command
			submit_file.write(submit_str)

		elif opts.nest:
			# TODO: Make submit file for lalapps_nest_pipe
			command = 'lalapps_nest_pipe --coinc-triggers %s --run-path %s --ini-file %s --condor-submit --ignore-science-mode --coherence-test' % (arg, os.getcwd(), opts.nest_ini_file)
			retcode = os.execute(command)
			if retcode != 0:
				sys.exit("Failed to call lalapps_nest_pipe, check arguments.")

