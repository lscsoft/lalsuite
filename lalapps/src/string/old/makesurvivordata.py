#!/usr/bin/python

import sys
import os
import re
import time
from types import *
from optparse import OptionParser

from matplotlib.patches import Patch
from pylab import *
from pylal import readMeta
from pylal import viz

# grab command line options
parser = OptionParser(version="%prog CVS $Id$")
parser.add_option("-s", "--gps-start-time", metavar="SECONDS", default=None, help="start from GPS time SECONDS")
parser.add_option("-d", "--duration", metavar="SECONDS", default=None, help="duration")
options, args = parser.parse_args()
del parser

if options.gps_start_time:
	options.gps_start_time = int(options.gps_start_time)
if options.duration:
	options.duration = int(options.duration)

# interpret the remaining command line arguments as trigger file names
triggers = readMeta.metaDataTable(args, "sngl_burst")

# extract some arrays from the trigger table
peak_time = triggers.mkarray("peak_time")

i=1

os.system('rm -rf plotdata/')
os.system('mkdir plotdata/')

for time in peak_time:
	print 'Producing plot data for trigger at: ', time
	seg_file=open('seg.txt',mode='r')
	seg_list=seg_file.readlines()
        #do the slides
	for seg in seg_list:
		[crap,Tstart,Tend,duration]=seg.split(None,4)
		if time > int(Tstart) and time < int(Tend):
			
			H1cache_filename='cache/H-H1_RDS_C01_LX-'+Tstart+'-'+Tend+'.cache'
			os.system('~/lscsoft/lalapps/src/string/lalapps_StringSearch --bw-flow 27 --sample-rate 4096 \
			--bank-lowest-hifreq-cutoff 75 --bank-freq-start 50 --threshold 4 --frame-cache '+H1cache_filename+ \
			' --channel-name H1:LSC-STRAIN --gps-start-time '+str(int(time)-7)+' --gps-end-time '+str(int(time)+8) \
			+' --settling-time 0.5 --pad 4 --short-segment-duration 2 --cusp-search --print-data > plotdata/H1data'+str(i)+'.txt')

			H2cache_filename='cache/H-H2_RDS_C01_LX-'+Tstart+'-'+Tend+'.cache'
			os.system('~/lscsoft/lalapps/src/string/lalapps_StringSearch --bw-flow 27 --sample-rate 4096 \
			--bank-lowest-hifreq-cutoff 75 --bank-freq-start 50 --threshold 4 --frame-cache '+H2cache_filename+ \
			' --channel-name H2:LSC-STRAIN --gps-start-time '+str(int(time)-7)+' --gps-end-time '+str(int(time)+8) \
			+' --settling-time 0.5 --pad 4 --short-segment-duration 2 --cusp-search --print-data > plotdata/H2data'+str(i)+'.txt')
			
			L1cache_filename='cache/L-L1_RDS_C01_LX-'+Tstart+'-'+Tend+'.cache'
			os.system('~/lscsoft/lalapps/src/string/lalapps_StringSearch --bw-flow 27 --sample-rate 4096 \
			--bank-lowest-hifreq-cutoff 75 --bank-freq-start 50 --threshold 4 --frame-cache '+L1cache_filename+ \
			' --channel-name L1:LSC-STRAIN --gps-start-time '+str(int(time)-7)+' --gps-end-time '+str(int(time)+8) \
			+' --settling-time 0.5 --pad 4 --short-segment-duration 2 --cusp-search --print-data > plotdata/L1data'+str(i)+'.txt')

			i=i+1
