#!/usr/bin/env python

__prog__ = "wscan_bg_setup_log.py"
__title__ = "Set up a log file with the start time of the run"

##############################################################################

import os, sys, subprocess, shutil
from optparse import *
import ConfigParser
import time

from glue import gpstime
from pylal import git_version

##############################################################################

######################## OPTION PARSING  #####################################
usage = """usage: %prog [options]
"""
parser = OptionParser(usage, version=git_version.verbose_msg)

parser.add_option("-s","--start-run",action="store_true",default=False,\
	help="creates a log file and write the start time of the run inside it")

parser.add_option("-t","--terminate-run",action="store_true",default=False,\
        help="write the end time of the run inside the log file")

parser.add_option("-o","--output-path",action="store",type="string",\
	default="",help="path where the log file must be written")

parser.add_option("-n","--log-name",action="store",type="string",\
	default="",help="name of the log file to be recorded")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#############################################################################

if opts.start_run:
	if os.path.isdir(opts.output_path):
		if os.path.isfile(opts.output_path + "/latest-run.log"):
			os.remove(opts.output_path + "/latest-run.log")
		logfile = open(opts.output_path + "/latest-run.log","w")
		logfile.write("# omega run started at:\n")
		logfile.write(str(gpstime.GpsSecondsFromPyUTC(time.time())) + "\n")
		logfile.close()
		if opts.log_name:
			shutil.copyfile(opts.output_path + "/latest-run.log",opts.output_path + "/" + opts.log_name)
	else:
		print >> sys.stderr, "Path " + opts.output_path + " is not a valid directory" 
		sys.exit(0)

if opts.terminate_run:
	if os.path.isdir(opts.output_path):
		if not os.path.isfile(opts.output_path + "/latest-run.log"):
			print >> sys.stderr, "File " + opts.output_path + "/latest-run.log was not found!"
			sys.exit(1)
		else:
			logfile = open(opts.output_path + "/latest-run.log","a")
			logfile.write("# ended at:\n")
			logfile.write(str(gpstime.GpsSecondsFromPyUTC(time.time())) + "\n")
			logfile.close()
			if opts.log_name:
				shutil.copyfile(opts.output_path + "/latest-run.log",opts.output_path + "/" + opts.log_name)
	else:
		print >> sys.stderr, "Path " + opts.output_path + " is not a valid directory"
		sys.exit(0)

