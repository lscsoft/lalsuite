# Copyright (C) 2013 Lindy Blackburn, Reed Essick and Ruslan Vaulin
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import os
#import glob
import sys
from sys import stdin

import logging
import traceback
import ConfigParser

import subprocess

from ligo.gracedb.rest import GraceDb

from laldetchar.idq import idq
#from laldetchar.idq import reed as idq
from laldetchar.idq import idq_gdb_utils

import time

from laldetchar import git_version

# get_LVAdata_from_stdin will be depricate starting from ER6, we might need it only for pre-ER6 testing.
#from ligo.lvalert.utils import get_LVAdata_from_stdin
# json will be used instead of get_LVAdata_from_stdin. We will use it as default here.
import json


from optparse import OptionParser

#===================================================================================================

__author__ = \
    'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

description = \
    """ This program is launched by lvalert_listen upon receiving the alert about a new event submitted to GraceDB. \
It gets necessary information about the event and runs idq tasks that upload relevant data quality information to GraceDB."""

#===================================================================================================

parser = OptionParser(version='Name: %%prog\n%s'
                               % git_version.verbose_msg,
                               usage='%prog [options]',
                               description=description)
parser.add_option(
    '-c',
    '--config-file',
    dest='config_file',
    help='configuration file',
    metavar='FILE',
    default='idq_gdb.ini',
    )
	
parser.add_option(
    '-l',
    '--log-file',
    default='idq_gdb.log',
    type='string',
    help='log file')

parser.add_option(
    "", "--no-robot-cert",
    default=False,
    action="store_true")
	
(options, args) = parser.parse_args()

#===================================================================================================

### read global configuration file
config = ConfigParser.SafeConfigParser()
config.read(options.config_file)

### get general settings from config file
ifo = config.get('general', 'ifo')
classifiers = config.get('general', 'classifiers').split(' ')
usertag = config.get('general', 'usertag')
idq_dir = config.get('general', 'idqdir')

realtimedir = config.get('general', 'realtimedir')

if usertag:
    usertag = "_%s"%usertag

### get script names
tables_script = config.get('executables','idq_gdb_glitch_tables')
timeseries_script = config.get('executables','idq_gdb_timeseries')

#=================================================
### local gdb directory where script will run and store its output
idq_gdb_main_dir = config.get('general', 'main_gdb_dir')
if not os.path.exists(idq_gdb_main_dir):
    os.makedirs(idq_gdb_main_dir)

### cd into this directory
os.chdir(idq_gdb_main_dir)
cwd = os.getcwd()

#=================================================
### read gracedb id for event and other attributes from stdin

### get_LVAdata_from_stdin is used pre-ER6, see comments at the top.
#streamdata = get_LVAdata_from_stdin(stdin, as_dict=True)
eventdata = json.loads(stdin.read())

gdb_id = eventdata['uid']
alert_type = eventdata['alert_type']
description = eventdata['description']

usertag = "%s%s"%(gdb_id, usertag)

### out/err files
tables_out = "%s/glitch-tables_%s.out"%(idq_gdb_main_dir, usertag)
tables_err = "%s/glitch-tables_%s.err"%(idq_gdb_main_dir, usertag)

timeseries_out = "%s/timeseries_%s.out"%(idq_gdb_main_dir, usertag)
timeseries_err = "%s/timeseries_%s.err"%(idq_gdb_main_dir, usertag)

#=================================================
# LOGIC TO SKIP EVENTS BASED ON MESSAGE DATA
#=================================================
if not alert_type=='new': ### only process new events
    sys.exit(0)

#=================================================
### logging setup

logger = idq.setup_logger('idq_gdb_logger', "%s_%s"%(gdb_id, options.log_file), sys.stdout, format='%(asctime)s %(message)s')

sys.stdout = idq.LogFile(logger)
sys.stderr = idq.LogFile(logger)

#=================================================
### process the event in full
logger.info("New event. GraceDB ID %s." % (gdb_id))
logger.info("Alert type: %s. Description: %s." % (alert_type,description))

### set up robot cert if needed
if not options.no_robot_cert:
    ### unset ligo-proxy just in case
    if os.environ.has_key("X509_USER_PROXY"):
        del os.environ['X509_USER_PROXY']

    ### get cert and key from ini file
    robot_cert = config.get('ldg_certificate', 'robot_certificate')
    robot_key = config.get('ldg_certificate', 'robot_key')

    ### set cert and key
    os.environ['X509_USER_CERT'] = robot_cert
    os.environ['X509_USER_KEY'] = robot_key


### initialize instance of gracedb interface
if config.has_option("general","gdb_url"):
    gdb_url = config.get('general', 'gdb_url')
    gracedb = GraceDb( gdb_url )
else:
    gdb_url = None
    gracedb = GraceDb()

### connect to gracedb and get event gps time
try: 
    gdb_entry = json.loads(gracedb.event(gdb_id).read())
except:
    traceback.print_exc()
    logger.info("    Error: Connection to GraceDB failed!")
    logger.info("    Exiting.")
    sys.exit(1)

#========================
# get parameters about event type from gracedb
#========================
group = gdb_entry['group']
pipeline = gdb_entry['pipeline']
if gdb_entry.has_key('search'):
    search = gdb_entry['search']
    event_type = "%s_%s_%s"%(group, pipeline, search)
else:
    event_type = "%s_%s"%(group, pipeline)

### cast to lower case to match config sections
event_type = event_type.lower()

#========================
# digest event type and pull the correct params from config file
#========================
if event_type not in config.sections():
    logger.info("Warning: event type not found. Defaulting to \"default_event\" settings in %s"%options.config_file)
    event_type = 'default_event'

logger.info("event_type : %s"%event_type)

time_before = config.getfloat(event_type, 'time_before')
time_after = config.getfloat(event_type, 'time_after')

plotting_time_before = config.getfloat(event_type, 'plotting_time_before')
plotting_time_after = config.getfloat(event_type, 'plotting_time_after')

max_wait = config.getfloat(event_type,'max_wait')
delay = config.getfloat(event_type, 'delay') ### the amount of extra time we wait for jobs to finish

#========================
# get start and end time for our look-up routines
#========================
event_gps_time = float(gdb_entry['gpstime'])

gps_start = event_gps_time - time_before
gps_end = event_gps_time + time_after

plotting_gps_start = event_gps_time - plotting_time_before
plotting_gps_end = event_gps_time + plotting_time_after

logger.info("Started searching for iDQ information within [%.3f, %.3f] at %s"%(gps_start, gps_end, ifo))
gracedb.writeLog(gdb_id, message="Started searching for iDQ information within [%.3f, %.3f] at %s"%(gps_start, gps_end, ifo))

#=================================================
# LOGIC for waiting for idq data 
#=================================================
### figure out if we need to wait for time to pass
wait = gps_end - (idq.nowgps()+delay)
if wait > 0:
    logger.info("waiting %.2f seconds until we pass gps=%.3f"%(wait, gps_end))
    time.sleep(wait)

### now we need to parse the realtime log to figure out where the realtime job is
realtime_logname = config.get("general","realtime_log")
logger.info("parsing %s to extract idq-realtime state"%realtime_logname)

realtime_log = open(realtime_logname, "r") ### open realtime log for reading
realtime_log.seek(0, 2) ### go to end of file

### wait until realtime has passed gps_end+delay
past, dead, timed_out = idq.block_until(gps_end, realtime_log, max_wait=max_wait, timeout=2*max_wait)

if past:
    logger.info("found realtime stride starting after t=%.3f"%(gps_end))

elif timed_out:
    logger.info("WARNING: could not find a recent enough stride in %s after searching for %.1f seconds. Realtime process may be behind"%(realtime_logname, 2*max_wait))
    gracedb.writeLog(gdb_id, message="WARNING: iDQ data quality coverage was not complete and no new information was reported after waiting %.1f seconds. Data quality information at %s may not be complete."%(2*max_wait, ifo))

else: # dead
    logger.info("WARNING: no new iDQ information was reported to %s after waiting %.1f seconds. Realtime process may be dead."%(realtime_logname, max_wait))
    gracedb.writeLog(gdb_id, message="WARNING: iDQ data quality coverage was not complete and no new information was reported after waiting %.1f seconds. Data quality information at %s may not be complete."%(max_wait, ifo))

#=================================================
# launch the actual processes
#=================================================
for classifier in classifiers:

    ### run idq-gdb-glitch-tables for each classifier
    logger.info("    Begin: executing idq-gdb-glitch-tables for " + classifier + " ...")

#    tables_cmd = "%s -s %f -e %f -g %s --ifo %s -c %s -i %s -t %s"%(tables_script, gps_start, gps_end, gdb_id, ifo, classifier, idq_dir, usertag)
    tables_cmd = "%s -s %f -e %f -g %s --ifo %s -c %s -i %s -t %s"%(tables_script, gps_start, gps_end, gdb_id, ifo, classifier, realtimedir, usertag)
    if config.has_option("general","gdb_url"):
        tables_cmd += " --gdb-url %s"%config.get('general','gdb_url')
    for option, value in config.items("gdb-glitch-tables"):
        tables_cmd += " --%s %s"%(option, value)

    logger.info('    Submitting glitch-tables script with the following options')
    logger.info(tables_cmd)

    tables_out_file = open(tables_out, 'a')
    tables_err_file = open(tables_err, 'a')
    proc = subprocess.Popen(tables_cmd.split(), stdout=tables_out_file, stderr=tables_err_file, cwd=cwd)
    tables_out_file.close()
    tables_err_file.close()
    exit_status = proc.wait() # block!

#    exit_status = idq_gdb_utils.execute_gdb_glitch_tables(str(gps_start), str(gps_end), gdb_id, ifo, classifier, config, idq_dir, config.get('executables', 'idq_gdb_glitch_tables'), usertag = str(gdb_id)+'_'+usertag)

    if exit_status != 0:
        logger.info("    WARNING: idq-gdb-glitch-tables failed for " + classifier)
        gracedb.writeLog(gdb_id, message="FAILED: iDQ glitch tables for " + classifier + " at " + ifo)
        gch_xml = None
    else:
        logger.info("    Done: idq-gdb-glitch-tables for " + classifier + ".")
        gch_xml = idq.gdb_xml(idq_gdb_main_dir, classifier, ifo, "_%s"%usertag, gps_start, gps_end-gps_start) ### compute name

    ### run idq-gdb-timeseries for each classifier
    logger.info("    Begin: executing idq-gdb-timeseries for " + classifier + " ...")

#    timeseries_cmd = "%s -s %f -e %f --gps %f -g %s --ifo %s -c %s -i %s -t %s --plotting-gps-start %f --plotting-gps-end %f"%(timeseries_script, gps_start, gps_end, event_gps_time, gdb_id, ifo, classifier, idq_dir, usertag, plotting_gps_start, plotting_gps_end)
    timeseries_cmd = "%s -s %f -e %f --gps %f -g %s --ifo %s -c %s -i %s -t %s --plotting-gps-start %f --plotting-gps-end %f"%(timeseries_script, gps_start, gps_end, event_gps_time, gdb_id, ifo, classifier, realtimedir, usertag, plotting_gps_start, plotting_gps_end)

    if config.has_option("general","gdb_url"):
        timeseries_cmd += " --gdb-url %s"%config.get("general", "gdb_url")

#    if gch_xml: ### currently broken, uncomment when fixed!
#        timeseries_cmd += " --gch-xml %s"%gch_xml

#    if cln_xml: ### never built!
#        timeseries_cmd += " --cln-xml %s"%scln

    for (option,value) in config.items('gdb-time-series'):
        timeseries_cmd += " --%s %s"%(option, value)

    logger.info('    Submitting timeseries script with the following options')
    logger.info(timeseries_cmd)

    timeseries_out_file = open(timeseries_out, 'a')
    timeseries_err_file = open(timeseries_err, 'a')
    proc = subprocess.Popen(timeseries_cmd.split(), stdout=timeseries_out_file, stderr=timeseries_err_file, cwd=cwd)
    timeseries_out_file.close()
    timeseries_err_file.close()
    exit_status = proc.wait() # block !

#    exit_status = idq_gdb_utils.execute_gdb_timeseries(str(gps_start), str(gps_end), str(event_gps_time), gdb_id, ifo, classifier, config, idq_dir, config.get('executables', 'idq_gdb_timeseries'), usertag = str(gdb_id)+'_'+usertag, plotting_gps_start=plotting_gps_start, plotting_gps_end=plotting_gps_end) #, gch_xml=glob.glob("%s/*glitch*%s*.xml*"%(idq_gdb_main_dir, gdb_id))) ### un-comment when adding tables is not broken!

    if exit_status != 0:
        logger.info("    WARNING: idq-gdb-timeseries failed for " + classifier)
        gracedb.writeLog(gdb_id, message="FAILED: iDQ glitch-rank timeseries for " + classifier + " at " + ifo)
    else:
        logger.info("    Done: idq-gdb-timeseries for " + classifier + ".")

#=================================================
### print clean-up statements
logger.info("Finished searching for iDQ information within [%.3f, %.3f] at %s"%(gps_start, gps_end, ifo))
gracedb.writeLog(gdb_id, message="Finished searching for iDQ information within [%.3f, %.3f] at %s"%(gps_start, gps_end, ifo))

