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
import sys
import logging
import optparse
import subprocess
import ConfigParser
from ligo.gracedb.rest import GraceDb
from laldetchar.idq import idq
from laldetchar.idq import idq_gdb_utils
import traceback

from laldetchar import git_version
from sys import stdin
# get_LVAdata_from_stdin will be depricate starting from ER6, we might need it only for pre-ER6 testing.
#from ligo.lvalert.utils import get_LVAdata_from_stdin
# json will be used instead of get_LVAdata_from_stdin. We will use it as default here.
import json



__author__ = \
    'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

description = \
    """ This program is launched by lvalert_listen upon receiving the alert about a new event submitted to GraceDB. It gets necessary information about the event and runs idq tasks that upload relevant data quality information to GraceDB."""

parser = optparse.OptionParser(version='Name: %%prog\n%s'
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
	
(options, args) = parser.parse_args()


# read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(options.config_file)

# get general settings from config file

ifo = config.get('general', 'ifo')
classifiers = config.get('general', 'classifiers').split(' ')
usertag = config.get('general', 'usertag')
idq_dir = config.get('general', 'idqdir')

# local gdb directory where script will run and store its output

idq_gdb_main_dir = config.get('general', 'main_gdb_dir')
if not os.path.exists(idq_gdb_main_dir):
    os.makedirs(idq_gdb_main_dir)

# cd into this directory
os.chdir(idq_gdb_main_dir)

# read gracedb id for event and other attributes from stdin

# get_LVAdata_from_stdin is used pre-ER6, see comments at the top.
#streamdata = get_LVAdata_from_stdin(stdin, as_dict=True)
eventdata = json.loads(stdin.read())

gdb_id = eventdata['uid']
#gdb_id = 'T121250'
alert_type = eventdata['alert_type']
#alert_type = 'new'
description = eventdata['description']

# to appear in ER6
#group = event['Group']
#pipeline = event['Pipeline']
#search = event['Search']

# The above event attributes could be useful in the future if we decide to code up logic
# that differentiates between gracedb entries based on their type, search, description etc.

# for now we treat all events in the same way

# logging setup

logger = logging.getLogger('idq_gdb_logger')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s %(message)s')
hdlr1 = logging.StreamHandler(sys.stdout)
hdlr1.setFormatter(formatter)
hdlr1.setLevel(logging.INFO)
logger.addHandler(hdlr1)
hdlr2 = logging.FileHandler(gdb_id+'_'+options.log_file)
hdlr2.setFormatter(formatter)
hdlr2.setLevel(logging.INFO)
logger.addHandler(hdlr2)

# redirect stdout and stderr into logger

sys.stdout = idq.LogFile(logger)
sys.err = idq.LogFile(logger)




# check that this is a new event. Exit if it is not ( not further processing is needed).
if not alert_type == 'new':
   sys.exit(0)

logger.info("New event. GraceDB ID %s." % (gdb_id))
logger.info("Alert type: %s. Description: %s." % (alert_type,description))

# make new directory for this event
#event_dir = gdb_id
#if not os.path.exists(event_dir):
#    os.makedirs(event_dir)

# We are using robot certificate, no need for finding and validating proxy certificate
# We only need to set environment to use the robot certificate.

#unset ligo-proxy just in case

del os.environ['X509_USER_PROXY']

# get cert and key from ini file

robot_cert = config.get('ldg_certificate', 'robot_certificate')
robot_key = config.get('ldg_certificate', 'robot_key')

# set cert and key

os.environ['X509_USER_CERT'] = robot_cert
os.environ['X509_USER_KEY'] = robot_key


# initialize instance of gracedb interface
gracedb = GraceDb()

# connect to gracedb and get event gps time
try: 
    gdb_entry = json.loads(gracedb.event(gdb_id).read())
except:
    traceback.print_exc()
    logger.info("    Error: Connection to GraceDB failed!")
    logger.info("    Exiting.")
    sys.exit(1)

event_gps_time = gdb_entry['gpstime']

# set the gps time range around the candidate for which idq data will be requested
# here we might want to code up some logic distinguishing event/search types e.g. Burst from CBC triggers
gps_start = event_gps_time - float(config.get('event', 'time_before'))
gps_end = event_gps_time + float(config.get('event', 'time_after'))


# run idq-gdb-timeseries for each classifier
for classifier in classifiers:
    logger.info("    Begin: executing idq-gdb-timeseries for " + classifier + " ...")
    exit_status = idq_gdb_utils.execute_gdb_timeseries(str(gps_start), str(gps_end), str(event_gps_time), gdb_id, ifo, classifier, config, idq_dir, config.get('executables', 'idq_gdb_timeseries'), usertag = str(gdb_id)+'_'+usertag)

    if exit_status != 0:
        logger.info("    WARNING: idq-gdb-timeseries failed for " + classifier)
        gracedb.writeLog(gdb_id, message="iDQ glitch-rank timeseries task failed for " + classifier + " at " + ifo)
    else:
        logger.info("    Done: idq-gdb-timeseries for " + classifier + ".")


	
# run idq-gdb-glitch-tables for each classifier
for classifier in classifiers:
    logger.info("    Begin: executing idq-gdb-glitch-tables for " + classifier + " ...")
    exit_status = idq_gdb_utils.execute_gdb_glitch_tables(str(gps_start), str(gps_end), gdb_id, ifo, classifier, config, idq_dir, config.get('executables', 'idq_gdb_glitch_tables'), usertag = str(gdb_id)+'_'+usertag)

    if exit_status != 0:
        logger.info("    WARNING: idq-gdb-glitch-tables failed for " + classifier)
        gracedb.writeLog(gdb_id, message="iDQ glitch tables task failed for " + classifier + " at " + ifo)
    else:
        logger.info("    Done: idq-gdb-glitch-tables for " + classifier + ".")





















