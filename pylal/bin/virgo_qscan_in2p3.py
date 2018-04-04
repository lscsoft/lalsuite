#!/usr/bin/env python
"""
Manager program to run specific parts of the followup on the CCIN2P3 cluster
"""

__author__ = 'Damir Buskulic <buskulic@lapp.in2p3.fr>'

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, copy, math, random
import socket, time
import re, string
from optparse import *
from time import strftime
import tempfile
import ConfigParser
import urlparse
import urllib
from UserDict import UserDict
#sys.path.append('/archive/home/buskulic/opt/s5_1yr_followup_20080131/lalapps/lib/python2.4/site-packages/lalapps')
from pylal import git_version

##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################

######################## OPTION PARSING  #####################################
usage = """
usage: %prog [options]
"""
parser = OptionParser(usage, version=git_version.verbose_msg)

parser.add_option("-o", "--output-directory",action="store",type="string",\
    help="output result directory")
parser.add_option("-T", "--times-file-scan",action="store",type="string",\
    help="file containing a list of times")
parser.add_option("-t", "--times-file-scanlite",action="store",type="string",\
    help="file containing a list of times")
parser.add_option("-S", "--configuration-file-scan",action="store",type="string",\
    help="configuration file for the scan process")
parser.add_option("-s", "--configuration-file-scan-seismic",action="store",type="string",\
    help="configuration file for the seismic scan process")
parser.add_option("-L", "--configuration-file-scanlite",action="store",type="string",\
    help="configuration file for the scanlite process")
parser.add_option("-l", "--configuration-file-scanlite-seismic",action="store",type="string",\
    help="configuration file for the seismic scanlite process")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#########  READING TIMES FILE AND LAUNCHING BATCH QSCANS SCRIPTS  ############

depIfoDir = './'
if not opts.times_file_scan:
   if os.path.exists(depIfoDir+'/TIMES/qscan_times.txt'):
      depQscanTFile =  open(depIfoDir+'/TIMES/qscan_times.txt','r')
else:
   depQscanTFile =  open(depIfoDir+'/'+opts.times_file_scan,'r')

if not opts.times_file_scanlite:
   if os.path.exists(depIfoDir+'/TIMES/background_qscan_times.txt'):
      depQscanLiteTFile =  open(depIfoDir+'/TIMES/background_qscan_times.txt','r')
else:
   depQscanLiteTFile =  open(depIfoDir+'/'+opts.times_file_scanlite,'r')

if not opts.output_directory:
   outputDir = 'RESULTS'
else:
   outputDir = opts.output_directory

if not opts.configuration_file_scan:
   depConfigScan = 'CONFIG/fg-rds-qscan_config.txt'
else:
   depConfigScan = opts.configuration_file_scan

if not opts.configuration_file_scan_seismic:
   depConfigSeismicScan = 'CONFIG/fg-seismic-qscan_config.txt'
else:
   depConfigSeismicScan = opts.configuration_file_scan_seismic

if not opts.configuration_file_scanlite:
   depConfigScanLite = 'CONFIG/bg-rds-qscan_config.txt'
else:
   depConfigScanLite = opts.configuration_file_scanlite

if not opts.configuration_file_scanlite_seismic:
   depConfigSeismicScanLite = 'CONFIG/bg-seismic-qscan_config.txt'
else:
   depConfigSeismicScanLite = opts.configuration_file_scanlite_seismic

# Save the commands in a file so that one can relaunch them in case of problems

timeCommands = strftime('%Y-%m-%d_%H:%M:%S')
outCommands = open('commands-'+timeCommands+'.out','w')
outCommands.write('List of commands submitted for qscans')
outCommands.write(' at time :'+timeCommands+'\n\n')

if os.path.exists(depIfoDir+'/TIMES/qscan_times.txt'):
   print '***'
   print '*** Foreground qscans'
   print '***'
   qscanLines = depQscanTFile.readlines()
   for qscanTimeRaw in qscanLines:
      qscanTime = qscanTimeRaw.rstrip('\n')
      print 'Launching foreground qscan for time '+qscanTime
      qscanCommand = './SCRIPTS/qsub_wscan.sh '+qscanTime+' '+depConfigScan+' '+outputDir+'/results_fg-rds-qscan/'+str(float(qscanTime))+' @foreground@'
      print '      command : '+qscanCommand
      outCommands.write('command for time '+qscanTime+'\n')
      outCommands.write(qscanCommand+'\n\n')
      os.system(qscanCommand)
      qscanCommand = './SCRIPTS/qsub_wscan.sh '+qscanTime+' '+depConfigSeismicScan+' '+outputDir+'/results_fg-seismic-qscan/'+str(float(qscanTime))+' @foreground-seismic@'
      print '      command : '+qscanCommand
      outCommands.write('command for time '+qscanTime+'\n')
      outCommands.write(qscanCommand+'\n\n')
      os.system(qscanCommand)

if os.path.exists(depIfoDir+'/TIMES/background_qscan_times.txt'):
   print '***'
   print '*** Background qscans (qscanlite)'
   print '***'
   qscanLines = depQscanLiteTFile.readlines()
   for qscanTimeRaw in qscanLines:
      qscanTime = qscanTimeRaw.rstrip('\n')
      print 'Launching background qscan (qscanlite) for time '+qscanTime
      qscanCommand = './SCRIPTS/qsub_wscanlite.sh '+qscanTime+' '+depConfigScanLite+' '+outputDir+'/results_bg-rds-qscan/'+str(float(qscanTime))
      print '      command : '+qscanCommand
      outCommands.write('command for time '+qscanTime+'\n')
      outCommands.write(qscanCommand+'\n\n')
      os.system(qscanCommand)
      qscanCommand = './SCRIPTS/qsub_wscanlite.sh '+qscanTime+' '+depConfigSeismicScanLite+' '+outputDir+'/results_bg-seismic-qscan/'+str(float(qscanTime))
      print '      command : '+qscanCommand
      outCommands.write('command for time '+qscanTime+'\n')
      outCommands.write(qscanCommand+'\n\n')
      os.system(qscanCommand)

outCommands.close()

##########################################################################
sys.exit(0)
