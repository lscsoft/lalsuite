#!/usr/bin/env python
"""
Manager program to run specific parts of the followup on the CCIN2P3 cluster
"""

__author__ = 'Romain Gouaty  <gouaty@lapp.in2p3.fr>'

##############################################################################
# import standard modules

import sys
import tarfile
import os, shutil
from subprocess import *
import copy
import re
import exceptions
import glob
import fileinput
import linecache
import string
from optparse import *
from types import *
import operator
from UserDict import UserDict
sys.path.append('@PYTHONLIBDIR@')

##############################################################################
# import the modules we need from GLUE/LAL/LALAPPS/PYLAL

from pylal import git_version
from pylal import stfu_pipe

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


parser.add_option("-i","--qscan-input-file",action="store",type="string",\
    help="input file containing the qscan results. It must be a tar.gz file")

parser.add_option("","--qscan-cache-background",action="store",type="string", \
    help="qscan cache file for background")

parser.add_option("","--qscan-cache-foreground",action="store",type="string", \
    help="qscan cache file for foreground")

parser.add_option("","--remote-ifo",action="store",type="string",\
    help="name of the interferometer for which the qscans were performed \
    remotely (example: V1)")

parser.add_option("","--qscan-type-list",action="store",type="string",\
    help="types of qscan to be analysed (example: \"rds-qscan,seismic-qscan\")")

nd_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#################################
# Sanity check of input arguments

if not opts.qscan_input_file:
  print >> sys.stderr, "No input file specified!!"
  print >> sys.stderr, "Use --qscan-input-file FILE to specify location."
  sys.exit(1)

if not opts.qscan_cache_background and not opts.qscan_cache_foreground:
  print >> sys.stderr, "No cache file specified!!"
  print >> sys.stderr, "Use at least one of --qscan-cache-background FILE or --qscan-cache-foreground FILE to specify where the qscans results need to be copied."
  sys.exit(1)

if not opts.remote_ifo:
  print >> sys.stderr, "No ifo specified!!"
  print >> sys.stderr, "Use --remote-ifo to specify the interferometer (example \"V1\")."
  sys.exit(1)

if not opts.qscan_type_list:
  print >> sys.stderr, "No qscan type specified!!"
  print >> sys.stderr, "Use --qscan-type-list to specify the list of types (example \"rds-qscan,seismic-qscan\")."
  sys.exit(1)


#################################
# CORE OF THE PROGRAM

qscanShortTypeList = [qscanType.strip() for qscanType in opts.qscan_type_list.split(",")]
qscanTypeList = []

# check that at least one of the cache files is valid
if (opts.qscan_cache_background and not os.path.exists(opts.qscan_cache_background) and (opts.qscan_cache_foreground and not os.path.exists(opts.qscan_cache_foreground) ) ):
  print >> sys.stderr, "None of the qscan-cache-background or qscan-cache-foreground are valid!!"
  sys.exit(1)

# get the list of qscans in the cache files
if opts.qscan_cache_foreground and os.path.exists(opts.qscan_cache_foreground):
  for qscan_type in qscanShortTypeList:
    type_description = "WPIPELINE_FG_"
    if qscan_type == "rds-qscan":
      type_description += "RDS"
    elif qscan_type == "seismic-qscan":
      type_description += "SEIS"
    exec("fg_"+qscan_type.replace('-','_') + "List = stfu_pipe.getParamsFromCache(opts.qscan_cache_foreground,\""+type_description+"\",opts.remote_ifo)")
    qscanTypeList.append("fg-"+qscan_type)
else:
  if not opts.qscan_cache_foreground:
    print >> sys.stderr, "Foreground qscans won't be processed because option qscan-cache-foreground is not provided."
  else:
    print >> sys.stderr, "File " + opts.qscan_cache_foreground + " could not be found!!"

if opts.qscan_cache_background and os.path.exists(opts.qscan_cache_background):
  for qscan_type in qscanShortTypeList:
    type_description = "WPIPELINE_BG_" 
    if qscan_type == "rds-qscan":
      type_description += "RDS"
    elif qscan_type == "seismic-qscan":
      type_description += "SEIS"
    exec("bg_"+qscan_type.replace('-','_') + "List = stfu_pipe.getParamsFromCache(opts.qscan_cache_background,\""+type_description+"\",opts.remote_ifo)")
    qscanTypeList.append("bg-"+qscan_type)
else:
  if not opts.qscan_cache_background:
    print >> sys.stderr, "Background qscans won't be processed because option qscan-cache-background is not provided."
  else:
    print >> sys.stderr, "File " + opts.qscan_cache_background + " could not be found!!"


# Check for the existence of the qscan-input-file and try to uncompress it
if os.path.exists(opts.qscan_input_file):
  if tarfile.is_tarfile(opts.qscan_input_file):
    tar = tarfile.open(opts.qscan_input_file,"r:gz")
    file_list = tar.getnames()
    for file in file_list:
      tar.extract(file)
    tar.close()
  else:
    print >> sys.stderr, "File " + opts.qscan_input_file + " is not a valid tar file!!"
    sys.exit(1)
else:
  print >> sys.stderr, "File " + opts.qscan_input_file + " could not be found!!"
  sys.exit(1)

for qscan_type in qscanTypeList:
  # define the input path of the qscan directories to be copied
  result_path = opts.qscan_input_file.split(".tar.gz")[0]+ "/RESULTS/results_" + qscan_type + "/"
  # check whether the qscan input directory is empty or not. If it is empty we don't need to do anything
  if not os.listdir(result_path):
    print >> sys.stderr, "Directory "+result_path+" does not contain any result, it will be ignored"
    continue
  else:
    for qscan in eval(qscan_type.replace('-','_') + "List"):
      qscan_result_path = result_path + qscan[0].strip('/').split('/')[-1]
      # check if directory exists before trying to move it
      if os.path.exists(qscan_result_path):
        # do not overwrite an existing output directory
        if os.path.exists(qscan[0]):
          print >> sys.stderr, "Directory " + qscan[0] + " already exists, cannot be overwritten with new qscan results"
        else:
          shutil.move(qscan_result_path, qscan[0]) 
          #print "\n Copying file " + qscan_result_path + " to " + qscan[0]
      else:
        print >> sys.stderr, "Directory " + qscan_result_path + " could not be found!!"
        #sys.exit(1)


# Do some cleaning in the local directory before exiting
if os.path.exists(opts.remote_ifo + "_qscans_results"):
  shutil.rmtree(opts.remote_ifo + "_qscans_results")


