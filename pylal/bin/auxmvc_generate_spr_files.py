#!/usr/bin/python

# Copyright (C) 2011 Rahul Biswas, Ruslan Vaulin, Kari Hodge
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

import os
import sys 
import string
from optparse import *
import glob
from commands import *
import subprocess
from pylal import auxmvc_utils
from pylal import git_version
import copy
import numpy
import random
usage = """
        This code generates training/evaluation files for MVSC classifier.
        """

def RoundRobin(glitch_list, clean_list, number):
    primary_glitch_set=glitch_list[number]
    primary_clean_set=clean_list[number]
    secondary_glitch_set_indices=[i for i in range(len(glitch_list)) if i != number]
    secondary_clean_set_indices=[i for i in range(len(clean_list)) if i != number]

     
    secondary_clean_set=clean_list[secondary_clean_set_indices[0]]
    for i in secondary_glitch_set_indices[1:]:
       secondary_clean_set=numpy.concatenate((secondary_clean_set,clean_list[i]))
    
    secondary_glitch_set=glitch_list[secondary_glitch_set_indices[0]]
    for i in secondary_glitch_set_indices[1:]:
       secondary_glitch_set=numpy.concatenate((secondary_glitch_set,glitch_list[i]))


    return  primary_clean_set, primary_glitch_set, secondary_clean_set, secondary_glitch_set

###########


def GenerateKWAuxGlitchTriggers(files):

   """Reads in the kw1-35_glitches_training sets and stores them in the memory 
   """
   KWAuxGlitchTriggers=auxmvc_utils.ReadKWAuxTriggers(files)
   return KWAuxGlitchTriggers


###########


def GenerateKWAuxCleanTriggers(files):

   """Reads in the kw1-35_signal_training sets and stores them in the memory 
   """
   KWAuxRandomTriggers = auxmvc_utils.ReadKWAuxTriggers(files)
   KWAuxCleanTriggers = auxmvc_utils.get_clean_samples(KWAuxRandomTriggers)
   return KWAuxCleanTriggers


def ApplyTimeWindow(AuxTriggers, time_window=0.1):
  """
  Apply time window keeping for each trigger only those auxiliary channel samples that statisfy dt <= time_window.
  """
  
  variables = list(AuxTriggers.dtype.names)
  channels = [var.split("_dt")[0] for var in variables if "_dt" in var]
  
  for row in AuxTriggers:
    for channel in channels:
      if abs(row[channel+"_dt"]) > 0.1:
        row[channel+"_signif"] = 0.0
        row[channel+"_dt"] = 0.0
        row[channel+"_dur"] = 0.0
        row[channel+"_freq"] = 0.0
        row[channel+"_npts"] = 0.0
        
  return AuxTriggers
        

##########

def parse_command_line():

  """
  Parser function dedicated
  """
  parser = OptionParser(version=git_version.verbose_msg) 
  parser.add_option("-c","--clean-paramsfile", help="file with events of class zero")
  parser.add_option("-g","--glitch-paramsfile", help="file with events of class one")
  parser.add_option( "", "--apply-time-window", action = "store_true", default = False, help ="Apply time window (typically narrower than used in creating samples)")
  parser.add_option("","--time-window",default=0.1,type="float",help="Time window to applied on samples (in seconds).") 
  parser.add_option("","--channels", help="file with the list of channels to be used in the analysis. If not given, all channels in the input data are used.")
  parser.add_option("","--unsafe-channels", help="file with the list of unsafe channels.")
  parser.add_option("-n","--roundrobin-number",default=10,type="int",help="number of round-robin training/testing sets to make")
  parser.add_option("","--dq-cats",action="store", type="string",default="ALL", help="Generate DQ veto categories" )
  parser.add_option("","--exclude-variables",action="store", type="string", default=None, help="Comma separated list of variables that should be excluded from MVSC parameter list" )
  parser.add_option("","--max-clean-samples",default=None,type="int",help="Maximum number of clean samples that wil be used in training")
  parser.add_option("","--max-glitch-samples",default=None,type="int",help="Maximum number of glitch samples that wil be used in training")  
  parser.add_option("","--output-tag",action="store",type="string", default=None, metavar=" OUTPUTTAG",\
      help="The output files will be named according to OUTPUTTAG" )

  (options,args) = parser.parse_args()

  return options, sys.argv[1:]


opts, args = parse_command_line()

###########

if not opts.clean_paramsfile:
  print >>sys.stderr, \
      "Must specify a clean triggers paramater text file"
  print >>sys.stderr, "Enter 'generate_mvsc_files.py --help' for usage"
  sys.exit(1)

if not opts.glitch_paramsfile:
  print >>sys.stderr, \
      "Must specify a glitch triggers paramater text file"
  print >>sys.stderr, "Enter 'generate_mvsc_files.py --help' for usage"
  sys.exit(1)



if opts.clean_paramsfile or opts.glitch_paramsfile is True:
 
  clean_paramsFile=[opts.clean_paramsfile]
  glitch_paramsFile=[opts.glitch_paramsfile]

    


  KWAuxCleanTriggers=GenerateKWAuxCleanTriggers(clean_paramsFile)
  print "Read in clean samples"
  
  # apply time window if required
  if opts.apply_time_window:
    KWAuxCleanTriggers = ApplyTimeWindow(KWAuxCleanTriggers, time_window=opts.time_window)
 
  if opts.channels:
    #construct list of channels to be excluded
    exclude_channels = []
    all_channels = []
    for var in KWAuxCleanTriggers.dtype.names:
      if "_sig" in var:
        all_channels.append(var.split("_sig")[0])
    include_channels = open(opts.channels).readlines()
    include_channels = [ch.strip() for ch in include_channels]
 
    for channel in all_channels:
      if not channel in include_channels:
        exclude_channels.append(channel)

  else:
    exclude_channels = []
   
  if opts.unsafe_channels:
    # read list of unsafe channels from file
    unsafe_file = open(opts.unsafe_channels,"r")
    for line in unsafe_file:
      exclude_channels.append(line.split("\n")[0])
    unsafe_file.close()

   
 
  KWAuxCleanTriggers= auxmvc_utils.FilterKWAuxTriggers(KWAuxCleanTriggers, opts.exclude_variables, exclude_channels)

  KWAuxCleanTriggers = auxmvc_utils.ShuffleKWAuxTriggers(KWAuxCleanTriggers)
  print "Clean triggers has been shuffled"



  # read in glitch samples
  KWAuxGlitchTriggers = GenerateKWAuxGlitchTriggers(glitch_paramsFile)
  
  # apply time window if required
  if opts.apply_time_window:
    KWAuxGlitchTriggers = ApplyTimeWindow(KWAuxGlitchTriggers, time_window=opts.time_window)  

  KWAuxGlitchTriggers = auxmvc_utils.FilterKWAuxTriggers(KWAuxGlitchTriggers, opts.exclude_variables, exclude_channels) 
  print "Read in glitch samples"
  
  KWAuxGlitchTriggers = auxmvc_utils.ShuffleKWAuxTriggers(KWAuxGlitchTriggers)
  print "Glitch triggers has been shuffled"

     
  dq_cats=opts.dq_cats.split(",")
  #if opts.exclude_variables :
  #  exclude_variables_list = opts.exclude_variables.split(",")
  #else:
  #  exclude_variables_list = None

  if opts.roundrobin_number:
    List_of_Clean_KW_Sets_cats = auxmvc_utils.split_array(KWAuxCleanTriggers, Nparts =int(opts.roundrobin_number))
    print "Splited clean samples into " + str(opts.roundrobin_number) + " parts."
  for cat in dq_cats:
 
    KW_Glitch_Triggers_cats=auxmvc_utils.getKWAuxTriggerFromDQCAT(KWAuxGlitchTriggers, cat)  
    print "Got glitches from " + cat + " CBC DQ category"  
    if opts.roundrobin_number:

      List_of_Glitch_KW_Sets_cats = auxmvc_utils.split_array(KW_Glitch_Triggers_cats, Nparts = int(opts.roundrobin_number))
      print "Splited glitch samples into " + str(opts.roundrobin_number) + " parts."
       
      for i in range(len(List_of_Glitch_KW_Sets_cats)):
    
        Primary_Clean_set_cats, Primary_Glitch_set_cats, Secondary_Clean_set_cats, Secondary_Glitch_set_cats=RoundRobin(List_of_Glitch_KW_Sets_cats, List_of_Clean_KW_Sets_cats,i)
        
        # set the limit on the size of glitch training sample if necessary
        if opts.max_glitch_samples:          
          if len(Secondary_Glitch_set_cats) > opts.max_glitch_samples:
            Secondary_Glitch_set_cats = Secondary_Glitch_set_cats[random.sample(numpy.arange(len(Secondary_Glitch_set_cats)), opts.max_glitch_samples)]   

        # set the limit on the size of clean training set if necessary
        if opts.max_clean_samples:
          if len(Secondary_Clean_set_cats) > opts.max_clean_samples:
            Secondary_Clean_set_cats = Secondary_Clean_set_cats[random.sample(numpy.arange(len(Secondary_Clean_set_cats)), opts.max_clean_samples)]
        print "Constructed " + str(i) +" round robin set"
        MVSC_evaluation_set_cats=auxmvc_utils.ConvertKWAuxToMVSC(KWAuxGlitchTriggers = Primary_Glitch_set_cats, KWAuxCleanTriggers = Primary_Clean_set_cats, ExcludeVariables = ["DQ2","DQ3","DQ4"])
        print "Converted evaluation KW to MVSC triggers"
        MVSC_training_set_cats=auxmvc_utils.ConvertKWAuxToMVSC(KWAuxGlitchTriggers = Secondary_Glitch_set_cats, KWAuxCleanTriggers = Secondary_Clean_set_cats, ExcludeVariables = ["DQ2","DQ3", "DQ4"])
        print "Converted training KW to MVSC triggers"
   
        output_evaluation=cat + "_" + opts.output_tag + "_set_" + str(i) + "_" + "evaluation.pat"
        auxmvc_utils.WriteMVSCTriggers(MVSC_evaluation_set_cats, output_filename = output_evaluation, Classified = False) 
                
        print "Wrote " + output_evaluation + " file to disk"

        output_training=cat + "_" + opts.output_tag + "_set_" + str(i) + "_" + "training.pat"
        auxmvc_utils.WriteMVSCTriggers(MVSC_training_set_cats, output_filename = output_training, Classified = False)
        print "Wrote " + output_training + " file to disk"










