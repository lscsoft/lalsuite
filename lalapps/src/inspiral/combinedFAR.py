#!/usr/bin/python

import sys
import os
from optparse import *
import re
import exceptions
import glob
from types import * #!!!what does this do?
from operator import itemgetter

from glue import lal
from glue import segments #!!!do I need this?
from glue import segmentsUtils #!!!do I need this?
from glue.ligolw import ligolw
from glue.ligolw import table as tab #!!!why as tab?
from glue.ligolw import lsctables
from glue.ligolw import utils
import glue.iterutils
from pylal import CoincInspiralUtils
from pylal import SnglInspiralUtils
from pylal import InspiralUtils
from pylal import SearchSummaryUtils

import matplotlib
matplotlib.use('Agg')
from pylab import *
from pylal import viz
from numpy import histogram

import itertools
import numpy
import operator

usage =  """
"""

def columnList():
  """
  A function to return a list of columns in the SimInspiralTable in the order
  expected.
  """
  column = []
  column.append('process_id')
  column.append('ifo')
  column.append('search')
  column.append('channel')
  column.append('end_time')
  column.append('end_time_ns')
  column.append('end_time_gmst')
  column.append('impulse_time')
  column.append('impulse_time_ns')
  column.append('template_duration')
  column.append('event_duration')
  column.append('amplitude')
  column.append('eff_distance')
  column.append('coa_phase')
  column.append('mass1')
  column.append('mass2')
  column.append('mchirp')
  column.append('mtotal')
  column.append('eta')
  column.append('tau0')
  column.append('tau2')
  column.append('tau3')
  column.append('tau4')
  column.append('tau5')
  column.append('ttotal')
  column.append('psi0')
  column.append('psi3')
  column.append('alpha')
  column.append('alpha1')
  column.append('alpha2')
  column.append('alpha3')
  column.append('alpha4')
  column.append('alpha5')
  column.append('alpha6')
  column.append('beta')
  column.append('f_final')
  column.append('snr')
  column.append('chisq')
  column.append('chisq_dof')
  column.append('sigmasq')
  column.append('rsqveto_duration')
  column.append('Gamma0')
  column.append('Gamma1')
  column.append('Gamma2')
  column.append('Gamma3')
  column.append('Gamma4')
  column.append('Gamma5')
  column.append('Gamma6')
  column.append('Gamma7')
  column.append('Gamma8')
  column.append('Gamma9')
  column.append('event_id')
  return column


def parse_command_line():
  """
  Parser function dedicated
  """

  parser = OptionParser( usage = usage, version = "%prog CVS $Id$ " )
  # following are related to file input and output naming
  parser.add_option( "-g", "--glob", action = "store", type = "string", \
        default = None, metavar = " GLOB", \
        help = "glob of CORSE files to read" )
  parser.add_option( "-R", "--min-rate", action = "store", type = "float", \
      default = None, metavar = "MINRATE", \
      help = "The minimum combined IFAR needed to output to xml for followup" )
  parser.add_option( "-I", "--cache-file", \
        help = "read CORSE file names from cache input file; " + \
        "Currently does not work.")
  parser.add_option( "-P", "--output-path", action = "store", type = "string", \
        default = "", metavar = "PATH", \
        help = "path where the figures should be stored" )
  parser.add_option( "-v", "--verbose", action = "store_true", default = False, \
        help = "print information to stdout" )
  parser.add_option( "", "--gps-start-time", action = "store", type = "int", \
        default = None, metavar = "GPSSTARTTIME", \
        help = "gps start time used in the figure and output file names" )
  parser.add_option( "", "--gps-end-time", action = "store", type = "int", \
        default = None, metavar = "GSPENDTIME", \
        help = "gps end time used in the figure and output file names" )
#  parser.add_option( "", "--num-slides", action = "store", type = "int", \
#        default = None, metavar = "NUMSLIDES", \
#        help = "The number of time slides that were performed." )
  parser.add_option( "", "--ifo-times", action = "store", default = None, metavar = "IFOS", \
        help = "puts the ifo times for which the plots were made into output file name." )
  parser.add_option( "", "--ifo-tag", action = "store", type = "string", \
        default = None, metavar = "IFOTAG", \
        help = "the ifo tag used in the name of the figure (e.g. SECOND_H1H2L1)" )
  parser.add_option( "-u", "--user-tag", action = "store", type = "string", \
        default = None, metavar = "USERTAG", \
        help = "user tag used in the name of the figures" )
  # following are ifar plot specific options
  parser.add_option( "", "--summary-file-path", action = "store",
        type = "string", default = None,
        help = "Directory to find corse summary-file. If none specified, " + \
        "plotifar will look in the same location as the relevant corse file." )
  parser.add_option( "", "--time-correct-file", metavar = "T_COR_FIL",\
        action = "store", type= "string", default = None,
        help="coire or corse  summary file that contains the amount of " + \
        "all_data analyzed time. " + \
        "Dividing the analysed time for the globbed files' data-type by " + \
        "this corrects the far if it was " + \
        "analyzed from an all_data septime file. ")


  (options,args) = parser.parse_args()
  #check if required options specified and for self-consistency
  if not options.glob or options.cache_file:
    raise ValueError, "--glob or --cache-file must be specified"
  if not options.output_path:
    raise ValueError, "--output-file must be specified"
  return options, sys.argv[1:]
  if not options.num_slides:
    raise ValueError, "--num-slides must be specified"


opts, args = parse_command_line()
opts.output_file = opts.output_path + '/' + opts.ifo_times.upper() + \
    '-CORSE_' + opts.user_tag + '_COMBINED_IFAR_CAT_3-'\
    + str(opts.gps_start_time) \
    + '-' + str(opts.gps_end_time) + '.xml.gz'
opts.output_file_loudest = opts.output_path + '/' + opts.ifo_times.upper() + \
    '-CORSE_' + opts.user_tag + '_COMBINED_IFAR_LOUDEST_CAT_3-'\
    + str(opts.gps_start_time) \
    + '-' + str(opts.gps_end_time) + '.xml.gz'
corsefiles = []
if opts.glob is not None:
  for gl in opts.glob.split(" "):
    corsefiles.extend(glob.glob(gl))
elif opts.cache_file is not None:
  # currently not a working feature; just print warning message and exit
  print >> sys.stderr, "--cache-file option currently not available; " + \
        "please  use --glob option instead."
  sys.exit(1)
if not corsefiles:
  print >> sys.stderr, "No corse files could be found. Check input args."
  sys.exit(1)

maxBkgFAN = {}
minBkgFAN = {}
modFAN = {}
massbin = {}
coincifos = {}
no_bkg = {}
no_frgnd = {}
no_bkg_frgnd = {}
coincT = {}
warn_msg = ""
for thisfile in corsefiles:
  if opts.summary_file_path: # use alternate path for summfile
    altfile = opts.summary_file_path + os.path.basename(thisfile)
    summfile = glob.glob(altfile.rstrip('.xml.gz') + '.txt')
  else:
    summfile = glob.glob(thisfile.rstrip('.xml.gz') + '.txt')
  # check if corse file has a corresponding summary file
  if not summfile:
    print >>sys.stderr, "A summary file for %s could not be found." %(thisfile)
    sys.exit(1)
  # get needed info from summary file
  file = open(summfile[0], 'r')
  for line in file:
    # get coincidence type (used for labeling in plots); will have regardless
    # of whether or not there are foreground triggers
    if line.startswith( 'coincident ifos:' ):
      # if no background triggers mark the corsefile for later removal and
      # break from reading the summary file 
      if line.split()[2] == 'no_background':
        no_bkg[thisfile] = True
      else:
        coincifos[thisfile] = line.split()[2]
    # get Number of background triggers
    elif line.startswith( 'number of reconstructed slide coincidences:' ):
      NbkgCoinc = float( line.split()[5] )
    elif line.startswith( 'number of reconstructed zero-lag coincidences' ):
      if float( line.split(':')[1] ) == 0.: # no foreground
        if thisfile in no_bkg:
          # no foreground or background: remove from the no_bkg dictionary,
          # mark the corsefile for later removal and break from
          # reading the summary file
          no_bkg_frgnd[thisfile] = corsefiles.index(thisfile)
          del no_bkg[thisfile]
          break
        else: # just have no foreground; add to no foreground dict
          no_frgnd[thisfile] = True
    # get foreground time analyzed; used maxbkgFAN as well as for normalizing
    elif line.startswith( 'amount of time analysed for triggers' ):
      pass
    # get background time analyzed
    elif line.startswith( 'amount of background time analyzed' ):
      BkgTime = float( line.split()[5] )
    # get mass-bin (used for labeling)
    # this is an inverted dictionary, i.e., mass-bin is referenced by the
    # masses and the elements are a list of the files that are in that mass
    # bin
    elif line.startswith( 'mass-bin:' ):
      mass = line.split(":")[1].rstrip('\n').lstrip()
      mass = mass.replace( '_','-' )
      if mass in massbin:
        massbin[ mass ].append(thisfile)
      else:
        massbin[ mass ] = []
        massbin[ mass ].append(thisfile)
  file.close()
  if opts.time_correct_file:
    corrfile = open(opts.time_correct_file,'r')
    for line in corrfile:
      if line.startswith( 'amount of time analysed' ):
        FrgrndTime = float( line.split()[6] )
    corrfile.close()
  else:
    raise ValueError, "Must specify --time-correct-file"

  # calculate min/max BkgFANs
  minBkgFAN[thisfile] = FrgrndTime/BkgTime
  if thisfile not in no_bkg_frgnd: # i.e., has background or foreground
    if thisfile not in no_bkg:
      maxBkgFAN[thisfile] = NbkgCoinc * FrgrndTime/BkgTime
  # apply correction factor
    modFAN[thisfile] = 1.
# End: loop over corsefiles

# remove files that had no bkg from corsefiles list and add their names to
# warn_msg
if no_bkg_frgnd:
  warn_msg = 'No foreground or background in files:\n'
  for idx in sorted(no_bkg_frgnd.items(), key=itemgetter(1), reverse=True):
    warn_msg = warn_msg + ' ' + os.path.basename(corsefiles.pop(idx[1])) + '\n'
  # check if still have a corsefiles list; if all the files that were globbed
  # don't have foreground and background, just make a generic plot with
  # warn_msg on it; this avoids future errors
  if not corsefiles:
    warn_msg = warn_msg + 'These were all the globbed files.'
    sys.exit(0)

coincStat = CoincInspiralUtils.coincStatistic("far")
for thisfile in corsefiles:
  insptrigs = SnglInspiralUtils.ReadSnglInspiralFromFiles( [thisfile] )
  coincT[ thisfile ] = CoincInspiralUtils.coincInspiralTable( insptrigs, coincStat )
  coincT[ thisfile ].sort() # sort by descending FAN
  # if this file has no_bkg, but does have foreground, get the ifo coincident
  # type from the first foreground trigger
  if thisfile in no_bkg:
    coincifos[thisfile] = coincT[ thisfile ][0].get_ifos()[0]
#for thisfile in corsefiles:
#  if NormTime[corsefiles[0]] != NormTime[thisfile]:
#    print >> sys.stderr, "Can't combine experiments with " + \
#      "different analysis times."
#    sys.exit( 1 )
maxFANs = [] # for storing max FAN of bkg (the dict is hard to sort by value)
FANc = [] # for storing the combined FANs of foreground triggers
zero_fanc = []
for thisfile in corsefiles:
  if thisfile in maxBkgFAN:maxFANs.append( maxBkgFAN[thisfile] )
  if coincT[thisfile].sngl_table: # if have foreground trigs
    for coinc in coincT[thisfile]:
      if coinc.stat == 0.: # mark and set to minbkgFAN[thisfile]
        zero_fanc.append(minBkgFAN[thisfile])
        FANc.append( minBkgFAN[thisfile])
      else:
        FANc.append( coinc.stat*modFAN[thisfile] )
FANc.sort(reverse=True) # order from weakest to strongest
maxFANs.sort(reverse=True)
FANc = array(FANc) # convert to array
maxFANs = array(maxFANs)
if zero_fanc: zero_fanc = array(zero_fanc)
# following works by starting from highest fan values and moving down (like
# moving left to right on an IFAN plot)
for ii in range(0,len(FANc)):
  for jj in range(1, len(maxFANs)): # cycle through bkg fans, skipping first one
    if FANc[ii] > maxFANs[jj]: # find the largest bkg fan < this foreground fan
      FANc[ii] = FANc[ii] * (jj) # multiply by number of active categories
      for kk in range(jj, len(maxFANs)):
        FANc[ii] = FANc[ii] + maxFANs[kk] # add bkg fans of inactive categories
      break # go to next fan in FANc
    elif jj == len(maxFANs)-1: # all categories active
      FANc[ii] = FANc[ii] * len(maxFANs)

combinedTrigs = lsctables.New( lsctables.SnglInspiralTable , columns=[])
loudestTrig = lsctables.New( lsctables.SnglInspiralTable , columns=[])

if opts.min_rate:
  minIFAN = opts.min_rate/(FrgrndTime/3.15567360E7)
  maxFAN = 1/minIFAN

for column in columnList():
  combinedTrigs.appendColumn(column)
  loudestTrig.appendColumn(column)

loudestTrigTemp = []
loudestTrigFAR = 99999999999.
for thisfile in corsefiles:
  insptrigs = SnglInspiralUtils.ReadSnglInspiralFromFiles( [thisfile] )
  if insptrigs:
    for trig in insptrigs:
      if trig.alpha != 0.:
        trig.alpha = trig.alpha * modFAN[thisfile]
        for jj in range(1, len(maxFANs)):
          if trig.alpha > maxFANs[jj]:
            trig.alpha = trig.alpha * (jj)
            for kk in range(jj, len(maxFANs)):
              trig.alpha = trig.alpha + maxFANs[kk]
            break # go to next fan in FANc
          elif jj == len(maxFANs)-1: # all categories active
            trig.alpha = trig.alpha * len(maxFANs)
      combinedTrigs.append(trig)
      if opts.min_rate:
        if trig.alpha < (maxFAN):	  
          loudestTrig.append(trig)
	elif len(loudestTrig) == 0:
	  if trig.alpha < (loudestTrigFAR-0.000001):
            loudestTrigTemp = [trig]
            loudestTrigFAR = trig.alpha
          elif not trig.alpha > (loudestTrigFAR+0.000001):
            loudestTrigTemp.append(trig)

if opts.min_rate:
  if len(loudestTrig) == 0:
    for trig in loudestTrigTemp:
      loudestTrig.append(trig)
  print >> sys.stdout, 'The time analyzed was' , FrgrndTime
  print >> sys.stdout, 'Accepting triggers with IFAN >' , minIFAN
  print >> sys.stdout, 'Accepting triggers with FAN <' , maxFAN
  print >> sys.stdout, 'Accepted' , len(loudestTrig) , ' SINGLE triggers.'

flag = False
integ = 0
while not flag:
  try:
    outputFile=utils.load_filename(corsefiles[integ],\
       gz=corsefiles[integ].endswith('.gz'))
    origtbl = tab.get_table(outputFile, lsctables.SnglInspiralTable.tableName)
    flag = True
  except:
    integ = integ + 1
parent = origtbl.parentNode
parent.replaceChild(combinedTrigs,origtbl)
utils.write_filename(outputFile,opts.output_file, \
    gz = opts.output_file.endswith('gz'))
if opts.min_rate:
  for trig in loudestTrig:
    trig.alpha = trig.alpha / (FrgrndTime/3.15567360E7)
  newtbl = tab.get_table(outputFile, lsctables.SnglInspiralTable.tableName)
  parent = newtbl.parentNode
  parent.replaceChild(loudestTrig,newtbl)
  utils.write_filename(outputFile,opts.output_file_loudest, \
      gz = opts.output_file.endswith('gz'))

#searchSumm = lsctables.New(lsctables.SearchSummaryTable)
#summVal = lsctables.New(lsctables.SummValueTable)
#
#outputFile = open(opts.output_file,"w")
#outputFile.write('''<?xml version='1.0' encoding='utf-8' ?> \n''')
#outputFile.write('''<!DOCTYPE LIGO_LW SYSTEM "http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt"><LIGO_LW> \n''')
#searchSumm.write(outputFile)
#summVal.write(outputFile)
#combinedTrigs.write(outputFile)
#outputFile.write('''</LIGO_LW>''')
