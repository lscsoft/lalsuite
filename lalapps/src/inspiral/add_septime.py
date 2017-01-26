#!/usr/bin/python

import sys
from optparse import *
from types import *

import matplotlib
matplotlib.use('Agg')
from pylab import *
rc('text', usetex=True)

##############################################################################
usage = """
usage: %prog [options] 

"""

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage, version="%prog CVS $Id$ " )

  # options related to input and output
  parser.add_option("","--input-file",action="store",type="string",\
      default=None, metavar=" FILE",help="ANTIME files to read" )

  parser.add_option("","--output-file",action="store",type="string",\
      default=None,metavar=" FILE",\
      help="file to write to")

  parser.add_option("","--num-slides",action="store",type="int",default=0,\
      metavar=" NUM_SLIDES",help="number of time slides performed" )
  
  parser.add_option("","--set-zero-lag-time-to",action="store",type="int",default=0,\
      metavar=" NUM_SLIDES",help="sets zero lag to a unit time (e.g. 1 year) in seconds" )

 
  (options,args) = parser.parse_args()


  return options, sys.argv[1:]


# ============================================================================
# -- get command line arguments
opts, args = parse_command_line()

if not opts.input_file:
  print >>sys.stderr, "Must specify a file in --input-file to read"
  print >>sys.stderr, "Enter 'add_septime --help' for usage"
  sys.exit(1)


###################################
# glob the list of files to read in
if opts.input_file is not None:
  allfiles = []
  fp = open(opts.input_file, 'r')
  for file in fp.readlines():
    allfiles.append(file.replace('\n',''))
  fp.close()
  if len(allfiles) < 1:
    print >>sys.stderr, "The glob for " + opts.glob + " returned no files" 
    sys.exit(1)

##################################
# setup structure to hold the analyzed times
analyzedTimes = {}
slides = range(opts.num_slides + 1)
slides.extend(range(-opts.num_slides, 0))

for slide in slides:
  analyzedTimes[str(slide)] = 0.0

##################################
# loop over the input files and add up the time analyzed
for file in allfiles:
  fp = open(file, 'r')
  times = {}
  for line in fp.readlines():
    slide,time = line.split(" ")
    times[slide] = float(time)
  fp.close()

  for key in analyzedTimes.keys():
    analyzedTimes[key] += times[key]

##################################
# plot the output
if opts.output_file:
  plotName = (opts.output_file.split("."))[0] + '.png'

  figure(1)
  for key in analyzedTimes.keys():
    if key == '0':
      color='r'
    else:
      color='b'
    bar(left=float(key), height=analyzedTimes[key], width=1, color=color,
        align="center")
    hold(True)
  ylabel('Analyzed Time (sec)', size='x-large')
  xlabel('Slide Number', size='x-large')
  savefig(plotName)

##################################

if opts.set_zero_lag_time_to:
  analyzedTimes["0"] = opts.set_zero_lag_time_to

#  write output file
if opts.output_file:
  fp = open(opts.output_file, 'w')
  for slide in slides:
    fp.write('%i %f\n' % (slide,analyzedTimes[str(slide)]))
  fp.close()
else:
  for slide in slides:
    print slide,analyzedTimes[str(slide)]

