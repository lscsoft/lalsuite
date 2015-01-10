#!/usr/local/bin/python2.3
"""
This script has been written to extract the GPS start time
of SFT files from the filenames of SFTs given a directory.
It also returns the total observation time and observation
span plus the number of SFTs in the directory.  If a timestamps
filename is supplied then it produces a timestamps file.
"""

# import required modules
import sys
import os
import getopt
import re
import string
import tempfile
import math
import ConfigParser

# program usage
def usage():
  msg = """\
Usage: GenerateFreqMeshFile [options]

  -h, --help               display this message
  -D, --datadir            the location of the SFT files
  -t, --tsft               the time baseline of the SFT files
  -s, --stampsfile         the file to output the timestamps to (optional)
  -o, --outfile            the file to output the dataset parameters to
  """
  print >> sys.stderr, msg

# parse the command line options to figure out what we should do
shortop = "hD:t:s:o:"
longop = [
"help",
"datadir=",
"tsft=",
"stampsfile=",
"outfile=",
]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

# default options
datadir = None
tsft = None
stampsfile = None
do_stamps = 0
outfile = None

# process options
for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-D", "--datadir"):
    datadir = a
  elif o in ("-t", "--tsft"):
    tsft = (int)(a)
  elif o in ("-s", "--stampsfile"):
    stampsfile = a
    do_stamps = 1
  elif o in ("-o", "--outfile"):
    outfile = a
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

if not datadir:
  print >> sys.stderr, "No data directory specified."
  print >> sys.stderr, "Use --datadir PATH to specify location."
  sys.exit(1)

if not tsft:
  print >> sys.stderr, "No SFT time baseline specified."
  print >> sys.stderr, "Use --tsft INT4 to specify a source."
  sys.exit(1)

if not outfile:
  print >> sys.stderr, "No data parameters outfile specified."
  print >> sys.stderr, "Use --outfile FILE to specify a value."
  sys.exit(1)

###############################################################################
# loop over each file in the SFT directory

# extract the filelist
filelist = os.listdir(datadir)

# initialise the min and max times
earliest = (int)(re.search('\d{9}', filelist[0]).group())
latest = (int)(re.search('\d{9}', filelist[0]).group())
nfiles = 0

# define the timestamps list if we need to
if do_stamps:
  stamp_time = []

# loop over the SFTs in filelist
for filename in filelist:
  
  # this finds any 9 digit numeric sequences in a filename and stores it in time
  time = (int)(re.search('\d{9}', filename).group())

  # find the earliest time
  if (time<earliest):
    earliest = time

  # find the latest time
  if (time>latest):
    latest = time

  # if we are making a stamps file then store time for sorting
  if do_stamps:
    stamp_time.append(time)

  # increment the file counter
  nfiles = nfiles + 1

###############################################################################
# now calculate the dataparams and output to file

start = earliest
end = latest + tsft
tobs = tsft * nfiles
tspan = end - start

# output to the output file
try: of = open(outfile, "w")
except:
  print "ERROR : unable to open output file %s" %(outfile)
  sys.exit(1)

# output to file
of.write("[dataparams]\n")
of.write("start = %d\n" %(start))
of.write("end = %d\n" %(end))
of.write("tobs = %d\n" %(tobs))
of.write("tspan = %d\n" %(tspan))
of.write("nsft = %d\n" %(nfiles))
of.write("datadir = %s\n" %(datadir))

of.close()

###############################################################################
# If making a stamps file then open that file for output

# output to the stamps file
if do_stamps:
  try: of = open(stampsfile, "w")
  except:
    print "ERROR : unable to open stamps file %s" %(stampsfile)
    sys.exit(1)

# sort the stamps
stamp_time.sort()

for i in stamp_time:
  of.write("%d 0\n" %(i))

of.close()

sys.exit(0)

