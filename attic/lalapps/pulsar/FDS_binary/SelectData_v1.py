#!/usr/local/bin/python2.3
"""
This script has been written to scan the results from a sensitivity
analysis of an SFT dataset.  It selects the most sensitive result
and copies the appropriate most sensitive SFT dataset to another location.
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
Usage: SelectData [options]

  -h, --help                    display this messag
  -f, --sourcefile              the name of the source file used
  -s, --source                  the name of the source 
  -S, --sensitivitydir          the location of the sensitivity results
  -D, --fulldatadir             the location of the full dataset
  -O, --optdatadir              the location where the optimal dataset is to be stored
  -o, --optdataparamsfile       the name of the output file containing the optimal datset parameters
  -x, --justcopy                set this flag if you just want to copy the full data as the optimal data 
  """
  print >> sys.stderr, msg

# parse the command line options to figure out what we should do
shortop = "hS:f:s:D:O:o:x"
longop = [
"help",
"sourcefile=",
"source=",
"sensitivitydir=",
"fulldatadir=",
"optdatadir=",
"optdataparamsfile=",
"justcopy"
]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

# default options
sourcefile = None
source = None
sensitivitydir = None
optimaldir = None
fulldatadir = None
optdataparamsfile = None
justcopy = 0

# process options
for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-f", "--sourcefile"):
    sourcefile = a
  elif o in ("-s", "--source"):
    source = a
  elif o in ("-S", "--sensitivitydir"):
    sensitivitydir = a
  elif o in ("-D", "--fulldatadir"):
    fulldatadir = a
  elif o in ("-t", "--optdatadir"):
    optdatadir = a
  elif o in ("-s", "--optdataparamsfile"):
    optdataparamsfile = a
  elif o in ("-x", "--justcopy"):
    justcopy = 1 
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

if not sourcefile:
  print >> sys.stderr, "No sourcefile specified."
  print >> sys.stderr, "Use --sourcefile FILE to specify location."
  sys.exit(1)

if not source:
  print >> sys.stderr, "No source specified."
  print >> sys.stderr, "Use --source STRING to specify name."
  sys.exit(1)

if ((sensitivitydir==None)&(justcopy==0)):
  print >> sys.stderr, "No sensitivity results directory specified."
  print >> sys.stderr, "Use --sensitivitydir PATH to specify location."
  sys.exit(1)

if not fulldatadir:
  print >> sys.stderr, "No full data directory specified."
  print >> sys.stderr, "Use --fulldatadir PATH to specify location."
  sys.exit(1)

if not optdatadir:
  print >> sys.stderr, "No optimal data directory location specified."
  print >> sys.stderr, "Use --optdatadir PATH to specify a location."
  sys.exit(1)

if not optdataparamsfile:
  print >> sys.stderr, "No optimal data parameter file specified."
  print >> sys.stderr, "Use --optdataparamsfile FILE to specify a file."
  sys.exit(1)

###############################################################################
# loop over each file in the sensitivity results directory

if not justcopy:

  # extract the filelist
  filelist = os.listdir(sensitivitydir)
  
  opt_result = 10000
  opt_line = " "
  
  # loop over the files in filelist
  for filename in filelist:
    
    # open each file
    try: fp = open(sensitivitydir + '/' + filename, 'r')
    except:
      print "ERROR : strange, can't open sensitivity results file %s" %(filename)
      sys.exit(1)
      
    # loop over each line in this file  
    for line in fp:
      temp = line.rstrip().split()
      if (temp[0]!="#"):
        sens = math.log10(float(temp[8]))
        if (sens<opt_result):
          opt_result = sens
          opt_line = line

    fp.close()

  # now we extract the relevant parameters from that line
  temp = opt_line.rstrip().split()
  start = (int)(temp[0])
  end = (int)(temp[1])

else:

  # if we are just copying the fulldataset to the optimal dataset
  start = 0
  end  = 999999999
  opt_line = "NO SENSITIVITY RESULTS : the full data set has been copied to the optimal data set"

#################################################################################
# extract the filelist from the full data directory
fulldatafilelist = os.listdir(fulldatadir)

# now loop over the files in the filelist
for filename in fulldatafilelist:

  # extract the time from the filename
  time = (int)(re.search('\d{9}', filename).group())

  # if the file falls into the optimal time range then soft link it
  if ((time>=start)&(time<end)):
    try: os.symlink(fulldatadir + '/' + filename, optdatadir + '/' + filename)
    except:
      print "ERROR : unable to make soft link from full data to optimal data"
      sys.exit(1)

#################################################################################
# now output the optimal dataset parameters to file

try: fp = open(optdataparamsfile,'w')
except:
  print "ERROR : failed to open optimal parameter output file %s" %(optdataparamsfile)
  sys.exit(1)

fp.write("# optimal data set parameters file\n")
fp.write("# ------------------------------------------------------------------------------------------------------------------------------------------- #\n")
fp.write("# sourcefile = %s\n" %(sourcefile))
fp.write("# source = %s\n" %(source))
fp.write("# full data directory = %s\n" %(fulldatadir))
fp.write("# ------------------------------------------------------------------------------------------------------------------------------------------- #\n")
fp.write("# tstart          tend            Sh_1            Sh_2            Sh_wtd          A               B               T_obs           Q\n")
fp.write("#\n")
fp.write("%s\n" %(opt_line))

fp.close()


sys.exit(0)

