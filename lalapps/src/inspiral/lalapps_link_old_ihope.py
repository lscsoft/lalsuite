"""

Try to set up files from an old run within the directory hierarchy of
a new run

"""

__version__ = "$Revision$"
__date__ = "$Date$"
__name__ = "plotinspiral"
__Id__ = "$Id$"
__title__ = "Inspiral Plots"


import sys
import os
from optparse import *
import re
import exceptions
import glob
from types import *

from glue import lal

#################################################################
# help message
usage = """\
%prog [options]
------------------------------------------------------------------------------
  
"""

#################################################################
"""
Parser function dedicated
"""
parser = OptionParser( usage=usage, \
    version= "%prog CVS\n" +
    "$Id$\n" +
    "$Name$\n")

#action related
parser.add_option("-a","--tmpltbank",action="store_true",default=False,\
    help="look for template bank files that you can make links to" )
parser.add_option("-b","--inspiral-first",action="store_true",default=False,\
    help="look for first inspiral files that you can make links to)" )
parser.add_option("","--inspiral-second",action="store_true",default=False,\
    help="look for second stage inspiral files that you can make links to)" )
parser.add_option("","--thinca-first",action="store_true",default=False,\
    help="look for first thinca files that you can make links to)" )
parser.add_option("","--thinca-second",action="store_true",default=False,\
    help="look for second stage thinca files (CAT_1 vetoes) that you " +
     "can make links to)" )
parser.add_option("","--thinca-second-veto",action="store_true",default=False,\
    help="look for second stage thinca files (CAT_2/3/4 vetoes) that you " +
    "can make links to)" )
parser.add_option("","--inspinj",action="store_true",default=False,\
    help="look for inspinj files that you can make links to)" )
parser.add_option("","--inca",action="store_true",default=False,\
    help="look for inca files that you can make links to)" )
parser.add_option("","--sire-first",action="store_true",default=False,\
    help="look for first sire files that you can make links to)" )
parser.add_option("","--sire-second",action="store_true",default=False,\
    help="look for second stage sire files (CAT_1 vetoes) that you " +
    "can make links to)" )
parser.add_option("","--sire-second-veto",action="store_true",default=False,\
    help="look for second stage sire files (CAT_2/3/4 vetoes) that you " +
    "can make links to)" )
parser.add_option("","--coire-first",action="store_true",default=False,\
    help="look for first coire files that you can make links to)" )
parser.add_option("","--coire-second",action="store_true",default=False,\
    help="look for second stage coire files (CAT_1 vetoes) that you " +
    "can make links to)" )
parser.add_option("","--coire-second-veto",action="store_true",default=False,\
    help="look for second stage coire files (CAT_2/3/4 vetoes) that you " +
    "can make links to)" )
parser.add_option("","--trigbank",action="store_true",default=False,\
    help="look for trigbank files that you can make links to)" )
parser.add_option("-c","--make-links",action="store_true",default=False,\
    help="make the links" )
#input
parser.add_option("-o", "--old-cache-file",action="store",type="string",\
    default="", metavar="OLD_CACHE_FILE_NAME", help="read filenames from "\
    "this old cache file. This is the cache file that contains the files "\
    "you want to link to.")
parser.add_option("-n", "--append-cache",action="store",type="string",\
    default="", metavar="NEW_CACHE_FILE_NAME", help="append the newly "\
    "linked files to this cache file." )
parser.add_option("-r", "--run-names",action="store",\
    type="string",default="", metavar="RUN_NAMES",\
    help="names for the different type of runs. Should "\
    "be a comma separated list with no spaces. For instance, "\
    "--run-names playground,full_data,inj001,inj002,inj003")
(opts,args) = parser.parse_args()

# load cache
if opts.old_cache_file is not None:
  cache = lal.Cache.fromfile(open(opts.old_cache_file))
else:
  print "Must specify a cache file"
  exit(1)

if not opts.run_names:
  print "Must specify the names of the runs. For example, --run-names "\
    "playground,inj001,inj002"
  exit(1)

# Split the names of the runs into a list. On the command line,
# the list should be comma-separated. (playground,inj001,etc.)
runs = opts.run_names.split(",")

# Make a list of all of the file types that you'd like to link
# Example: ["TMPLTBANK_PLAYGROUND", "INSPIRAL_FIRST_FULL_DATA"]
filetypes = []

for run in runs:
  runCaps = run.upper()
  if opts.tmpltbank:
    filetypes.append("TMPLTBANK_" + runCaps)
  if opts.inspiral_first:
    filetypes.append("INSPIRAL_FIRST_" + runCaps)
  if opts.inspiral_second:
    filetypes.append("INSPIRAL_SECOND_*" + runCaps)
  if opts.thinca_first:
    filetypes.append("THINCA*FIRST_" + runCaps)
  if opts.thinca_second:
    filetypes.append("THINCA*SECOND_*" + runCaps)
  if opts.thinca_second_veto:
    filetypes.append("THINCA*SECOND_*" + runCaps + "*CAT*VETO")
  if opts.sire_first:
    filetypes.append("SIRE*FIRST_" + runCaps)
  if opts.sire_second:
    filetypes.append("SIRE*SECOND_*" + runCaps)
  if opts.sire_second_veto:
    filetypes.append("SIRE*SECOND_*" + runCaps + "*CAT*VETO")
  if opts.coire_first:
    filetypes.append("COIRE*FIRST_" + runCaps) 
  if opts.coire_second:
    filetypes.append("COIRE*SECOND_*" + runCaps)
  if opts.coire_second_veto:
    filetypes.append("COIRE*SECOND_*" + runCaps + "*CAT*VETO")
  if opts.trigbank:
    filetypes.append("TRIGBANK*" + runCaps)
  if opts.inspinj:
    filetypes.append("INJECTIONS*" + runCaps)
  if opts.inca:
    filetypes.append("INCA_" + runCaps)

# Start a list that will contain the name of every linked file.
# The files in this list will later be appended to the new cache file.
filelist = []

# Loop over filetypelist (each item in filetypelist is a sieve pattern).
# Simultaneously loop over runs, which has the corresponding 
# directory for each sieve pattern. In order to match up the directory
# names, the items in "runs" must be repeated (hence the multiplication
# by len(filetypes)).

directoryList = []
for run in runs:
  for i in range(0,len(filetypes)//len(runs)):
    directoryList.append(run)

for filetype,filedir in zip(filetypes,directoryList):
  tmpcache = cache.sieve(description=filetype,exact_match=True)
  # Check if files from the old cache exist on disk and place 
  # them in the found and missed lists.
  (tmpfound,tmpmissed) = tmpcache.checkfilesexist()
  if len(tmpmissed)>0:
    print "Warning: These files do not exist on disk, but are listed in the cache."
    for missingfile in tmpmissed:
      print "Missing file: " + str(missingfile)
  if len(tmpfound)==0:
    print "Nothing matched the sieve pattern: "+filetype
  # Make the links and append the new file name to filelist
  if opts.make_links and len(tmpfound)>0:
    print "Linking the " + filetype + " files..."
    for file in tmpfound.pfnlist():
      os.symlink(file,filedir+"/"+os.path.basename(file))
      filelist.append(filedir+"/"+os.path.basename(file))

# Append the newly linked files to the cache file 
if opts.append_cache:
  outcache = cache.from_urls(filelist)
  outcache.tofile(open(opts.append_cache,"a"))
