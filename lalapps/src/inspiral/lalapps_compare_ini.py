"""
Compare two ini files and print out their differences
"""

import sys
from optparse import *
import ConfigParser

usage = """usage: %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-c", "--config-file",action="store",type="string",\
    metavar=" FILE",help="use configuration file FILE")
parser.add_option("-r", "--ref-file",action="store",type="string",\
    metavar=" FILE",help="use reference file FILE")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

# Set the configuration file to cp
cp = ConfigParser.ConfigParser()
cp.read(opts.config_file)

# Set the reference file to cpref
cpref = ConfigParser.ConfigParser()
cpref.read(opts.ref_file)

# Get the sections in the config file
configsections = set(cp.sections())

# Get the sections in the reference file
refsections = set(cpref.sections())

# Make a set of sections that appear in the reference file, but not the config file
missingsections = refsections-configsections

# Make a set of sections that appear in the config file, but not the reference file
extrasections = configsections-refsections

# Make a set of sections that appear in both the config and reference file
commonsections = configsections & refsections

# Print a list of sections that do not appear in the config file
for section in missingsections:
  print "MISSING items in " + opts.config_file + ": [" + section + "] and options " + str(cpref.items(section))

# Print a list of sections that appear in the config file, but not in the reference
for section in extrasections:
  print "EXTRA items in " + opts.config_file + ": [" + section + "] and options " + str(cp.items(section))

# Loop over the sections that appear in both files

for section in commonsections:
  # Make a dictionary that contains the option and option values
  configdict = dict(cp.items(section))
  refdict = dict(cpref.items(section))

  # Make a list of options that appear in the reference file under this section
  refopts = set(cpref.options(section))
  # Make a list of options that appear in the config file under this section
  configopts = set(cp.options(section))

  # Find the missing, extra and common options for this section
  missingopts = refopts - configopts
  extraopts = configopts - refopts
  commonopts = refopts & configopts

  # Print a list of options missing from the config file
  for opt in missingopts:
    print "MISSING option in " + opts.config_file + " in [" + section + "]: " + opt + " = " + refdict[opt]  

  # Print a list of extra options in the config file that aren't in the reference
  for opt in extraopts:
    print "EXTRA option in " + opts.config_file + " in [" + section + "]: "+ opt + " = " + configdict[opt]

  # Compare the option values between the two files and print a FAILURE message if they don't match
  for opt in commonopts:
    configvalue = configdict[opt]
    refvalue = refdict[opt]
    if configvalue != refvalue:
      print "FAILED MATCH: [" + section +"]: " + opt + " ; REFERENCE = " + refvalue + ", CONFIG FILE = " + configvalue
