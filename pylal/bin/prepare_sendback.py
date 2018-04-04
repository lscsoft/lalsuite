#!/usr/bin/env python
"""
Manager program to run specific parts of the followup on the CCIN2P3 cluster
"""

__author__ = 'Damir Buskulic <buskulic@lapp.in2p3.fr>'

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, copy, math, random
from os import *
from optparse import *
import ConfigParser
from UserDict import UserDict
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

parser.add_option("-o", "--output-file",action="store",type="string",\
    help="name of the tar.gz output file, omit the .tar.gz extension. Leave blank for default")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#########  cleaning result files, remove all image .png files  ############
#########     (the image files should have been copied in the  ############
#########      web directory and are not lost)                 ############

depIfoDir = './'

if not opts.output_file:
   outputFile = 'V1_qscans_results'
else:
   outputFile = opts.output_directory

print 'Cleaning RESULTS directory, this may take a while...'

os.system('find '+depIfoDir+'RESULTS -name "*.png" | xargs rm -f')

os.system('cd '+depIfoDir+'../; cp -r V1_qscans_config '+outputFile+'; tar zcvf '+outputFile+'.tar.gz '+outputFile)

print '***'
print '***  Prepared the file '+outputFile+'.tar.gz,'
print '***  you can send it back for further followup analysis'
print '***'

##########################################################################
sys.exit(0)
