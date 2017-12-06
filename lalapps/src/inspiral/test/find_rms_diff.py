#!/usr/bin/python
"""
Run a unit test
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
#from types import *

from pylal import SnglInspiralUtils

#################################################################
# help message
usage = """\
%prog [options]
"""

def getrmsdiff(ref_file_glob, test_file_glob, columns):

  refFiles = []
  testFiles = []

  #print "I'm in the function and the globs are %s and %s" %(ref_file_glob,test_file_glob)

  refFiles.extend(glob.glob(ref_file_glob))
  testFiles.extend(glob.glob(test_file_glob))

  refTrigs = SnglInspiralUtils.ReadSnglInspiralFromFiles(refFiles)
  testTrigs = SnglInspiralUtils.ReadSnglInspiralFromFiles(testFiles)

  if len(refTrigs) != len(testTrigs):
    #print "unequal number of triggers: %d %d" % (len(refTrigs), len(testTrigs))
    meandiff = ["NaN"]
  else:
    meandiff = []
    for c in columns:
      refcolumn = refTrigs.get_column(c)
      testcolumn = testTrigs.get_column(c)
      tmpmeandiff = sqrt( mean ( pow(refcolumn-testcolumn,2.0)))
      meandiff.append(tmpmeandiff)
    #print "%s %s" % (c,meandiff)

  return len(refTrigs), len(testTrigs), meandiff

def filloutput(reflen, testlen, cols, rmsdiff, type, ifo):

  # Name is either "triggers" or "templates" depending on the file type
  if type == "TMPLTBANK" or type == "TRIGBANK":
    name = "templates"
  else:
    name = "triggers"

  if reflen != testlen:
    outfile.write("%s %s: There are %i %s in the reference and %i %s in the test. No comparison possible.\n" %(ifo, type, reflen, name, testlen, name))
    pffile.write("F")
  else:
    outfile.write("%s %s: %i %s\n" % (ifo, type, reflen, name))
    pffile.write(".")
    for param,value in zip(cols,rmsdiff):
      outfile.write("%s %s\n" %(param,value))
      if value < 1.e-3:
        pffile.write(".")
      else:
        pffile.write("F")

parser = OptionParser( usage=usage, \
    version= "%prog CVS\n" +
    "$Id$\n" +
    "$Name$\n")

parser.add_option("-r","--ref-dir",action="store",type="string",\
    default=None, metavar=" GLOB",help="reference directory to read" )
parser.add_option("-t","--test-dir",action="store",type="string",\
    default=None, metavar=" GLOB",help="test directory to read" )
parser.add_option("-u","--unit-test",action="store_true",default=False,\
    help="Specify this option if a unit-test is being performed." \
    "If this option is used, you may not specify the ref-dir "\
    "or test-dir")
parser.add_option("-o","--output-file",action="store",type="string",\
    metavar=" OUTPUT",help="output file name")

(opts,args) = parser.parse_args()

import matplotlib
matplotlib.use('Agg')
from pylab import *

if (opts.ref_dir == None or opts.test_dir == None) and not opts.unit_test:
  print "You must supply a reference directory and a test directory or "\
  "give the unit_test option which will set the directories for you."
  sys.exit(1)

if (opts.ref_dir or opts.test_dir) and opts.unit_test:
  print "You cannot give the ref or test directory when the "\
  "unit test is specified."
  sys.exit(1)

if opts.output_file == None:
  output = sys.stdout
else:
  output = opts.output_file

######################################################
# SET THE DIRECTORIES

if opts.ref_dir:
  refdir = opts.ref_dir
else:
  refdir = "/scratch3/jclayton/run_cbc_s5_1yr_20070129/testinj_lvmonth2/866088014-868721414/anotherunittest/"

if opts.test_dir:
  testdir = opts.test_dir
else:
  testdir = os.getcwd()

######################################################
# OPEN THE OUTPUT FILE
outfile = open(output,'w')

outfile.write("REFERENCE DIRECTORY = %s\n" %(refdir))
outfile.write("TEST DIRECTORY = %s\n" %(testdir))
outfile.write("\n")

# OPEN THE PASS/FAIL FILE
pffile = open("passfail.txt","w")

#######################################################
# COMPARE THE TEMPLATE BANK FILES

outfile.write("**********************************\n")
outfile.write("Compare the template bank files...\n")

reftmpltH1glob = "%s/H1-TMPLTBANK_INJ*" %(refdir)
reftmpltH2glob = "%s/H2-TMPLTBANK_INJ*" %(refdir)
reftmpltL1glob = "%s/L1-TMPLTBANK_INJ*" %(refdir)

testtmpltH1glob = "%s/H1-TMPLTBANK_INJ*" %(testdir)
testtmpltH2glob = "%s/H2-TMPLTBANK_INJ*" %(testdir)
testtmpltL1glob = "%s/L1-TMPLTBANK_INJ*" %(testdir)

columns = ["mass1","mass2","mchirp"]

lenH1ref_tmplt, lenH1test_tmplt, diff_H1tmplt = getrmsdiff(reftmpltH1glob,testtmpltH1glob,columns)
lenH2ref_tmplt, lenH2test_tmplt, diff_H2tmplt = getrmsdiff(reftmpltH2glob,testtmpltH2glob,columns)
lenL1ref_tmplt, lenL1test_tmplt, diff_L1tmplt = getrmsdiff(reftmpltL1glob,testtmpltL1glob,columns)

filloutput(lenH1ref_tmplt, lenH1test_tmplt, columns, diff_H1tmplt, "TMPLTBANK", "H1")
filloutput(lenH2ref_tmplt, lenH2test_tmplt, columns, diff_H2tmplt, "TMPLTBANK", "H2")
filloutput(lenL1ref_tmplt, lenL1test_tmplt, columns, diff_L1tmplt, "TMPLTBANK", "L1")

#######################################################
# COMPARE THE INSPIRAL FIRST FILES

outfile.write("***********************************\n")
outfile.write("Compare the inspiral first files...\n")

refinspiral1H1glob = "%s/H1-INSPIRAL_FIRST_INJ001*gz" %(refdir)
refinspiral1H2glob = "%s/H2-INSPIRAL_FIRST_INJ001*gz" %(refdir)
refinspiral1L1glob = "%s/L1-INSPIRAL_FIRST_INJ001*gz" %(refdir)

testinspiral1H1glob = "%s/H1-INSPIRAL_FIRST_INJ001*gz" %(testdir)
testinspiral1H2glob = "%s/H2-INSPIRAL_FIRST_INJ001*gz" %(testdir)
testinspiral1L1glob = "%s/L1-INSPIRAL_FIRST_INJ001*gz" %(testdir)

columns = ["mass1","mass2", "mchirp", "snr"]

lenH1ref_inspiral1, lenH1test_inspiral1, diff_H1inspiral1 = getrmsdiff(refinspiral1H1glob,testinspiral1H1glob,columns)
lenH2ref_inspiral1, lenH2test_inspiral1, diff_H2inspiral1 = getrmsdiff(refinspiral1H2glob,testinspiral1H2glob,columns)
lenL1ref_inspiral1, lenL1test_inspiral1, diff_L1inspiral1 = getrmsdiff(refinspiral1L1glob,testinspiral1L1glob,columns)

filloutput(lenH1ref_inspiral1, lenH1test_inspiral1, columns, diff_H1inspiral1, "INSPIRAL FIRST", "H1")
filloutput(lenH2ref_inspiral1, lenH2test_inspiral1, columns, diff_H2inspiral1, "INSPIRAL FIRST", "H2")
filloutput(lenL1ref_inspiral1, lenL1test_inspiral1, columns, diff_L1inspiral1, "INSPIRAL FIRST", "L1")

#######################################################
# COMPARE THE THINCA FIRST FILES

outfile.write("***********************************\n")
outfile.write("Compare the thinca first files...\n")

refthinca1glob = "%s/H1H2L1-THINCA_FIRST*" %(refdir)

testthinca1glob = "%s/H1H2L1-THINCA_FIRST*" %(testdir)

columns = ["mass1","mass2","mchirp","snr"]

lenref_thinca1, lentest_thinca1, diff_thinca1 = getrmsdiff(refthinca1glob,testthinca1glob,columns)

filloutput(lenref_thinca1, lentest_thinca1, columns, diff_thinca1, "THINCA FIRST", "H1H2L1")

#######################################################
# COMPARE THE TRIGBANK FILES

outfile.write("**********************************\n")
outfile.write("Compare the trig bank files...\n")

reftrigbankH1glob = "%s/H1-TRIGBANK_H1H2L1_*" %(refdir)
reftrigbankH2glob = "%s/H2-TRIGBANK_H1H2L1_*" %(refdir)
reftrigbankL1glob = "%s/L1-TRIGBANK_H1H2L1_*" %(refdir)

testtrigbankH1glob = "%s/H1-TRIGBANK_H1H2L1_*" %(testdir)
testtrigbankH2glob = "%s/H2-TRIGBANK_H1H2L1_*" %(testdir)
testtrigbankL1glob = "%s/L1-TRIGBANK_H1H2L1_*" %(testdir)

columns = ["mass1","mass2","mchirp"]

lenH1ref_trigbank, lenH1test_trigbank, diff_H1trigbank = getrmsdiff(reftrigbankH1glob,testtrigbankH1glob,columns)
lenH2ref_trigbank, lenH2test_trigbank, diff_H2trigbank = getrmsdiff(reftrigbankH2glob,testtrigbankH2glob,columns)
lenL1ref_trigbank, lenL1test_trigbank, diff_L1trigbank = getrmsdiff(reftrigbankL1glob,testtrigbankL1glob,columns)

filloutput(lenH1ref_trigbank, lenH1test_trigbank, columns, diff_H1trigbank, "TRIGBANK", "H1")
filloutput(lenH2ref_trigbank, lenH2test_trigbank, columns, diff_H2trigbank, "TRIGBANK", "H2")
filloutput(lenL1ref_trigbank, lenL1test_trigbank, columns, diff_L1trigbank, "TRIGBANK", "L1")

#######################################################
# COMPARE THE INSPIRAL SECOND FILES

outfile.write("***********************************\n")
outfile.write("Compare the inspiral second files...\n")

refinspiral2H1glob = "%s/H1-INSPIRAL_SECOND_H1H2L1_*" %(refdir)
refinspiral2H2glob = "%s/H2-INSPIRAL_SECOND_H1H2L1_*" %(refdir)
refinspiral2L1glob = "%s/L1-INSPIRAL_SECOND_H1H2L1_*" %(refdir)

testinspiral2H1glob = "%s/H1-INSPIRAL_SECOND_H1H2L1_*" %(testdir)
testinspiral2H2glob = "%s/H2-INSPIRAL_SECOND_H1H2L1_*" %(testdir)
testinspiral2L1glob = "%s/L1-INSPIRAL_SECOND_H1H2L1_*" %(testdir)

columns = ["mass1","mass2", "mchirp","snr","chisq"]

lenH1ref_inspiral2, lenH1test_inspiral2, diff_H1inspiral2 = getrmsdiff(refinspiral2H1glob,testinspiral2H1glob,columns)
lenH2ref_inspiral2, lenH2test_inspiral2, diff_H2inspiral2 = getrmsdiff(refinspiral2H2glob,testinspiral2H2glob,columns)
lenL1ref_inspiral2, lenL1test_inspiral2, diff_L1inspiral2 = getrmsdiff(refinspiral2L1glob,testinspiral2L1glob,columns)

filloutput(lenH1ref_inspiral2, lenH1test_inspiral2, columns, diff_H1inspiral2, "INSPIRAL SECOND", "H1")
filloutput(lenH2ref_inspiral2, lenH2test_inspiral2, columns, diff_H2inspiral2, "INSPIRAL SECOND", "H2")
filloutput(lenL1ref_inspiral2, lenL1test_inspiral2, columns, diff_L1inspiral2, "INSPIRAL SECOND", "L1")

#######################################################
# COMPARE THE THINCA SECOND FILES

outfile.write("***********************************\n")
outfile.write("Compare the thinca second files...\n")

refthinca2glob = "%s/H1H2L1-THINCA_SECOND*" %(refdir)

testthinca2glob = "%s/H1H2L1-THINCA_SECOND*" %(testdir)

columns = ["mass1","mass2", "mchirp", "snr","chisq"]

lenref_thinca2, lentest_thinca2, diff_thinca2 = getrmsdiff(refthinca2glob,testthinca2glob,columns)

filloutput(lenref_thinca2, lentest_thinca2, columns, diff_thinca2, "THINCA SECOND", "H1H2L1")

outfile.write("*************************************\n")

pffile.write("\n")

