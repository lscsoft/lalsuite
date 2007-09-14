#!/usr/bin/python 
__author__ = 'Thomas Cokelaer < Thomas.Cokelaer@astro.cf.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'


import sys
import os
import optparse
from optparse import OptionParser
import time


name = 'dagStatus'

def help():
  print name + '\n' 
#+ __version__ + \'n' + __date__
  print __author__
  print """
  [Purpose]
  Parse a dag file and dagman.out files to obtain a brief summary of the dag status

  [Usage]

  -h, --help : this help
  -D, --dag-ihope-file : the input dag file
  -I, --ini-ihope-file : the ini file to parse the filename
   """


# let us parse the command line
parser = OptionParser('%prog [options]\n' +__version__+'\n' + __author__ + '\n\n[Purpose] Parse one or several dagman.out file(s) to obtain a brief status summary of the dag run(s).\n\n[Example] First, you can parse a dag tht contains several dag such as ihope.dag : python %prog --ihope-dag-file ihope.dag\n If no sub-dag files are found within the file provide, the dagman.out is parsed %prog --dag-file hipe.dag')
parser.add_option( "--dag-file", dest='scandag', default='ihope.dag', metavar='DAGFILE', help="The dag file to parse")

(options, args) = parser.parse_args()


""" There are 2 possibilities. 
    First, the file provided is a dag file that ends in '.dag' that
    contains other dag calls such as in ihope
    Second, it is a pure dag file which might not even contain a dag by itself.

"""
print 'Parsing ' +options.scandag + '...'
time.sleep(1)
dagFile = open(options.scandag, 'r')
line = dagFile.readline()
filenames = []
while line:
  if line.find('SCRIPT') >=0:
    if line.endswith('\n'):
      line = line.strip()
    thisline = line.split(' ')
    # we assume that the last field is the one we want. 
    this= thisline[len(thisline)-1]

    filenames.append(this)
    ## we can check this assumption by searching for the string .dag
    if this.find('dag')>0:
      print 'Found ' + this 
  line = dagFile.readline()
dagFile.close()

if len(filenames)==0:
  print 'No dag files found in ' + options.scandag
  print 'Assuming that ' + options.scandag + 'is the dag file you want to parse'
  filenames.append(options.scandag)
  
print 'Parsing the dag files for status...'
### we found some dag files hopefully

time.sleep(1)


for i in xrange(0,len(filenames), 1):
  tag = filenames[i] 
  filename = tag +'.dagman.out'
  print filename
  print '=='+ tag +' status =='
  print '--------------------------'
  os.system('tail -n 100  '+filename+' | grep -v macro | grep -v Note | grep -v Event| grep -v Number | grep -v Node | grep -v Of | grep -v Submit | grep -v node | grep -v submit | grep -v assigned | tail -n 3 - ')
