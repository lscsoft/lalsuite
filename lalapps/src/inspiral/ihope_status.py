#!/usr/bin/python 
__author__ = 'Thomas Cokelaer < Thomas.Cokelaer@astro.cf.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'
__name__ = 'ihope_status'


import sys
import os
import optparse
from optparse import OptionParser
import time


usage = """\
%prog [options]
----------------
 brief status of the dag
   
 Purpose - parse one or several dagman.out file(s) to obtain a brief  summary 
 of the dag status.
 
 Example - first, you can parse a dag that contains several dag such as ihope.dag 

         >>> python %prog --ihope-dag-file ihope.dag
         
         - second,  if no sub-dag files are found within the file provide, the 
         dagman.out is parsed %prog --dag-file hipe.dag')

"""
 
# let us parse the command line

parser = OptionParser(usage=usage, \
   version="%prog CVS\n" + 
   "$Id$\n" +
   "$Name$\n")
   

parser.add_option( "--dag-file", dest='scandag', default='ihope.dag', metavar='DAGFILE', help="The dag file to parse")
parser.add_option( "--parse-n-lines", dest='nlines', default='1000', metavar='NLINES', help="The number of lines to parse (starting from the end)")
parser.add_option( "--failed-status", dest='failed', default='1000', metavar='NLINES', help="Give a summary of the failed jobs")

(options, args) = parser.parse_args()


""" There are 2 possibilities. 
    First, the file provided is a dag file that ends in '.dag' that
    contains other dag calls such as in ihope.
    Second, it is a pure dag file which might not even contain a dag
    by itself, and in such case we search for the .dagman.out file.
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
  
print 'Parsing the dag files for status...\n'
### we found some dag files hopefully

time.sleep(1)


for i in xrange(0,len(filenames), 1):
  tag = filenames[i] 
  filename = tag +'.dagman.out'
  print '-->  '+ tag +' status '
  print '------------------------------------------------------------------------'
  os.system('tail -n ' + options.nlines + ' '  + filename +' | grep -v macro | grep -v Note | grep -v Event| grep -v Number | grep -v Node | grep -v Of | grep -v Submit | grep -v node | grep -v submit | grep -v assigned | grep -v seconds | grep -v failed |tail -n 3 - ')
  print '\n\n'

  os.system('grep ERROR '+filename + '> ' + filename + '.status')
  
  




