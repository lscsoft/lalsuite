"""
ihope_status.in - give a brief dag status

$Id $

Reads in a dag file and looks for sub dags. The status of each sub dag is
reported as well as the total number of jobs, the number of submitted jobs
and the progress.

This is specificly designed for use with lalapps_ihope.

"""
__author__ = 'Thomas Cokelaer < Thomas.Cokelaer@astro.cf.ac.uk>, David Mckechan'
__date__ = '$Date$'
__version__ = '$Revision$'
__name__ = 'ihope_status'

import os
import sys
from optparse import OptionParser
from time import sleep


##############################################################################
usage = """\
%prog [options]
---------------- 
  Specifically designed for use with lalapps_ihope uber-dags.
   
 Purpose - parse  several dags and their dagman.out file(s) to obtain a 
 concise summary of the dag status.
 
 Example - 
         >>> python %prog --dag-file ihope.dag
        

"""

#
# =============================================================================
#
#                     Parse Command Line
#
# =============================================================================
#

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage )
  parser.add_option( "--dag-file", default='ihope.dag', metavar='DAGFILE', \
                   help="The dag file to parse")
  comand_line = sys.argv[1:]
  (options,args) = parser.parse_args()

  return options, sys.argv[1:]

# ============================================================================
# -- get command line arguments
options, args = parse_command_line()


##############################################################################
def get_status(filename, totals, debug=False):
  """
  parse a dag filem search for failed/done jobs
  @param filename: name of a dag file to parse
  @param totals: cumulitive number of jobs in each category
  @return status_int

  """
  tab = "    "
  dag_status = 0
  num = 0
  output = 0

  try:
    thisdag = open(filename + '.dagman.out', 'r')
  except:
    dag_status = -1

  if dag_status == 0:
    line = "dummy"
    while line:
      line = thisdag.readline()
      if line.find('Done')>=0:
        output = tab + line      
        line = thisdag.readline()
        output = output + tab + line
        #read the numbers
        line = thisdag.readline()
        #numbers = line.split(" ")
        num = line
        #number_of_undone_jobs = int(numbers[len(numbers)-1])
        output = output + tab + line
      if line.find('EXITING WITH STATUS 0')>=0:
        dag_status = 1

    # count the number of jobs
    if num != 0:  
      num=num.split()
      num = num[2:]
      for i in xrange( len( num ) ):
        totals[i] += int(num[i])

  elif dag_status == -1:
    try:
      thisdag = open(filename, 'r')
    except:
      dag_status = -2
    
    if dag_status == -1:
      line = "dummy"
      while line:
        line = thisdag.readline()
        if line.find('.sub')>=0:
          totals[7] += 1 

  if output != 0:
    print output[:-1]

  return dag_status


#
# =============================================================================
#
#                    Main
#
# =============================================================================
#
 

""" There are 2 possibilities. 
    First, the file provided is a dag file that ends in '.dag' that
    contains other dag calls such as in ihope.
    Second, it is a pure dag file which might not even contain a dag
    by itself, and in such case we search for the .dagman.out file.
"""

print 'Parsing ' + options.dag_file + '...'
sleep(1.5)

# Look for sub dags
dagFile = open(options.dag_file, 'r')
filenames = []
lines = dagFile.readlines()
for line in lines:
  if "SUBDAG" in line:
    strs = line.split()
    if "DIR" in strs:
      fname = strs[-1]
    else:
      fname = ''
    for string in strs:
      if ".dag" in string:
        filenames.append(os.path.join(fname,string))
  if "pegasus-plan" in line:
    strs = line.split()
    for string in strs:
      if ".dax" in string:
        string = string.replace('/dax/','/')
        string = string.replace('.dax','-0.dag')
        filenames.append(string)

dagFile.close()

if len(filenames)==0:
  print 'No dag files found in ' + options.dag_file
  print 'Assuming that ' + options.dag_file + ' is the dag file you want to parse'
  sub_dags = 0
else:
  sub_dags = len( filenames )
  print 'Found', str( sub_dags ), 'subdags'
done_dags = 0  

print 'Parsing the dag files for status...\n'
### we found some dag files hopefully

sleep(1.5)
totals = [0,0,0,0,0,0,0,0]

#parse all the dags one by one
for i in xrange(0,len(filenames), 1):
  filename = filenames[i]

  print '-------------------------------------------------------------------------'
  # How many jobs are completed/failed ?
  print "Parsing " + filenames[i]

  status = get_status(filename, totals)
  
  if status==1: 
    done_dags += 1 
    print "COMPLETE :)"
  elif status==0:
    print "incomplete :("
  elif status==-1:
   print "dag not yet started!"
  elif status==-2:
   print "dag pending and not yet planned!"


# Print totals
print "  --------------------------- "
print "      Done    =", totals[0]
#print " Pre     =", totals[1]
print "      Queued  =", totals[2]
#print " Post    =", totals[3]
print "      Ready   =", totals[4]
print "      Unready =", totals[5]
print "      Failed  =", totals[6]
print "  ----------------------------------- "
print "      Completed Jobs = " + str( totals[0] )
print "      Submitted Jobs = " + str( sum( totals[:7] ) )
print "      Total Jobs     = " + str( sum(totals[:]) )
if sub_dags != 0:
  print "      Sub-dags       =", str( done_dags ) + "/" + str( sub_dags )
print "  ----------------------------------- "

# Confirm if completed
dagFile = open(options.dag_file + '.dagman.out', 'r')
uber = 0
lines = dagFile.readlines()
for line in lines:
  if "EXITING WITH STATUS 0" in line:
    uber = 1
if uber == 1:
  print "  ihope status... COMPLETED!"
else:
  print "  ihope status... incomplete."
print ""
  




