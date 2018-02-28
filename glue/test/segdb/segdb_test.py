#!/usr/bin/env python



"""This script test ligolw_segment_query and ligolw_segments_from_cats:

For ligolw_segment_query, the test runs against E13 data at the time of writing.
1. run test query "ligolw_segments_from_cats --segment-url=https://segdb.ligo.caltech.edu --gps-start-time 924821634 --gps-end-time 924828992 --veto-file=/H1H2-CBC_E13_ONLINE-923682800-2419200.xml --separate-categories"
2. get segment start_time, end_time out from the result test_ligolw_segment_query.xml and put them in a temp time file "segScript"
3. diff segScript against the validated correct results in the "correct_ligolw_segment_query_results.txt" 

For ligolw_segments_from_cats, the test runs against E13 data at the time of writing
1. run test query "ligolw_segments_from_cats --segment-url=https://segdb.ligo.caltech.edu --gps-start-time 924821634 --gps-end-time 924828992 --veto-file=/H1H2-CBC_E13_ONLINE-923682800-2419200.xml --separate-categories". This command returns 8 xml files.
2. loop in the 12 xml files to get segment start_time and end_time and put them in a temp time file, for example, result_H1CAT1
3. in the loop, diff the temp time file against its correspond validated results in, for example, H1CAT1 
"""


from glue import segments
import commands
import sys
import os


#--------------------------------------------------------------------------------
#    Test ligolw_segment_query without explicit versions
#--------------------------------------------------------------------------------
print "Testing ligolw_segment_query against E13 data (without versions)..."

# run the testing ligolw_segment_query command and generate the result xml file
com = "ligolw_segment_query --segment-url=https://segdb.ligo.caltech.edu --gps-start-time 924821632 --gps-end-time 924921632 --include-segments H1:DMT-SCIENCE --exclude-segments H1:DMT-BADGAMMA --query-segments | ligolw_print -t segment -c start_time -c end_time -d ' ' > segScript"
a = commands.getstatusoutput(com)
if a[0] == 0:
  pass
else:
  print "Error executing command to generate result xml file"
  sys.exit(1)


# diff result file from ligolw_segment_query and from database
com = 'diff correct_ligolw_segment_query_results.txt segScript'
a = commands.getstatusoutput(com)
if a[0] == 0:
  print "Test pass"
  print
else:
  print "Test fail"
  print a[1]

os.remove('segScript')

#--------------------------------------------------------------------------------
#    Test ligolw_segment_query with explicit versions
#--------------------------------------------------------------------------------
print "Testing ligolw_segment_query against E13 data (with versions)..."

# run the testing ligolw_segment_query command and generate the result xml file
com = "ligolw_segment_query --segment-url=https://segdb.ligo.caltech.edu --gps-start-time 924821632 --gps-end-time 924921632 --include-segments H1:DMT-SCIENCE:1 --exclude-segments H1:DMT-BADGAMMA:1 --query-segments | ligolw_print -t segment -c start_time -c end_time -d ' ' > segScript"
a = commands.getstatusoutput(com)
if a[0] == 0:
  pass
else:
  print "Error executing command to generate result xml file"
  sys.exit(1)


# diff result file from ligolw_segment_query and from database
com = 'diff correct_ligolw_segment_query_results.txt segScript'
a = commands.getstatusoutput(com)
if a[0] == 0:
  print "Test pass"
  print
else:
  print "Test fail"
  print a[1]

os.remove('segScript')

#---------------------------------------------------------------------------------
#     Test ligolw_segments_from_cats
#---------------------------------------------------------------------------------
print
print "Testing ligolw_segments_from_cats against E13 data..."
print "        It may take a while ..."

# run ligolw_segments_from_cats and get 12 result files back
com = "ligolw_segments_from_cats --segment-url=https://segdb.ligo.caltech.edu --gps-start-time 924821634 --gps-end-time 924828992 --veto-file=H1H2-CBC_E13_ONLINE-923682800-2419200.xml --separate-categories"
a = commands.getstatusoutput(com)
if a[0] == 0:
  pass
else:
  print "Error executing ligolw_segments_from_cats command"
  sys.exit(1)

ret = 0
for i in ['H1', 'H2']:
  for c in [1,2,3,4]: # loop in categories
    # get the segment start and end time from the result xml file, and put them in a temp time file
    result_file_name = i + '-VETOTIME_CAT' + str(c) + '-924821634-7358.xml' 
    com = "cat " + result_file_name + " | ligolw_print -t segment -c start_time -c end_time -d ' ' > " + "result_" + i + "CAT" + str(c)
    a = commands.getstatusoutput(com)
    if a[0] != 0:
      print "Error execute command to get segment start and end time from result xml file"
      sys.exit(0)
  
    # diff result file agaisnt the correct results
    com = 'diff ' + i + 'CAT' + str(c) + " result_" + i + "CAT" + str(c) 
    a = commands.getstatusoutput(com)
    if a[0]!=0:
      print "Error diff time file %s from %s" % ("result_" + i + "CAT" + str(c), result_file_name)
      ret = 1 
    # remove the temp result file and the time file
    os.remove("result_" + i + "CAT" + str(c))
    os.remove(result_file_name)

if ret == 0:
  print "test pass"

sys.exit(0)
