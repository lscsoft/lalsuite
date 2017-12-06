#!/usr/bin/python
"""
FstatFollowUp.py - Uses the coherent ComputeFStatistic_v2 code to follow up on a given parameter space.
"""

#Imports necessary modules
import os,sys,getopt

for this_IFO in ['H1','H2','L1']:
  listFile = open(''.join(['TotalListOnNodes_',this_IFO]),'w')
  for node in range(1,501):
    if (node != 277):
      print str(node)
      os.sys.stdout.flush()
      try:
        dir_string = ''.join(['/data/node',str(node),'/frames/S5/sfts/'])
        if this_IFO in ['H1','H2']:
          dir_string = ''.join([dir_string,'LHO/'])
        elif this_IFO in ['L1']:
          dir_string = ''.join([dir_string,'LLO/'])
        list_of_dir = os.listdir(dir_string)
        for directory in list_of_dir:
          if 'C03hoft' in directory:
            new_dir_string = ''.join([dir_string,directory,'/'])
            file_names = os.listdir(new_dir_string)
            for name in file_names:
              if this_IFO in name:
                listFile.write(''.join([new_dir_string,name,'\n']))
      except:
        print "Failed on " + str(node)
        os.sys.stdout.flush()
  listFile.close()
