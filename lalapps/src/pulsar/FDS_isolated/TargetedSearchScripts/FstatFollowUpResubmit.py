#!/usr/bin/python
"""
FstatFollowUpSubmit.py - Submits the actual job in the scheduler universe
"""

#Imports necessary modules
import os,sys,getopt,time

print sys.version

################################
#Prints the usage of this script
def usage():
  msg = """\
Usage: FstatFollowUp.py [options]

  -h, --help               display this message

  -l, --log                Log directory
  -f, --file               File to check exists and determines whether to run
  -d, --directory          Directory to call condor_dag_submit from
"""
  print >> sys.stderr, msg

################################################
#This section reads in the command line arguments
shortop = "hl:f:d:"
longop = [
  "help",
  "log=",
  "file=",
  "directory=",
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:],shortop,longop)
except getopt.GetoptError:
  print getopt.GetoptError
  usage()
  sys.exit(1)

for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-l","--log"):
    log_dir = str(a)
  elif o in ("-f","--file"):
    dag_file = str(a)
  elif o in ("-d","--directory"):
    start_directory = str(a)
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

time.sleep(50)

os.chdir(start_directory)

if os.path.exists(dag_file):
  command_line = ' '.join(['condor_submit_dag -outfile_dir',log_dir,dag_file])
  os.system(command_line)

