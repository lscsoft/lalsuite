#!/usr/bin/python
"""
FstatFollowUpSubmit.py - Submits the actual job in the scheduler universe
"""

#Imports necessary modules
import os,sys,getopt,time

################################
#Prints the usage of this script
def usage():
  msg = """\
Usage: FstatFollowUpSubmit.py [options]

  -h, --help               display this message
  -d, --directory      Directory to run dag from

  -l, --log            Log directory
  -f, --file          Dag file to run
"""
  print >> sys.stderr, msg

################################################
#This section reads in the command line arguments
shortop = "hl:d:f:"
longop = [
  "help",
  "log=",
  "directory=",
  "file=",
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
  elif o in ("-d","--directory"):
    start_directory = str(a)
  elif o in ("-f","--file"):
    dag_file = str(a)
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

os.chdir(start_directory)
print "Full name found " + dag_file + " is " + str(os.path.exists(dag_file))
print "Relative name found is " + str(os.path.exists(os.path.basename(dag_file)))


time.sleep(50)
loops = 0
while loops < 20:
  if os.path.exists(dag_file):
    print "Full name found " + dag_file + " is " + str(os.path.exists(dag_file))
    print "Relative name found is " + str(os.path.exists(os.path.basename(dag_file)))
    print "I found it"
    command_line = ' '.join(['condor_submit_dag -outfile_dir',log_dir,dag_file])
    print command_line
    os.system(command_line)
    sys.exit(0)                    
  else:
    print str(loops)
    time.sleep(5)
    loops = loops+1
    
    
sys.exit(-1)

