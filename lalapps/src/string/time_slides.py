#!/usr/bin/env python
"""
stand alone time-slides script 
X. Siemens
"""

# import standard modules to the python path
import sys, os, shutil, math,random
import getopt, re, string,popen2
import ConfigParser
sys.path.append('/usr/lib/python')

# Function usage
def usage():
  msg = """\
Usage: time_slides.py [options]
  -h, --help               display this message
  -p  --params-file        Search parameter configuration file         
"""
  print >> sys.stderr, msg


# ------------------- parse the command line options -------------------- #

# initialise command line argument variables
params_file = None

shortop = "hp:"
longop = [
   "help",
   "params-file=",
   ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)
  
for o, a in opts:
  if o in ("-p", "--params-file"):
    params_file = a      
  elif o in ("-h", "--help"):
    usage()
    sys.exit(0)
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

if not params_file:
  print >> sys.stderr, "No parameter file specified."
  print >> sys.stderr, "Use --params-file to specify it."
  sys.exit(1)


# -------------------------------------------------------------------------------- #

# ------------------------- read configuration file ------------------------------ # 

cp = ConfigParser.ConfigParser()
cp.read(params_file)

# read time-slide parameters
slide_time = int(cp.get('pipeline','slide-time'))
slide_time_ns =  int(cp.get('pipeline','slide-time-ns'))
nslides  =  int(cp.get('pipeline','number-slides'))

# get chunk lengths from the values in the ini file
pad = int(cp.get('string', 'pad'))
length = int(cp.get('pipeline', 'segment-length'))+2*pad
short_seg_lenth = int(cp.get('string', 'short-segment-duration'))

# get coincidence time windows
dtHL = int(cp.get('pipeline','dtHL'))
dtHH = int(cp.get('pipeline','dtHH'))

# ----------- Check existence of merged trigger files and segment file  ---------- #

if not os.path.exists('H1triggers.xml'):
  print >> sys.stderr, "Can't find merged H1 trigger file (H1triggers.xml)"
  sys.exit(1)
if not os.path.exists('H2triggers.xml'):
  print >> sys.stderr, "Can't find merged H2 trigger file (H2triggers.xml)"
  sys.exit(1)
if not os.path.exists('L1triggers.xml'):
  print >> sys.stderr, "Can't find merged L1 trigger file (L1triggers.xml)"
  sys.exit(1)
if not os.path.exists('seg.txt'):
  print >> sys.stderr, "Can't find segment file (seg.txt)"
  sys.exit(1)
  
# Set up slides directory
# name of local work directory
os.system('rm -rf ./slides/')
os.mkdir('./slides/')
i=1
while i< nslides+1:
  os.mkdir('./slides/slideP'+str(i))
  os.mkdir('./slides/slideM'+str(i))
  i=i+1

# ------------------------- Foreground computation of H1 survivors -------------- # 

#H1 survivors in coincidence with H2:
print 'RUNNING: ~/lscsoft/lalapps/src/power/lalapps_burca --output-dir slides/ --ifo-a H1triggers.xml --ifo-b \
H2triggers.xml --ignore-playground --ignore-tfcomparison --dt '+str(dtHH)+' --no-repeats --amplitude-cut 1'
os.system('~/lscsoft/lalapps/src/power/lalapps_burca --output-dir slides/ --ifo-a H1triggers.xml --ifo-b \
H2triggers.xml --ignore-playground --ignore-tfcomparison --dt '+str(dtHH)+' --no-repeats --amplitude-cut 1')

#H1-H2 survivors in coincidence with L1:
print 'RUNNING: ~/lscsoft/lalapps/src/power/lalapps_burca --output-dir slides/ --ifo-a slides/H1-BURCA_H1H2_P_0_--1-0.xml \
--ifo-b L1triggers.xml --ignore-playground --ignore-tfcomparison --dt '+str(dtHL)+' --no-repeats'
os.system('~/lscsoft/lalapps/src/power/lalapps_burca --output-dir slides/ --ifo-a slides/H1-BURCA_H1H2_P_0_--1-0.xml \
--ifo-b L1triggers.xml --ignore-playground --ignore-tfcomparison --dt '+str(dtHL)+' --no-repeats')

# -------------------------------------------------------------------------------- #

# ------------------------- Background Computation ------------------------------- # 

seg_file=open('seg.txt',mode='r')
seg_list=seg_file.readlines()

for seg in seg_list:
  [crap,Tstart,Tend,duration]=seg.split(None,4)
  if int(duration) > length:
    Ti=int(Tstart)+pad+short_seg_lenth/4
    Tf=int(Tend)-pad-short_seg_lenth/4
    print 'RUNNING: ~/lscsoft/lalapps/src/power/lalapps_burca --start-time '+ str(Ti)+' --stop-time '+str(Tf)+\
    ' --output-dir slides/ --ifo-a slides/H1-BURCA_H1H2_P_0_--1-0.xml --ifo-b L1triggers.xml --ignore-playground \
    --ignore-tfcomparison --dt '+str(dtHL)+' --no-repeats \
    --slide-time '+str(slide_time)+' --slide-time-ns '+str(slide_time_ns)+' --number-slides '+str(nslides) 
    os.system('~/lscsoft/lalapps/src/power/lalapps_burca --start-time '+ str(Ti)+' --stop-time '+str(Tf)+\
    ' --output-dir slides/ --ifo-a slides/H1-BURCA_H1H2_P_0_--1-0.xml --ifo-b L1triggers.xml --ignore-playground \
    --ignore-tfcomparison --dt '+str(dtHL)+' --no-repeats \
    --slide-time '+str(slide_time)+' --slide-time-ns '+str(slide_time_ns)+' --number-slides '+str(nslides))
    i=1
    while i< nslides+1:
      tSlide=1e9*slide_time+slide_time_ns
      dirnameM='./slides/slideM'+str(i)
      dirnameP='./slides/slideP'+str(i)
      os.system('mv ./slides/H1-BURCA_H1L1_M*'+str(int(i*tSlide))+'* '+dirnameM)
      os.system('mv ./slides/H1-BURCA_H1L1_P*'+str(int(i*tSlide))+'* '+dirnameP)
      i=i+1

# -------------------------------------------------------------------------------- #
