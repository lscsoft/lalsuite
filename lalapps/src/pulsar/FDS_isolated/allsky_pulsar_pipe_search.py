#!/usr/bin/env python2.2
"""
pulsar_pipe.in - standalone pulsar pipeline driver script
X. Siemens
"""

# import standard modules to the python path
import sys, os, shutil, math
import getopt, re, string
import ConfigParser
sys.path.append('/usr/lib/python2.2')

# Function usage
def usage():
  msg = """\
Usage: allsky_pulsar_pipe.in [options]

  -h, --help               display this message
  -j, --job-id             job ID #  (used to compute starting frequency of this particular job)
  -S  --starting-dir       Starting directory (location of ComputeFStatistic, ephemeris files etc...)
  -W  --local-work-dir     Local working directory (on the nodes)
  -p  --params-file        Search parameter configuration file         
"""
  print >> sys.stderr, msg


# ------------------- parse the command line options -------------------- #

# initialise command line argument variables
job_id = -1.0
starting_dir = None
local_work_dir = None
params_file = None

shortop = "hj:S:W:p:"
longop = [
   "help",
   "job-id=",
   "starting-dir=",
   "local-work-dir=",
   "params-file=",
   ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)
  
for o, a in opts:
  if o in ("-j", "--job-id"):
    job_id = int(a)
  elif o in ("-S", "--starting-dir"):
    starting_dir = a      
  elif o in ("-W", "--local-work-dir"):
    local_work_dir = a      
  elif o in ("-p", "--params-file"):
    params_file = a      
  elif o in ("-h", "--help"):
    usage()
    sys.exit(0)
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

if job_id == -1.0:
  print >> sys.stderr, "No job number specified."
  print >> sys.stderr, "Use --job-id to specify it."
  sys.exit(1)
if not starting_dir:
  print >> sys.stderr, "No starting directory specified."
  print >> sys.stderr, "Use --starting-dir to specify it."
  sys.exit(1)
if not local_work_dir:
  print >> sys.stderr, "No local working directory specified."
  print >> sys.stderr, "Use --local-work-dir to specify it."
  sys.exit(1)
if not params_file:
  print >> sys.stderr, "No search parameter file specified."
  print >> sys.stderr, "Use --params-file to specify it."
  sys.exit(1)

# -------------------------------------------------------------------------------- #

# ------------------------- read configuration file ------------------------------ # 

cp = ConfigParser.ConfigParser()
cp.read(params_file)

# read F-statistic search parameters
start_freq = cp.get('fstat-params','start_freq')
freq_band = cp.get('fstat-params','freq_band')
df =  cp.get('fstat-params','df')
ifo1 = cp.get('fstat-params','ifo1')
ifo2 = cp.get('fstat-params','ifo2')
a_search = cp.get('fstat-params','a_search')
d_search = cp.get('fstat-params','d_search')
data1 = cp.get('fstat-params','data1')
data2 = cp.get('fstat-params','data2')
Fth = cp.get('fstat-params','Fth')
wings = cp.get('fstat-params','wings')

# -------------------------------------------------------------------------------- #

# ------------------------- Set up working directory ----------------------------- # 

# name of local work directory
subdir=''.join([local_work_dir,'/','run.',str(job_id)])

# remove local work directory in case it exists
rm_subdir=''.join(['rm -rf ',subdir])
os.system(rm_subdir)

# make local work directory
try: os.mkdir(local_work_dir)
except OSError, err:
  import errno
  print "Warning:", err
os.mkdir(subdir)

# change to local working directory
os.chdir(subdir)

# define paths (as strings) to necessary executables and components
#     executables
cfstat=''.join([starting_dir,'/lalapps_ComputeFStatistic'])

#     ephemeris and timestamps
earth=''.join([starting_dir,'/earth00-04.dat'])
sun=''.join([starting_dir,'/sun00-04.dat'])

# copy execulables and components to working sub-directory 
shutil.copy(cfstat,subdir)
shutil.copy(earth,subdir)
shutil.copy(sun,subdir)

# -------------------------------------------------------------------------------- #

# ------------------------- run the F-statistic code ----------------------------- # 

# starting frequency for this job
freq = float(start_freq) + float(job_id) * float(freq_band) - float(wings)

# first ifo
ifo=ifo1
data=data1

# define command line for run on first ifo
cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(freq),'-b',str(float(freq_band)+2*float(wings)),\
                      '-I',ifo,'-r',df,a_search,d_search,\
                      '-D',data,'-E . -y 00-04 -F',Fth,'-o',''.join(['-',ifo,'-',str(freq)])])

#execute ComputeFStatistic code on first ifo
print 'running: ',cfstat_args
os.system(cfstat_args)
FstatsFileName1=''.join(['Fstats','-',ifo,'-',str(freq)])

#second ifo
ifo=ifo2
data=data2

# define command line for run on second ifo
cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(freq),'-b',freq_band,\
                      '-I',ifo,'-r',df,a_search,d_search,\
                      '-D',data,'-E . -y 00-04 -F',Fth,'-o',''.join(['-',ifo,'-',str(freq)])])
#execute ComputeFStatistic code on second ifo
print 'running: ',cfstat_args
os.system(cfstat_args)
FstatsFileName2=''.join(['Fstats','-',ifo,'-',str(freq)])

# gzip fstats files and copy them to starting dir
zip_fstats1=''.join(['gzip ',FstatsFileName1])
os.system(zip_fstats1)
zipped_fstats1=''.join([FstatsFileName1,'.gz'])
                                                                                                                                                                                                    
shutil.copy(zipped_fstats1,starting_dir)
                                                                                                                                                                                                    
zip_fstats2=''.join(['gzip ',FstatsFileName2])
os.system(zip_fstats2)
zipped_fstats2=''.join([FstatsFileName2,'.gz'])
                                                                                                                                                                                                    
shutil.copy(zipped_fstats2,starting_dir)
                                                                                                                                                                                                    
# -------------------------------------------------------------------------------- #

os.system(rm_subdir)


