#!/usr/bin/env python2
"""
findLoudest.py - given the Fstats-#-## files this py file determines
the loudest event every 20 Fstats files for each IFO
M.Alessandra Papa, May 2005, for the S2 analysis.

The starting_dir should have a subdir called Fstats_results/
that contains the input Fstats-#-## files.

The output is going to be a set of files called
Loudest.Fstats-#-## in the starting_dir
"""

# import standard modules to the python path
import sys, os, shutil, math,random
import getopt, re, string,popen2
import ConfigParser
sys.path.append('/usr/lib/python2')

pi=3.14159265358979323844

# Function usage
def usage():
  msg = """\
Usage: allsky_pulsar_pipe.py [options]

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

# read coincidence (polka) parameters
freq_window = cp.get('polka-params','freq_window')
coincidence_band = cp.get('polka-params','coincidence_band')

sNpolka = cp.get('mc-params','Npolka')
Npolka = int(sNpolka)


# -------------------------------------------------------------------------------- #

# ------------------------- Set up working directory ----------------------------- # 

# name of local work directory
subdir=''.join([local_work_dir,'/','findLoudestFstats.',str(job_id)])

#define this ere because it enters in confidence.data file name
freq=float(start_freq) + float(job_id) * float(coincidence_band) * Npolka

# remove local work directory in case it exists
rm_subdir=''.join(['rm -rf ',subdir])
os.system(rm_subdir)

# make local work directory
try: os.mkdir(local_work_dir)
except OSError, err:
  import errno
  print "Warning:", err
os.mkdir(subdir)
sys.stdout.flush()

# change to local working directory
os.chdir(subdir)


#=============

#=============


ifile=0
while ifile < Npolka:
  freq_dmp = float(start_freq)+float(job_id)*Npolka*float(coincidence_band)+ifile*float(coincidence_band)
  res_out=''.join([starting_dir,'/Fstats_results/Fstats-1-',str(freq)])
  if os.path.exists(res_out):
    sys.exit(0)
    print 'File ',res_out,' not found. Exiting ...'
  #endif the file exists
  res_out=''.join([starting_dir,'/Fstats_results/Fstats-2-',str(freq)])
  if os.path.exists(res_out):
    sys.exit(0)
    print 'File ',res_out,' not found. Exiting ...'
  #endif the file exists
  ifile=ifile+1
#end while


# -------------------------------------------------------------------------------- #

# -------------- extract and copy relevant data onto local dir  ------------------- # 
 
res1_out=''.join([starting_dir,'/Loudest.Fstats-1-',str(freq)])
res2_out=''.join([starting_dir,'/Loudest.Fstats-2-',str(freq)])


ifile=0
max1=0.0
max2=0.0
while ifile < Npolka:

  freq_dmp = float(start_freq)+float(job_id)*Npolka*float(coincidence_band)+ifile*float(coincidence_band)

########### first file
 
  gzpd_res_in=''.join(['Fstats-1-',str(freq_dmp),'.gz'])
  res_in=''.join(['Fstats-1-',str(freq_dmp)])
  gzpd_res_file=''.join([starting_dir,'/Fstats_results/',gzpd_res_in])
  if os.path.exists(gzpd_res_file)!= 1:
    print "could not find ",gzpd_res_files
    print "Proceeding to next band...."
    print "...."
    sys.exit(0)
  #endif os.path.exists

  #print subdir 
  shutil.copy(gzpd_res_file,subdir)
  #print 'Unzipping ',gzpd_res_in
  unzip_polka=''.join(['gunzip ',gzpd_res_in])
  os.system(unzip_polka)

  #sort_and_cat=''.join(['sort -nr -k 7,7 ',res_in,' | head -n 1 > dmp'])
  sort_and_cat=''.join(['sort -nr -k 7,7 ',res_in,' > dmp'])
  #print sort_and_cat
  os.system(sort_and_cat)
  #os.system('head -n1 dmp1 > dmp')

  print "======"
  print "ifile= ", ifile+1

  max_file=open('dmp',mode='r')
  #line_list=max_file.readlines()
  #length=len(line_list)
  #line=line_list[0]
  line=max_file.readline()
  [sf,sa,sd,sn,smu,sstd,sF]=line.split(None,7)
  print "1 ",sF
  if float(sF) > max1:
    max1=float(sF)
    fmax1=sf
    amax1=sa
    dmax1=sd
    nmax1=sn
    mumax1=smu
    stdmax1=sstd
    Fmax1=sF
  #endif
  max_file.close()
  

########### second file
 
 
  gzpd_res_in=''.join(['Fstats-2-',str(freq_dmp),'.gz'])
  res_in=''.join(['Fstats-2-',str(freq_dmp)])
  gzpd_res_file=''.join([starting_dir,'/Fstats_results/',gzpd_res_in])
  
  if os.path.exists(gzpd_res_file)!= 1:
    print "could not find ",gzpd_res_files
    print "Proceeding to next band...."
    print "...."
    sys.exit(0)
  #endif os.path.exists

  shutil.copy(gzpd_res_file,subdir)
  #print gzpd_res_file, subdir
  #print 'Unzipping ',gzpd_res_in
  unzip_polka=''.join(['gunzip ',gzpd_res_in])
  os.system(unzip_polka)
  
  sort_and_cat=''.join(['sort -nr -k 7,7 ',res_in,' > dmp'])
  os.system(sort_and_cat)
  #print sort_and_cat
  #print "====="

  max_file=open('dmp',mode='r')
  line=max_file.readline()
  [sf,sa,sd,sn,smu,sstd,sF]=line.split(None,7)
  print "2 ",sF
  if float(sF) > max2:
    max2=float(sF)
    fmax2=sf
    amax2=sa
    dmax2=sd
    nmax2=sn
    mumax2=smu
    stdmax2=sstd
    Fmax2=sF
  #endif
  max_file.close()

  ifile=ifile+1

  print " "

#end whileifile

print "Fmax 1 and 2 =",Fmax1,Fmax2

outfile=open(res1_out,'w')
print >>outfile,fmax1,amax1,dmax1,nmax1,mumax1,stdmax1,Fmax1
outfile.close()

outfile=open(res2_out,'w')
print >>outfile,fmax2,amax2,dmax2,nmax2,mumax2,stdmax2,Fmax2
outfile.close()


cleanup=''.join(['rm -rf ', subdir])
os.system(cleanup)
  
#print res1_out
#print res2_out
