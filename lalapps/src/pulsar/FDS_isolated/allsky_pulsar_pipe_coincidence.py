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

# read coincidence (polka) parameters
freq_window = cp.get('polka-params','freq_window')
alpha_window = cp.get('polka-params','alpha_window')
delta_window = cp.get('polka-params','delta_window')

# read chi-sq threshold params
Fth_chisq = cp.get('chisq-params','Fth_chisq')
chisq_points = cp.get('chisq-params','chisq_points')
ts1 = cp.get('chisq-params','ts1')
ts2 = cp.get('chisq-params','ts2')
sifo1 = cp.get('chisq-params','sifo1')
sifo2 = cp.get('chisq-params','sifo2')

# -------------------------------------------------------------------------------- #

# ------------------------- Set up working directory ----------------------------- # 

# name of local work directory
subdir=''.join([local_work_dir,'/','run.',str(job_id)])

# figure out name of results file that we're using to compute the upper limit
freq = float(start_freq) + float(job_id) * float(freq_band)

# names of Fstats files
fstats1_zipped=''.join(['Fstats-1-',str(freq),'.gz'])
fstats2_zipped=''.join(['Fstats-2-',str(freq),'.gz'])

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
mkdata=''.join([starting_dir,'/makefakedata_v2'])
fstatshape=''.join([starting_dir,'/lalapps_FstatShapeTest'])
makeinveto=''.join([starting_dir,'/makeInvetofile'])
polka=''.join([starting_dir,'/polka'])

#     ephemeris, timestamps and results file
earth=''.join([starting_dir,'/earth00-04.dat'])
sun=''.join([starting_dir,'/sun00-04.dat'])
times1=''.join([starting_dir,'/',ts1])
times2=''.join([starting_dir,'/',ts2])
fstats1_file=''.join([starting_dir,'/',fstats1_zipped])
fstats2_file=''.join([starting_dir,'/',fstats2_zipped])

# copy execulables and components to working sub-directory 
shutil.copy(cfstat,subdir)
shutil.copy(mkdata,subdir)
shutil.copy(fstatshape,subdir)
shutil.copy(makeinveto,subdir)
shutil.copy(polka,subdir)
shutil.copy(earth,subdir)
shutil.copy(sun,subdir)
shutil.copy(times1,subdir)
shutil.copy(times2,subdir)
shutil.copy(fstats1_file,subdir)
shutil.copy(fstats2_file,subdir)

# unzip fstats files:
unzip_fstats1=''.join(['gunzip ',fstats1_zipped])
os.system(unzip_fstats1)
FstatsFileName1=''.join(['Fstats-1-',str(freq)])


unzip_fstats2=''.join(['gunzip ',fstats2_zipped])
os.system(unzip_fstats2)
FstatsFileName2=''.join(['Fstats-2-',str(freq)])


# -------------------------------------------------------------------------------- #

# -------------------- run coincidence code (polka) ------------------------------ #

# name of output file for polka
polka_out=''.join(['polka_out',''.join(['-',str(freq)])])

polka_args = ' '.join(['./polka','-1',FstatsFileName1,'-2',FstatsFileName2,'-f',\
                       str(freq_window),'-a',str(alpha_window),'-d',\
                       str(delta_window),'-o',polka_out,\
                       '-3',FstatsFileName1,'-4',FstatsFileName2])

print 'running: ',polka_args
os.system(polka_args)

shutil.copy(polka_out,starting_dir)


# -------------------------------------------------------------------------------- #

# ----------------------------- chi-sq test -------------------------------------- #

# Need to read file and loop through coincident candidates; if F is sufficiently large
# need to perform a chi-sq test on one or both candidates

polka_file=open(polka_out,mode='r')

for line in polka_file:
  [sf1,sa1,sd1,sF1,sfa1,sf2,sa2,sd2,sF2,sfa2,sfa]=line.split(None,11)

  # initialise values of chi sq
  chisq1=-1.0
  chisq2=-1.0

  F1=float(sF1)
  if F1 > float(Fth_chisq):
    # 1) re-run ComputeFstat for a small band with -p (get pulsar parameters) option
    ifo=ifo1
    sifo=sifo1
    data=data1
    ts=ts1
    sa=sa1
    sd=sd1
    sf=sf1
    chi_fstart=float(sf)-float(chisq_points)*float(df) # starting frequency of search is chisq_points points to the left
    chi_freq_band=2.0*float(df)*float(chisq_points)     # band is 2*chisq_points 

          # define command line for run on first ifo
    cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',str(chi_freq_band),\
                          '-I',str(ifo),'-r',df,'-a',sa,'-d',sd,'-D',data,'-E . -y 00-04 -F',\
                          Fth,'-p'])
    print 'running: ',cfstat_args
    os.system(cfstat_args)
    
    # 2) run makeinvetofile makes the In.data file for makefakedata
    #    and checks that there's only one outlier
    makeinveto_args=' '.join(['./makeInvetofile','-f Fstats -p ParamMLE -o In.data -t',ts,\
                              '-l 1800.0 -n 20','-s',str(int(float(sf))-1),'-b 3.0'])
    print 'running: ',makeinveto_args
    os.system(makeinveto_args)

        # make signal directory to put sfts with fake data
    os.mkdir('signal')
    
    # 3) run makefakedata using Inveto.data file
    makefakedata_args=' '.join(['./makefakedata_v2','-i In.data -n signal/SFT -I',sifo,'-E .'])
    print 'running: ',makefakedata_args
    os.system(makefakedata_args)
    
    # 4) run ComputeFStat on fake signal data with signal only flag (-S)
    cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',str(chi_freq_band),\
                          '-I',str(ifo),'-r',str(df),'-a',sa,'-d',sd,'-D signal','-E . -y 00-04 -F 0.0',\
                          '-p -S'])
    print 'running: ',cfstat_args
    os.system(cfstat_args)    
        # remove signal directory
    os.system('rm -rf signal/')

    # 5) run Fstatshpetest and see if chisq test is passed
    fstatshape_args=' '.join(['./lalapps_FstatShapeTest -o FaFb00.001 -t FaFb01.001 > chisq.txt'])
    print 'running: ',fstatshape_args
    os.system(fstatshape_args)
  
    # 6) is chi sq test passed?
    chisq_file=open('chisq.txt',mode='r')

    for line in chisq_file:
      [crap,crap,crap,crap,schisq1]=line.split(None,5)
      print 'chi sq from data:',schisq1

    chisq1=float(schisq1)
    
    chisq_th=(4*(2*float(chisq_points)+1)-4)*(math.sqrt(F1)/10.0 +1)**2
    print 'theoretical chi sq:',chisq_th

    if chisq1 < chisq_th:
      chisq1_pass=1
    else: chisq1_pass=0

       #clean-up
    os.system('rm chisq.txt FaFb0* Fmax Fstats In.data ParamMLE')
    
  else: chisq1_pass=1

  F2=float(sF2)
  if F2 > float(Fth_chisq):
          # run chisq test on sub-candidate
    # 1) re-run ComputeFstat for a small band with -p (get pulsar parameters) option
    ifo=ifo2
    sifo=sifo2
    data=data2
    ts=ts2
    sa=sa2
    sd=sd2
    sf=sf2
    chi_fstart=float(sf)-float(chisq_points)*float(df) # starting frequency of search is chisq_points points to the left
    chi_freq_band=2.0*float(df)*float(chisq_points)     # band is 2*chisq_points 

          # define command line for run on first ifo
    cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',str(chi_freq_band),\
                          '-I',str(ifo),'-r',str(df),'-a',sa,'-d',sd,'-D',data,'-E . -y 00-04 -F',\
                          str(Fth),'-p'])
    print 'running: ',cfstat_args
    os.system(cfstat_args)
    
    # 2) run makeinvetofile makes the In.data file for makefakedata
    #    and checks that there's only one outlier
    makeinveto_args=' '.join(['./makeInvetofile','-f Fstats -p ParamMLE -o In.data -t',ts,\
                              '-l 1800.0 -n 20','-s',str(int(float(sf))-1),'-b 3.0'])
    print 'running: ',makeinveto_args
    os.system(makeinveto_args)

        # make signal directory to put sfts with fake data
    os.mkdir('signal')
    
    # 3) run makefakedata using Inveto.data file
    makefakedata_args=' '.join(['./makefakedata_v2','-i In.data -n signal/SFT -I',sifo,'-E .'])
    print 'running: ',makefakedata_args
    os.system(makefakedata_args)
    
    # 4) run ComputeFStat on fake signal data with signal only flag (-S)
    cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',str(chi_freq_band),\
                          '-I',str(ifo),'-r',str(df),'-a',sa,'-d',sd,'-D signal','-E . -y 00-04 -F 0.0',\
                          '-p -S'])
    print 'running: ',cfstat_args
    os.system(cfstat_args)    
        # remove signal directory
    os.system('rm -rf signal/')

    # 5) run Fstatshpetest and see if chisq test is passed
    fstatshape_args=' '.join(['./lalapps_FstatShapeTest -o FaFb00.001 -t FaFb01.001 > chisq.txt'])
    print 'running: ',fstatshape_args
    os.system(fstatshape_args)
  
    # 6) is chi sq test passed?
    chisq_file=open('chisq.txt',mode='r')

    for line in chisq_file:
      [crap,crap,crap,crap,schisq2]=line.split(None,5)
      print 'chi sq from data:',schisq2
    chisq2=float(schisq2)

    
    chisq_th=(4*(2*float(chisq_points)+1)-4)*(math.sqrt(F2)/10.0 +1)**2
    print 'theoretical chi sq:',chisq_th

    if chisq2 < chisq_th:
      chisq2_pass=1
    else: chisq2_pass=0

       #clean-up
    os.system('rm chisq.txt FaFb0* Fmax Fstats In.data ParamMLE')
    
  else: chisq2_pass=1
     

# Have both candidates passed the chi_sq veto (or not qualified to take it)?
  if chisq1_pass and chisq2_pass:
    # print coincident candidate to a file and exit    
    res_out=''.join(['results_out-',str(freq)])
    res_file=open(res_out,'w')
    results_string=' '.join([sf1,sa1,sd1,sF1,sfa1,str(chisq1),sf2,sa2,sd2,sF2,sfa2,str(chisq2),sfa])
    res_file.write(results_string)
    res_file.write("\n")
    res_file.close()
    shutil.copy(res_out, starting_dir)
    os.system(rm_subdir)
    sys.exit(0)
    
polka_file.close()

# -------------------------------------------------------------------------------- #

os.system(rm_subdir)

