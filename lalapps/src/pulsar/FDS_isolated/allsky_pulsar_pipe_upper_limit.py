#!/usr/bin/env python2.2
"""
pulsar_pipe.in - standalone pulsar pipeline driver script
X. Siemens
"""

# import standard modules to the python path
import sys, os, shutil, math,random
import getopt, re, string,popen2
import ConfigParser
sys.path.append('/usr/lib/python2.2')

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

# read mc threshold params
sNinj = cp.get('mc-params','Ninj')
Ninj=int(sNinj)
sh0 = cp.get('mc-params','h0')
h0=float(sh0)
sdh0 = cp.get('mc-params','dh0')
dh0=float(sdh0)
sc0 = cp.get('mc-params','c0')
c0=float(sc0)
stol = cp.get('mc-params','tol')
tol=float(stol)


# -------------------------------------------------------------------------------- #

# ------------------------- Set up working directory ----------------------------- # 

# name of local work directory
subdir=''.join([local_work_dir,'/','ul-run.',str(job_id)])

# figure out name of results file that we're using to compute the upper limit
freq = float(start_freq) + float(job_id) * float(freq_band)
res_in=''.join(['results_out-',str(freq)])

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
extract_data=''.join([starting_dir,'/extractSFTband'])
semi_F=''.join([starting_dir,'/lalapps_SemiAnalyticF'])
findSh=''.join([starting_dir,'/FindSh'])

#     ephemeris, timestamps and results file
earth=''.join([starting_dir,'/earth00-04.dat'])
sun=''.join([starting_dir,'/sun00-04.dat'])
times1=''.join([starting_dir,'/',ts1])
times2=''.join([starting_dir,'/',ts2])
res_file=''.join([starting_dir,'/',res_in])
fstats1_file=''.join([starting_dir,'/',fstats1_zipped])
fstats2_file=''.join([starting_dir,'/',fstats2_zipped])

# copy execulables and components to working sub-directory 
shutil.copy(cfstat,subdir)
shutil.copy(mkdata,subdir)
shutil.copy(fstatshape,subdir)
shutil.copy(makeinveto,subdir)
shutil.copy(polka,subdir)
shutil.copy(extract_data,subdir)
shutil.copy(semi_F,subdir)
shutil.copy(findSh,subdir)
shutil.copy(earth,subdir)
shutil.copy(sun,subdir)
shutil.copy(times1,subdir)
shutil.copy(times2,subdir)
shutil.copy(res_file,subdir)
shutil.copy(fstats1_file,subdir)
shutil.copy(fstats2_file,subdir)

# unzip fstats files:
unzip_fstats1=''.join(['gunzip ',fstats1_zipped])
os.system(unzip_fstats1)
fstats1=''.join(['Fstats-1-',str(freq)])

unzip_fstats2=''.join(['gunzip ',fstats2_zipped])
os.system(unzip_fstats2)
fstats2=''.join(['Fstats-2-',str(freq)])

# -------------------------------------------------------------------------------- #

# -------------- extract and copy relevant data onto local dir  ------------------- # 

# 1st ifo
os.mkdir('xdata1/')
data=data1

fmin=freq-2*float(freq_band)
band=5*float(freq_band)

# define command line first ifo
extract_args=' '.join(['./extractSFTband','-n xdata1/SFT','-d',data,'-N 20',
                       '-b',str(band),'-f',str(fmin)])

# extract data for first ifo
os.system(extract_args)

# 2nd ifo
os.mkdir('xdata2/')
data=data2

# define command line second ifo
extract_args=' '.join(['./extractSFTband','-n xdata2/SFT','-d',data,'-N 20',
                       '-b',str(band),'-f',str(fmin)])

# extract data for second ifo
os.system(extract_args)


fmin=freq-float(freq_band)
band=3*float(freq_band)

# -------------------------------------------------------------------------------- #

# --------------------------  Target false alarm  -------------------------------- #

results_file=open(res_in,mode='r')
for line in results_file:
  [dmp,dmp,dmp,dmp,dmp,dmp,dmp,dmp,dmp,dmp,dmp,dmp,sfa]=line.split(None,13)
results_file.close()
#target false alarm:
fa=float(sfa)

# -------------------------------------------------------------------------------- #

# --------------------------- MC injection loop  --------------------------------- # 
cont=1
confidence=0.0

inj_out=''.join(['Injection.data-',str(freq)])
inj_data_file=open(inj_out,mode='a')

while cont:
  i=0
  counter=0
  while i < Ninj:
    # Random inbjection parameters:
    cosi=random.uniform(-1.0,1.0)
    Ap=(1+cosi**2)/2 * h0
    Ac=cosi*h0
    psi=random.uniform(0.0,2*pi)
    phi01=random.uniform(0.0,2*pi)
    f0=random.uniform(freq,freq+float(freq_band))
    alpha=random.uniform(0.0,2*pi)
    cosdelta=random.uniform(-1.0,1.0)
    delta=math.acos(cosdelta)
    if delta > pi/2:
      delta = delta-pi

    phi02=2*pi*f0*(733803157-732489168)  # adjust phase in second injection by time elapsed
  
    # open and write in.data files
    # first ifo
    indata_file=open('In1.data',mode='w')
    print >>indata_file,'1800.0 ##TSFT'
    print >>indata_file,'20 ##NSFT'
    print >>indata_file,fmin,' ##fmin'
    print >>indata_file,band,' ##band'
    print >>indata_file,'-1.0 ##sigma_read_noise_flag'
    print >>indata_file,Ap,' ##Ap'
    print >>indata_file,Ac,' ##Ac'
    print >>indata_file,psi,' ##psi'
    print >>indata_file,phi01,' ##phi0'
    print >>indata_file,f0,' ##f0'
    print >>indata_file,delta,' ##latitude'
    print >>indata_file,alpha,' ##longitude'
    print >>indata_file,'0 ##SpinOrder'
    print >>indata_file,ts1,' ##TS'
    indata_file.close()

    #second ifo
    indata_file=open('In2.data',mode='w')
    print >>indata_file,'1800.0 ##TSFT'
    print >>indata_file,'20 ##NSFT'
    print >>indata_file,fmin,' ##fmin'
    print >>indata_file,band,' ##band'
    print >>indata_file,'-1.0 ##sigma_read_noise_flag'
    print >>indata_file,Ap,' ##Ap'
    print >>indata_file,Ac,' ##Ac'
    print >>indata_file,psi,' ##psi'
    print >>indata_file,phi02,' ##phi0'
    print >>indata_file,f0,' ##f0'
    print >>indata_file,delta,' ##latitude'
    print >>indata_file,alpha,' ##longitude'
    print >>indata_file,'0 ##SpinOrder'
    print >>indata_file,ts2,' ##TS'
    indata_file.close()

    # inject the fake signal into 1st ifo data
    os.mkdir('signal1')
    sifo=sifo1
    # run makefakedata using In1.data file
    makefakedata_args=' '.join(['./makefakedata_v2','-i In1.data -n signal1/SFT -I',sifo,'-D xdata1/','-E .'])
    os.system(makefakedata_args)

    # mismatched template variables
    alphaT=alpha+float(alpha_window)*random.uniform(-1.0,1.0)
    deltaT=delta+float(delta_window)*random.uniform(-1.0,1.0)
    fstart=f0-float(freq_window)

    # search for the signal 
    ifo=ifo1
    cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(fstart),'-b',str(2*float(freq_window)),\
                          '-I',str(ifo),'-r',str(df),'-a',str(alphaT),'-d',str(deltaT),
                          '-D signal1','-E . -y 00-04 -F 0.0','-o -1' ])
    os.system(cfstat_args)    

    fstats1_file=open('Fstats-1',mode='r')
    for line in fstats1_file:
      [sf1,dmp,dmp,dmp,dmp,dmp,sF1]=line.split(None,7)
    fstats1_file.close()
    
    # inject the fake signal into 2nd ifo data
    os.mkdir('signal2')
    sifo=sifo2
    # run makefakedata using In1.data file
    makefakedata_args=' '.join(['./makefakedata_v2','-i In2.data -n signal2/SFT -I',sifo,'-D xdata2/','-E .'])
    os.system(makefakedata_args)
    # search for the signal 
    ifo=ifo2
    cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(fstart),'-b',str(2*float(freq_window)),\
                          '-I',str(ifo),'-r',str(df),'-a',str(alphaT),'-d',str(deltaT),
                          '-D signal2','-E . -y 00-04 -F 0.0','-o -2'])
    os.system(cfstat_args)    
    fstats2_file=open('Fstats-2',mode='r')
    for line in fstats2_file:
      [sf2,dmp,dmp,dmp,dmp,dmp,sF2]=line.split(None,7)
    fstats2_file.close()

    # calculate theoretical values we expect from the noise and signal parameters
    findSh_1=' '.join(['./FindSh -D signal1 -b 0.03 -f',str(float(sf1)-0.015)])
    findSh_2=' '.join(['./FindSh -D signal2 -b 0.03 -f',str(float(sf2)-0.015)])
                                                                                                        
    myJob = popen2.Popen3(findSh_1)
    sSh1 = myJob.fromchild.readline()
    del myJob
    myJob = popen2.Popen3(findSh_2)
    sSh2 = myJob.fromchild.readline()
    del myJob
                                                                                                            
    semiF_args1=' '.join(['./lalapps_SemiAnalyticF -a',str(alpha),'-d',str(delta),'-Q',str(phi01/2.0),\
                          '-Y',str(psi),'-i',str(cosi),'-s',str(h0),'-N',str(float(sSh1)),'-T',ts1,\
                          '-t 1800.0 -n 20 -E . -D 1'])

    myJob = popen2.Popen3(semiF_args1)
    sF1th = myJob.fromchild.readline()
    del myJob
 
    semiF_args2=' '.join(['./lalapps_SemiAnalyticF -a',str(alpha),'-d',str(delta),'-Q',str(phi02/2.0),\
                          '-Y',str(psi),'-i',str(cosi),'-s',str(h0),'-N',str(float(sSh2)),'-T',ts2,\
                          '-t 1800.0 -n 20 -E . -D 2'])

    myJob = popen2.Popen3(semiF_args2)
    sF2th = myJob.fromchild.readline()
    del myJob

    # if the results are > the 20 and we are within our coincidence window in f0 then we run polka  
    if float(sF1) > 20.0 and float(sF2) > 20.0 and abs(float(sf1)-float(sf2)) < float(freq_window):  

      # run polka on output of search to compute significances
      polka_args = ' '.join(['./polka','-1 Fstats-1','-2 Fstats-2','-f',\
                           str(freq_window),'-a',str(alpha_window),'-d',\
                           str(delta_window),'-o polka_out',\
                           '-3',fstats1,'-4',fstats2])
      os.system(polka_args)
                                                                                                                                                                                           
      polka_file=open('polka_out',mode='r')
      for line in polka_file:
         [sf1,sa1,sd1,sF1,sfa1,sf2,sa2,sd2,sF2,sfa2,sfa]=line.split(None,11)
      polka_file.close()
                                                                                                                                                                                           
      # first check if false alarm is lower than for lodest candidate; if so proceed to chi-sq test
      if float(sfa) <= fa:

        F1=float(sF1)
        if F1 > float(Fth_chisq):
          # 1) re-run ComputeFstat for a small band with -p (get pulsar parameters) option
          ifo=ifo1
          sifo=sifo1
          data='signal1'
          ts=ts1
          sa=sa1
          sd=sd1
          sf=sf1
          # starting frequency of search is chisq_points points to the left
          chi_fstart=float(sf)-float(chisq_points)*float(df) 
          chi_freq_band=2.0*float(df)*float(chisq_points)     # band is 2*chisq_points 

          # define command line for run on first ifo
          cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',str(chi_freq_band),\
                                '-I',str(ifo),'-r',df,'-a',sa,'-d',sd,'-D',data,'-E . -y 00-04 -F',\
                                Fth,'-p'])
          os.system(cfstat_args)
    
          # 2) run makeinvetofile makes the In.data file for makefakedata
          #    and checks that there's only one outlier
          makeinveto_args=' '.join(['./makeInvetofile','-f Fstats -p ParamMLE -o In.data -t',ts,\
                                    '-l 1800.0 -n 20','-s',str(int(float(sf))-1),'-b 3.0'])
          os.system(makeinveto_args)

          # make signal directory to put sfts with fake data
          os.mkdir('signal')
    
          # 3) run makefakedata using Inveto.data file
          makefakedata_args=' '.join(['./makefakedata_v2','-i In.data -n signal/SFT -I',sifo,'-E .'])
          os.system(makefakedata_args)
    
          # 4) run ComputeFStat on fake signal data with signal only flag (-S)
          cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',str(chi_freq_band),\
                                '-I',str(ifo),'-r',str(df),'-a',sa,'-d',sd,'-D signal','-E . -y 00-04 -F 0.0',\
                                '-p -S'])
          os.system(cfstat_args)    
          # remove signal directory
          os.system('rm -rf signal/')

          # 5) run Fstatshpetest and see if chisq test is passed
          fstatshape_args=' '.join(['./lalapps_FstatShapeTest -o FaFb00.001 -t FaFb01.001 > chisq.txt'])
          os.system(fstatshape_args)
  
          # 6) is chi sq test passed?
          chisq_file=open('chisq.txt',mode='r')
          for line in chisq_file:
            [crap,crap,crap,crap,schisq1]=line.split(None,5)
          chisq_file.close()

          chisq1=float(schisq1)
          
          chisq_th=(4*(2*float(chisq_points)+1)-4)*(math.sqrt(F1)/10.0 +1)**2

          if chisq1 < chisq_th:
            chisq1_pass=1
          else:
            chisq1_pass=0
            print 'Candidate 1 did not survive chi sq test'

          # clean-up
          os.system('rm chisq.txt FaFb0* Fmax Fstats In.data ParamMLE')
    
        else: chisq1_pass=1

        F2=float(sF2)
        if F2 > float(Fth_chisq):
          # run chisq test on sub-candidate
          # 1) re-run ComputeFstat for a small band with -p (get pulsar parameters) option
          ifo=ifo2
          sifo=sifo2
          data='signal2'
          ts=ts2
          sa=sa2
          sd=sd2
          sf=sf2
          # starting frequency of search is chisq_points points to the left
          chi_fstart=float(sf)-float(chisq_points)*float(df) 
          chi_freq_band=2.0*float(df)*float(chisq_points)     # band is 2*chisq_points 

          # define command line for run on first ifo
          cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',str(chi_freq_band),\
                                '-I',str(ifo),'-r',str(df),'-a',sa,'-d',sd,'-D',data,'-E . -y 00-04 -F',\
                                str(Fth),'-p'])
          os.system(cfstat_args)
      
          # 2) run makeinvetofile makes the In.data file for makefakedata
          #    and checks that there's only one outlier
          makeinveto_args=' '.join(['./makeInvetofile','-f Fstats -p ParamMLE -o In.data -t',ts,\
                                    '-l 1800.0 -n 20','-s',str(int(float(sf))-1),'-b 3.0'])
          os.system(makeinveto_args)

          # make signal directory to put sfts with fake data
          os.mkdir('signal')
       
          # 3) run makefakedata using Inveto.data file
          makefakedata_args=' '.join(['./makefakedata_v2','-i In.data -n signal/SFT -I',sifo,'-E .'])
          os.system(makefakedata_args)
    
          # 4) run ComputeFStat on fake signal data with signal only flag (-S)
          cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',str(chi_freq_band),\
                                '-I',str(ifo),'-r',str(df),'-a',sa,'-d',sd,'-D signal','-E . -y 00-04 -F 0.0',\
                                '-p -S'])
          os.system(cfstat_args)    
          # remove signal directory
          os.system('rm -rf signal/')

          # 5) run Fstatshpetest and see if chisq test is passed
          fstatshape_args=' '.join(['./lalapps_FstatShapeTest -o FaFb00.001 -t FaFb01.001 > chisq.txt'])
          os.system(fstatshape_args)
  
          # 6) is chi sq test passed?
          chisq_file=open('chisq.txt',mode='r')
          for line in chisq_file:
            [crap,crap,crap,crap,schisq2]=line.split(None,5)
          chisq_file.close()
        
          chisq2=float(schisq2)
            
          chisq_th=(4*(2*float(chisq_points)+1)-4)*(math.sqrt(F2)/10.0 +1)**2

          if chisq2 < chisq_th:
            chisq2_pass=1
          else:
            chisq2_pass=0
            print 'Candidate 2 did not survive chi sq test'
            
          # clean-up
          os.system('rm chisq.txt FaFb0* Fmax Fstats In.data ParamMLE')
       
        else: chisq2_pass=1
     

        # Have both candidates passed the chi_sq veto (or not qualified to take it)?
        if chisq1_pass and chisq2_pass:
          counter=counter+1

        print sfa1,sfa2,float(sfa),fa,h0,float(sF1),2*float(sF1th),\
              float(sF2),2*float(sF2th),abs(float(sf1)-float(sf2)),\
              float(counter),float(i+1),float(counter)/float(i+1)
        print >>inj_data_file,sfa1,sfa2,float(sfa),fa,h0,float(sF1),2*float(sF1th),\
              float(sF2),2*float(sF2th),abs(float(sf1)-float(sf2)),\
              float(counter),float(i+1),float(counter)/float(i+1)

    print i+1,counter,float(counter)/float(i+1)
    os.system('rm -rf signal1 signal2')            
    i=i+1

  confidenceOLD=confidence  
  confidence=float(counter)/float(Ninj)
  # first ifo
  res_out=''.join(['Confidence.data-',str(freq)])
  outCdata_file=open(res_out,mode='a')
  print >>outCdata_file,h0,confidence
  outCdata_file.close()
  shutil.copy(res_out,starting_dir)
  
  if abs(confidence-c0) < tol:
    cont=0
  if confidence > c0:
    if confidenceOLD < c0:
      dh0=dh0/2
    h0=h0-dh0
  if confidence < c0:
    if confidenceOLD > c0:
      dh0=dh0/2
    h0=h0+dh0

inj_data_file.close()
shutil.copy(inj_out,starting_dir)


# -------------------------------------------------------------------------------- #
