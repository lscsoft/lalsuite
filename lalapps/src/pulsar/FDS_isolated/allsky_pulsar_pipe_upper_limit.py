#!/usr/bin/env python2
"""
pulsar_pipe.in - standalone pulsar pipeline driver script
X. Siemens
modified by M.Alessandra Papa starting in Feb 05
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
df =  cp.get('fstat-params','df')
ifo1 = cp.get('fstat-params','ifo1')
ifo2 = cp.get('fstat-params','ifo2')
a_search = cp.get('fstat-params','a_search')
d_search = cp.get('fstat-params','d_search')
data1 = cp.get('fstat-params','data1')
data2 = cp.get('fstat-params','data2')
sFth = cp.get('fstat-params','Fth')

# read coincidence (polka) parameters
freq_window = cp.get('polka-params','freq_window')
alpha_window = cp.get('polka-params','alpha_window')
delta_window = cp.get('polka-params','delta_window')
coincidence_band = cp.get('polka-params','coincidence_band')

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
sgps_start1 = cp.get('mc-params','gps_start1')
sgps_start2 = cp.get('mc-params','gps_start2')
gps_start1 = int(sgps_start1)
gps_start2 = int(sgps_start2)
sNpolka = cp.get('mc-params','Npolka')
Npolka = int(sNpolka)

# -------------------------------------------------------------------------------- #

# ------------------------- Set up working directory ----------------------------- # 

# name of local work directory
subdir=''.join([local_work_dir,'/','ul-run.',str(job_id)])

# check that we've got all the polka files that we need

ipolka=0
while ipolka < Npolka:
 freq_dmp = float(start_freq) + (float(job_id)+ipolka) * float(coincidence_band)
 gzpd_res_in=''.join(['polka_out-',str(freq_dmp),'.gz'])
 res_in=''.join(['polka_out-',str(freq_dmp)])
 gzpd_res_file=''.join([starting_dir,'/polka_results/',gzpd_res_in])
 ipolka=ipolka+1
 if os.path.exists(gzpd_res_file)!= 1:
  print "Not all necessary polka files are present"
  print "could not find ",gzpd_res_in
  print "Proceeding to next band...."
  print "...."
  sys.exit(0)

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

# define paths (as strings) to necessary executables and components
#     executables
cfstat=''.join([starting_dir,'/lalapps_ComputeFStatistic'])
mkdata=''.join([starting_dir,'/makefakedata_v2'])
fstatshape=''.join([starting_dir,'/lalapps_FstatShapeTestLAL'])
makeinveto=''.join([starting_dir,'/lalapps_makeInvetofile'])
polka=''.join([starting_dir,'/lalapps_polka'])
extract_data=''.join([starting_dir,'/lalapps_extractSFTband'])
semi_F=''.join([starting_dir,'/lalapps_SemiAnalyticF'])
findSh=''.join([starting_dir,'/lalapps_FindSh'])

#     ephemeris, timestamps and results file
earth=''.join([starting_dir,'/earth00-04.dat'])
sun=''.join([starting_dir,'/sun00-04.dat'])
times1=''.join([starting_dir,'/',ts1])
times2=''.join([starting_dir,'/',ts2])

# copy execulables and components to working sub-directory 
shutil.copy(cfstat,subdir)
shutil.copy(mkdata,subdir)
shutil.copy(fstatshape,subdir)
shutil.copy(makeinveto,subdir)
shutil.copy(extract_data,subdir)
shutil.copy(semi_F,subdir)
shutil.copy(findSh,subdir)
shutil.copy(earth,subdir)
shutil.copy(sun,subdir)
shutil.copy(times1,subdir)
shutil.copy(times2,subdir)



# copy over and unzip input polka files
#############################################
ipolka=0
while ipolka < Npolka:
 freq_dmp = float(start_freq) + (float(job_id)+ipolka) * float(coincidence_band)
 gzpd_res_in=''.join(['polka_out-',str(freq_dmp),'.gz'])
 res_in=''.join(['polka_out-',str(freq_dmp)])
 gzpd_res_file=''.join([starting_dir,'/polka_results/',gzpd_res_in])
 gzpd_res_file=''.join([starting_dir,'/polka_results/',gzpd_res_in])
 shutil.copy(gzpd_res_file,subdir)
 print 'Unzipping ',gzpd_res_in
 unzip_polka=''.join(['gunzip ',gzpd_res_in])
 os.system(unzip_polka)
 ipolka=ipolka+1
#############################################


# -------------------------------------------------------------------------------- #

# -------------- extract and copy relevant data onto local dir  ------------------- # 

# 1st ifo
os.mkdir('xdata1/')
data=data1

freq=float(start_freq) + float(job_id) * float(coincidence_band) * Npolka
fmin=freq-2.0*float(coincidence_band)
band=float(coincidence_band) * (Npolka + 4.0)

print 'freq,fmin,band=',freq,fmin,band

# define command line first ifo
#print e
extract_args=' '.join(['./lalapps_extractSFTband','-n xdata1/SFT','-d',data,'-N 20',
                       '-b',str(band),'-f',str(fmin)])

print extract_args

# extract data for first ifo
#print extract_args
os.system(extract_args)

# 2nd ifo
os.mkdir('xdata2/')
data=data2

# define command line second ifo
extract_args=' '.join(['./lalapps_extractSFTband','-n xdata2/SFT','-d',data,'-N 20',
                       '-b',str(band),'-f',str(fmin)])

# extract data for second ifo
print extract_args
os.system(extract_args)


# -------------------------------------------------------------------------------- #

# --------------------------  Target false alarm  -------------------------------- #

ipolka=0
fa=100.0
while ipolka < Npolka:
 freq_dmp = float(start_freq) + (float(job_id)+ipolka) * float(coincidence_band)
 res_in=''.join(['polka_out-',str(freq_dmp)])
 results_file=open(res_in,mode='r')
 line=results_file.readline()
 words=line.split()
 sfa=words[10]
 results_file.close()
 fa1=float(sfa)
 if fa1 < fa:
   fa=fa1
   freqFA=freq
 ipolka=ipolka+1
fa1=fa
fa=math.log(fa1)
print ' *  *  * ' 
print 'Loudest event was from band ',freqFA
print 'Target log false alarm is ',fa
print ' *  *  * '
#sys.exit(0)
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

    print '======================='
    print 'Injection # ',i
    # Random injection parameters:
    cosi=random.uniform(-1.0,1.0)
    Ap=(1+cosi**2)/2 * h0
    Ac=cosi*h0
    psi=random.uniform(0.0,2*pi)
    phi01=random.uniform(0.0,2*pi)
    f0=random.uniform(freq,freq+float(coincidence_band)*(Npolka-1))
    alpha=random.uniform(0.0,2*pi)
    cosdelta=random.uniform(-1.0,1.0)
    delta=math.acos(cosdelta)
    if delta > pi/2:
      delta = delta-pi

    phi02=2*pi*f0*(gps_start2-gps_start1)  # adjust phase in second injection by time elapsed


    #fmin=int(f0*10.0)/10.0-0.1
    #band=0.2


    #print 'Minimum freq of SFTs:',fmin,'Band:',band
    sys.stdout.flush()

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


##    os.system('cat In1.data')
##    print ' '
    
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

##    os.system('cat In2.data')

    # inject the fake signal into 1st ifo data
    os.mkdir('signal1')
    sifo=sifo1
    # run makefakedata using In1.data file
    makefakedata_args=' '.join(['./makefakedata_v2','-i In1.data -n signal1/SFT -I',sifo,'-D xdata1/','-E .'])
    #print makefakedata_args
    #print 'IFO',sifo,'. Injecting fake signal in real data'
    #print ' '
	
    os.system(makefakedata_args)

    # mismatched template variables
    theta=random.uniform(0.0,2*pi)
    AngularDistance=random.uniform(0.0,math.sqrt(float(alpha_window)**2+float(delta_window)**2)/2.0)
    alphaT=alpha+AngularDistance*math.cos(theta)
    deltaT=delta+AngularDistance*math.sin(theta)
    fstart=f0-float(freq_window)
    print 'f0,fstart=',f0,fstart

    # search for the signal 
    ifo=ifo1
    cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(fstart),'-b',str(2*float(freq_window)),\
                          '--expLALDemod -I ',str(ifo),'-r',str(df),'-a',str(alphaT),'-d',str(deltaT),
                          '-D signal1','-E . -y 00-04 -F 0.0','-o -1' ])
    #print cfstat_args
    #print 'IFO',sifo,'. Searching fake signal in real data'
    #print ' '

    os.system(cfstat_args)    


    #line_kounter=os.system['cat Fstats-1 | python countlines.py']

    fstats1_file=open('Fstats-1',mode='r')
    data=fstats1_file.readlines()
    line_kounter=len(data)
    fstats1_file.close()
    if line_kounter != 2:
      print 'IFO: ',str(ifo),'. Injection and search.'
      print 'Fstats with ',line_kounter,' lines. Exiting 1'
      print '========='
      sys.exit(0)

    fstats1_file=open('Fstats-1',mode='r')
    line=fstats1_file.readline()
    #print "reading..."
    [sf1,dmp,dmp,dmp,dmp,dmp,sF1]=line.split()
    #print sf1,sF1
    fstats1_file.close()

    if float(sf1) < float(freq) or float(sf1) > (float(freq)+float(coincidence_band)*Npolka):
      print 'discarding injection.... continuing....'
      os.system('rm -rf signal1')            
      continue
    
    # inject the fake signal into 2nd ifo data
    os.mkdir('signal2')
    sifo=sifo2
    # run makefakedata using In1.data file
    makefakedata_args=' '.join(['./makefakedata_v2','-i In2.data -n signal2/SFT -I',sifo,'-D xdata2/','-E .'])
    #print makefakedata_args
    #print 'IFO',sifo,'. Injecting fake signal in real data'
    #print ' '

    
    os.system(makefakedata_args)
    # search for the signal 
    ifo=ifo2
    cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(fstart),'-b',str(2*float(freq_window)),\
                          '--expLALDemod -I',str(ifo),'-r',str(df),'-a',str(alphaT),'-d',str(deltaT),
                          '-D signal2','-E . -y 00-04 -F 0.0','-o -2'])
    #print cfstat_args
    #print 'IFO',sifo,'. Searching fake signal in real data'
    #print ' '
    os.system(cfstat_args)    


    fstats2_file=open('Fstats-2',mode='r')
    data=fstats2_file.readlines()
    line_kounter=len(data)
    fstats2_file.close()      

    if line_kounter != 2:
      print 'IFO: ',str(ifo),'. Injection and search.'
      print 'Fstats with ',line_kounter,' lines. Exiting 1'
      print '========='
      sys.exit(0)

    fstats2_file=open('Fstats-2',mode='r')
    line=fstats2_file.readline()
    [sf2,dmp,dmp,dmp,dmp,dmp,sF2]=line.split()
    fstats2_file.close()
  

    # calculate theoretical values we expect from the noise and signal parameters
    findSh_1=' '.join(['./lalapps_FindSh -D signal1 -b 0.03 -f',str(float(sf1)-0.015)])
    findSh_2=' '.join(['./lalapps_FindSh -D signal2 -b 0.03 -f',str(float(sf2)-0.015)])
                                                                                                        
    myJob = popen2.Popen3(findSh_1)
    sSh1 = myJob.fromchild.readline()
    del myJob
    myJob = popen2.Popen3(findSh_2)
    sSh2 = myJob.fromchild.readline()
    del myJob
                                                                                                            
    semiF_args1=' '.join(['./lalapps_SemiAnalyticF -a',str(alpha),'-d',str(delta),'-Q',str(phi01/2.0),\
                          '-Y',str(psi),'-i',str(cosi),'-s',str(h0),'-N',str(float(sSh1)),'-T',ts1,\
                          '-t 1800.0 -n 20 -E . -D 1'])

    #print semiF_args1
    myJob = popen2.Popen3(semiF_args1)
    sF1th = myJob.fromchild.readline()
    del myJob
 
    semiF_args2=' '.join(['./lalapps_SemiAnalyticF -a',str(alpha),'-d',str(delta),'-Q',str(phi02/2.0),\
                          '-Y',str(psi),'-i',str(cosi),'-s',str(h0),'-N',str(float(sSh2)),'-T',ts2,\
                          '-t 1800.0 -n 20 -E . -D 2'])

    #print semiF_args2
    myJob = popen2.Popen3(semiF_args2)
    sF2th = myJob.fromchild.readline()
    del myJob


    F1=float(sF1)
    F2=float(sF2)
    Fth=float(sFth)
    f1=float(sf1)
    f2=float(sf2)

    # if the results are smaller than our threshold or they are
    # not coincident in frequency   
    if F1 < Fth or F2 < Fth or abs(f1-f2) > float(freq_window):  
      print '%Maximum event is too small or not freq. coincident:', F1, F2,  abs(f1-f2)
      os.system('rm -rf signal1')
      os.system('rm -rf signal2')            
      sys.stdout.flush()
      
    # if the results are greater than our threshold and they are
    # within the coincidence window in frequency
    # then we compute the false alarm probs. 
    if F1 > Fth and F2 > Fth and abs(f1-f2) < float(freq_window):  


      Ifa1=math.log((1.0+F1/2.0))-(F1/2.0)
      Ifa2=math.log((1.0+F2/2.0))-(F2/2.0)      
      Ifa=Ifa1+Ifa2

      print 'F values for the injections: ',F1,F2
      print 'corresponding log false alarms: ' ,Ifa1,Ifa2
      print 'joint log false alarm: ',Ifa
      print 'false alarm of loudest event: ', fa
      print ' '

      if Ifa > fa:
        print '%Joint false alarm is too large:',Ifa,'Target was:',fa
        os.system('rm -rf signal1')
        os.system('rm -rf signal2')            
        sys.stdout.flush()

      
      if Ifa <= fa:

        print '... continuing...'

        if F1 > float(Fth_chisq):

          #print 'F1 larger than chi2 thereshold'
          #print ' '
          # 1) re-run ComputeFstat for a small band with -p (get pulsar parameters) option
          ifo=ifo1
          sifo=sifo1
          data='signal1'
          ts=ts1
          #sa=sa1
          #sd=sd1
          sf=sf1
          sa=alphaT
          sd=deltaT
          # starting frequency of search is chisq_points points to the left
          chi_fstart=float(sf)-float(chisq_points)*float(df) 
          chi_freq_band=2.0*float(df)*float(chisq_points)     # band is 2*chisq_points 

          # define command line for run on first ifo
          cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),\
                                '-b',str(chi_freq_band),'--expLALDemod -I',str(ifo),\
                                '-r',df,'-a',str(sa),'-d',str(sd),'-D',data,'-E . -y 00-04 -F',\
                                str(Fth),'-p'])


          #print cfstat_args
          #print 'Running search again with parameter estimation'
          #print ' '
          os.system(cfstat_args)

          fstats1_file=open('Fstats',mode='r')
          data=fstats1_file.readlines()
          line_kounter=len(data)
          fstats1_file.close()
          if line_kounter != 2:
            print 'IFO: ',str(ifo)
            print 'Fstats file has ',line_kounter,' lines. Exiting 1'
            print '===='
            #os.system('cat Fstats')
            sys.exit(0)
 
    
          # 2) run makeinvetofile makes the In.data file for makefakedata
          #    and checks that there's only one outlier
          makeinveto_args=' '.join(['./lalapps_makeInvetofile','-f Fstats -p ParamMLE -o In.data -t',ts,\
                                    '-l 1800.0 -n 20','-s',str(int(float(sf))-1),'-b 3.0'])
          #print makeinveto_args
          #print 'making In.data file for model signal...'
          #print ' '
          os.system(makeinveto_args)

          # make signal directory to put sfts with fake data
          os.mkdir('signal')
    
          # 3) run makefakedata using Inveto.data file
          makefakedata_args=' '.join(['./makefakedata_v2','-i In.data -n signal/SFT -I',sifo,'-E .'])
          #print makefakedata_args
          os.system(makefakedata_args)
          #print 'making the model signal...'
          #print ' '
    
          # 4) run ComputeFStat on fake signal data with signal only flag (-S)
          cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',\
                                str(chi_freq_band),'--expLALDemod -I ',str(ifo),'-r',\
                                str(df),'-a',str(sa),'-d',str(sd),'-D signal',\
                                '-E . -y 00-04 -F 0.0','-p -S'])
          #print  cfstat_args
          #print 'searching the model signal to construct vetoe signal'
          #print ' '
          os.system(cfstat_args)
          
          # remove signal directory
          os.system('rm -rf signal/')

          # 5) run Fstatshpetest and see if chisq test is passed
          fstatshape_args=' '.join(['./lalapps_FstatShapeTestLAL -o FaFb00.001 -t FaFb01.001 > chisq1.txt'])
          #print fstatshape_args
          os.system(fstatshape_args)
          #print 'running Fstat shape vetoe test'
          #print ' '
          
          # 6) is chi sq test passed?


          chisq_file=open('chisq1.txt',mode='r')
          data=chisq_file.readlines()
          line_kounter=len(data)
          chisq_file.close()
          if line_kounter != 1:
            print 'IFO: ',str(ifo),', chi2 test.'
            print 'chisq1.txt file with ',line_kounter,' lines. Exiting 4'
            print '========='
            sys.exit(0)

          
          chisq_file=open('chisq1.txt',mode='r')
          line=chisq_file.readline()
          [crap,crap,sdof1,schisq1]=line.split()
          chisq_file.close()
          
          chisq1=float(schisq1)
          dof1=float(sdof1)    
          chisq_th=dof1*(math.sqrt(F1)/10.0 +1)**2
#          chisq_th=(4*(2*float(chisq_points)+1)-4)*(math.sqrt(F1)/10.0 +1)**2

          #print 'Injections chisq= ',chisq1
          #print 'Threshold at ', chisq_th
          #print ' '

          if chisq1 < chisq_th:
            chisq1_pass=1
          else:
            chisq1_pass=0
            #print '%Candidate 1 did not survive chi sq test'
            sys.stdout.flush()

          # clean-up
          os.system('rm chisq1.txt FaFb0* Fstats In.data ParamMLE')
    
        else: chisq1_pass=1

        #print 'pass ?', chisq1_pass
        #print ' '
        #print 'Other IFO:,', ifo2
        #print ' '

        F2=float(sF2)
        if F2 > float(Fth_chisq):
          # run chisq test on sub-candidate
          # 1) re-run ComputeFstat for a small band with -p (get pulsar parameters) option
          ifo=ifo2
          sifo=sifo2
          data='signal2'
          ts=ts2
          #sa=sa2
          #sd=sd2
          sa=alphaT
          sd=deltaT
          sf=sf2
          # starting frequency of search is chisq_points points to the left
          chi_fstart=float(sf)-float(chisq_points)*float(df) 
          chi_freq_band=2.0*float(df)*float(chisq_points)     # band is 2*chisq_points 

          # define command line for run on first ifo
          cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',\
                                str(chi_freq_band),'--expLALDemod -I ',str(ifo),'-r',\
                                str(df),'-a',str(sa),'-d',str(sd),'-D',data,'-E . -y 00-04 -F',\
                                str(Fth),'-p'])
          #print cfstat_args
          os.system(cfstat_args)
          #print 'Running search again with parameter estimation'
          #print ' '

          fstats2_file=open('Fstats',mode='r')
          data=fstats2_file.readlines()
          line_kounter=len(data)
          fstats2_file.close()
          if line_kounter != 2:
            print 'IFO: ',str(ifo)
            print 'Fstats file has ',line_kounter,' lines. Exiting 1'
            print '===='
            #os.system('cat Fstats')
            sys.exit(0)
 

          # 2) run makeinvetofile makes the In.data file for makefakedata
          #    and checks that there's only one outlier
          makeinveto_args=' '.join(['./lalapps_makeInvetofile','-f Fstats -p ParamMLE -o In.data -t',ts,\
                                    '-l 1800.0 -n 20','-s',str(int(float(sf))-1),'-b 3.0'])
          #print makeinveto_args
          os.system(makeinveto_args)
          #print 'making In.data file for model signal...'
          #print ' '


          # make signal directory to put sfts with fake data
          os.mkdir('signal')
       
          # 3) run makefakedata using Inveto.data file
          makefakedata_args=' '.join(['./makefakedata_v2','-i In.data -n signal/SFT -I',sifo,'-E .'])
          #print makefakedata_args
          os.system(makefakedata_args)
          #print 'making the model signal...'
          #print ' '
          
    
          # 4) run ComputeFStat on fake signal data with signal only flag (-S)
          cfstat_args=' '.join(['./lalapps_ComputeFStatistic','-f',str(chi_fstart),'-b',\
                                str(chi_freq_band),'--expLALDemod -I ',str(ifo),'-r',\
                                str(df),'-a',str(sa),'-d',str(sd),'-D signal',\
                                '-E . -y 00-04 -F 0.0','-p -S'])
          #print cfstat_args
          os.system(cfstat_args)
          #print 'searching the model signal to construct vetoe signal'
          #print ' '
          
          
          # remove signal directory
          os.system('rm -rf signal/')

          # 5) run Fstatshpetest and see if chisq test is passed
          fstatshape_args=' '.join(['./lalapps_FstatShapeTestLAL -o FaFb00.001 -t FaFb01.001 > chisq2.txt'])
          #print fstatshape_args
          os.system(fstatshape_args)
          #print 'running Fstat shape vetoe test'
          #print ' '

          # 6) is chi sq test passed?

          chisq_file=open('chisq2.txt',mode='r')
          data=chisq_file.readlines()
          line_kounter=len(data)
          chisq_file.close()
          if line_kounter != 1:
            print 'IFO: ',str(ifo),', chi2 test.'
            print 'chisq2.txt file with ',line_kounter,' lines. Exiting 5'
            print '========='
            sys.exit(0)

          
          chisq_file=open('chisq2.txt',mode='r')
          line=chisq_file.readline()
          [crap,crap,sdof2,schisq2]=line.split()
          chisq_file.close()
     

          chisq2=float(schisq2)
          dof2=float(sdof2)          
          chisq_th=dof2*(math.sqrt(F2)/10.0 +1)**2 
#          chisq_th=(4*(2*float(chisq_points)+1)-4)*(math.sqrt(F2)/10.0 +1)**2

          #print 'Injections chisq= ',chisq2
          #print 'Threshold at ', chisq_th
          #print ' '

          if chisq2 < chisq_th:
            chisq2_pass=1
          else:
            chisq2_pass=0
            #print '%Candidate 2 did not survive chi sq test'
            sys.stdout.flush()

          # clean-up
          os.system('rm chisq2.txt FaFb0* Fstats In.data ParamMLE')
       
        else: chisq2_pass=1
     
        #print 'pass ?', chisq2_pass
        #print ' '
             
        #Have both candidates passed the chi_sq veto (or not qualified to take it)?
        if chisq1_pass and chisq2_pass:
          counter=counter+1

        sys.stdout.flush()
        os.system('rm -rf signal1 signal2')            
        print f0,alpha,delta,float(sF1), 2*float(sF1th), float(sF2), 2*float(sF2th),abs(float(sf1)-float(sf2)),h0,i+1,counter,float(counter)/float(i+1)
        
      #end if Ifa <= fa  
    #end if injected signals are in coincidence  
    #increment i, independently of whether fa was big enough 
    i=i+1
    print 'counter,i= ',counter,i

  #end i-th injection, next injection
  confidenceOLD=confidence  
  confidence=float(counter)/float(Ninj)
      # first ifo
  res_out=''.join(['Confidence.data-',str(freq)])
  outCdata_file=open(res_out,mode='a')
  print >>outCdata_file,h0,confidence
  outCdata_file.close()
  shutil.copy(res_out,starting_dir)
   
  print 'confidence =',confidence  
  if abs(confidence-c0) < tol:
    if Ninj < 3000:
      Ninj=Ninj*10
      tol=tol/math.sqrt(10.0)
      cont=1
    else: cont=0
  if confidence > c0:
      if confidenceOLD < c0:
        dh0=dh0/2
      h0=h0-dh0
  if confidence < c0:
      if confidenceOLD > c0:
        dh0=dh0/2
      h0=h0+dh0

print 'Ninj,tolerance,dh0= ',Ninj,tol,dh0                
inj_data_file.close()
shutil.copy(inj_out,starting_dir)


# -------------------------------------------------------------------------------- #
