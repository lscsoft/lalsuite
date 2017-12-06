#!/usr/bin/env python2
"""
Upper_Limit_V2.py - Monte Carlo Upper Limit script - for use with targeted searches using ComputeFStatistic_v2
J. Betzwieser
"""

# Import standard modules to the python path
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
  -j, --job-id             job ID #
  -S  --starting-dir       Starting directory (location of ComputeFStatistic, ephemeris files etc...)
  -W  --local-work-dir     Local working directory (on the nodes)
  -p  --params-file        Search parameter configuration file         
"""
  print >> sys.stderr, msg

# ------------------- Parse the command line options -------------------- #

# Initialise command line argument variables
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
  if o in ("-j", "--job-id"): #An id marker for this particular job - can be traced back to determine frequency and other parameters
    job_id = int(a)
  elif o in ("-S", "--starting-dir"): #Directory in which the executables and necessary data files are located
    starting_dir = a      
  elif o in ("-W", "--local-work-dir"): #Directory on the node in which temporary files will be kept
    local_work_dir = a      
  elif o in ("-p", "--params-file"): #File containing necessary parameters of the search and MC injections
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
alpha = cp.get('fstat-params','alpha')
delta = cp.get('fstat-params','delta')
fixed_angles = cp.get('fstat-params','fixed_angle')
if fixed_angles == 'True':
  siota_mean = cp.get('fstat-params','iota_mean')
  iota_mean = float(siota_mean)
  siota_sigma = cp.get('fstat-params','iota_sigma')
  iota_sigma = float(siota_sigma)
  spsi_mean = cp.get('fstat-params','psi_mean')
  psi_mean = float(spsi_mean)
  spsi_sigma = cp.get('fstat-params','psi_sigma')
  psi_sigma = float(spsi_sigma)
sfreq = cp.get('fstat-params','freq')
freq = float(sfreq)

# read data params
sadd_band = cp.get('data-params','add_band')
add_band = float(sadd_band)
sjob_group_step = cp.get('data-params','job_group_step')
job_group_step = float(sjob_group_step)
sjob_band = cp.get('data-params','job_band')
job_band = float(sjob_band)

# read file names
loudest_dir = cp.get('files-names','loudest_dir')
loudest_files = cp.get('files-names','loudest_files')
sNfiles = cp.get('files-names','Nfiles')
Nfiles = float(sNfiles)

# read ifo parameters
sNifo = cp.get('ifo-params','Nifo')
Nifo = float(sNifo)
ifo1 = cp.get('ifo-params','ifo1')
ifo2 = cp.get('ifo-params','ifo2')
ifo3 = cp.get('ifo-params','ifo3')
data1 = cp.get('ifo-params','data1')
data2 = cp.get('ifo-params','data2')
data3 = cp.get('ifo-params','data3')

# read injection parameters
ref_time = cp.get('inj-params','ref_time')
sfdot_start = cp.get('inj-params','fdot_start')
fdot_start = float(sfdot_start)
sfdot_end = cp.get('inj-params','fdot_end')
fdot_end = float(sfdot_end)
sfdotdot_start = cp.get('inj-params','fdotdot_start')
fdotdot_start = float(sfdotdot_start)
sfdotdot_end = cp.get('inj-params','fdotdot_end')
fdotdot_end = float(sfdotdot_end)
sinj_band = cp.get('inj-params','inj_band')
inj_band = float(sinj_band)

# read search parameters
sdFreq = cp.get('search-params','dFreq')
dFreq = float(sdFreq)
sdfdot = cp.get('search-params','dfdot')
dfdot = float(sdfdot)
sdfdotdot = cp.get('search-params','dfdotdot')
dfdotdot = float(sdfdotdot)
sfreq_window = cp.get('search-params','freq_window')
freq_window = float(sfreq_window)
sfdot_window = cp.get('search-params','fdot_window')
fdot_window = float(sfdot_window)
sfdotdot_window = cp.get('search-params','fdotdot_window')
fdotdot_window = float(sfdotdot_window)

# read mc injection loop parameters
sNinj = cp.get('mc-params','Ninj')
Ninj = int(sNinj)
sNinjmax = cp.get('mc-params','Ninjmax')
Ninjmax = int(sNinjmax)
sh0 = cp.get('mc-params','h0')
h0 = float(sh0)
sdh0 = cp.get('mc-params','dh0')
dh0 = float(sdh0)
sc0 = cp.get('mc-params','c0')
c0 = float(sc0)

#Setting up some necessary starting parameters

node_id = job_id*job_group_step
if (job_group_step > 0):
  narrow_band_min = freq + job_id*job_band - (add_band/2)
  narrow_band_max = freq + (job_id+1)*job_band + (add_band/2)
  inj_fmin = freq + (job_id+(1/2))*job_band - (inj_band/2)
else:
  narrow_band_min = freq - (add_band/2)
  narrow_band_max = freq + job_band + (add_band/2)
  inj_fmin = freq + job_band/2 - (inj_band/2)

tol = 0.9/math.sqrt(Ninj)
checkpoint_start = 0
iteration = 0
injections_per_checkpoint = 20


# -------------------------------------------------------------------------------- #

# ------------------------- Check checkpoint and results files ------------------- # 

# Name of the results file
res_out=''.join([starting_dir,'/Confidence.data-',str(job_id)])
print '====='
print res_out

# Name of a temporary checkpoint file, stored in the starting directory, in case the job gets interrupted, it can be continued with a minimal loss of time
checkpoint_out =''.join([starting_dir,'/Checkpoint.data-',str(job_id)])
print '====='
print checkpoint_out

print '====='
cont=1
confidence=0.0



#Checks to see if the results file exists and if its already complete (i.e. %DONE)
if os.path.exists(res_out):
  print 'file exists! reading ...'
  Cdata_file=open(res_out,mode='r')
  Cdata_file.seek(-6,2)
  done=Cdata_file.read(5)
  endOfFile=Cdata_file.tell()
  Cdata_file.close()
  
  if done == "%DONE":
    print 'Confidence.data-',str(job_id),' already exists, thus exiting...\n '
    sys.exit(0)
  #Endif the file is complete.  


  #Open the results file and read in the last set of results to determine where to proceed from
  Cdata_file=open(res_out,mode='r')
  line_list=Cdata_file.readlines()
  length=len(line_list)
  print 'length=',length
  print line_list
  if length > 0:
    print 'reading line...'
    line=line_list[length-1]
    print line
    [siteration,sNinj,stol,sh0,sdh0,sconfidence]=line.split(None,6) 
    iteration = int(siteration)
    confidence=float(sconfidence)
    dh0=float(sdh0)
    h0=float(sh0)
    Ninj=int(sNinj)
    tol=0.9/math.sqrt(Ninj)

  #If more than 1 line, look at previous to determine where we should go from here
  if length > 1:
    line=line_list[length-2]
    [siteration,sNinjOLD,stolOLD,sh0OLD,sdh0OLD,sconfidenceOLD]=line.split(None,6)
    confidenceOLD=float(sconfidenceOLD)

    if abs(confidence-c0) < tol:  #If (current confidence - desired confidence) < tolerance, then increase number of injections and decrease tolerance
      if Ninj < Ninjmax: #If we're not yet at the desired number of injections, increase the injection factor by 2 (and modify tolerance for the new number of injections)
        Ninj=Ninj*2
        tol=tol/math.sqrt(2.0)
        cont=1
      else: cont=0 #If we are at or above the desired number of injections, then we're done and should not continue

    if confidence > c0: #If current confidence is above desired confidence
      if confidenceOLD < c0: #And if the previous confidence is below desired confidence, then we passed the optimal value and should reduce our step size
        dh0=dh0/2 
      h0=h0-dh0 #Decrease the h0 injection strength

    if confidence < c0: #If current confidence is below desired confidence
      if confidenceOLD > c0: #And the previous confidence is above desired confidence, then we passed the optimal value and should reduce our step size
        dh0=dh0/2
      h0=h0+dh0 #Increase the h0 injection strength 

  Cdata_file.close()
  #End of loading the Confidence.data file

#Print to screen indicating current number of injection, tolerance and incremental h0   
print 'Ninj,tolerance,dh0= ',Ninj,tol,dh0

#Checks to see if a checkpoint file exists - having stored the state of the injections midway through a cycle
if os.path.exists(checkpoint_out):
  print 'Checkpoint file exists! reading ...'
   

  #Load in the data from the checkpoint file and set the starting parameters
  CheckpointData_file=open(checkpoint_out,mode='r')
  line_list=CheckpointData_file.readlines()
  length=len(line_list)
  print 'length=',length
  print line_list
  if length > 0:
    print 'reading line...'
    line=line_list[length-1]
    print line
    [si,scounter,slast_checkpoint,siteration,sNinj,stol,sh0,sdh0]=line.split(None,7)
    starting_iteration = int(siteration) 
    if starting_iteration > iteration:  #Sanity check that this checkpoint file is in fact for injections after the last Confidence file results
      Ninj=int(sNinj)
      i=int(si)
      counter=int(scounter)
      last_checkpoint = int(slast_checkpoint)
      h0 = float(sh0)
      dh0 = float(sdh0)
      tol=0.9/math.sqrt(Ninj)
      checkpoint_start = 1
  CheckpointData_file.close()


# -------------------------------------------------------------------------------- #

# ------------------------- Check search results files --------------------------- # 

# check that we've got all the data files that we need

idata=0
while idata < job_group_step:

 # name of loudest output files from ComputeFStatistic_V2
 current_node_id = int(node_id) + int(idata)
 loudest_in_file=''.join([loudest_dir,loudest_files,str(current_node_id)])
 if os.path.exists(loudest_in_file) != 1:
  print "Not all necessary files are present"
  print "could not find ",loudest_in_file
  print "Proceeding to next band...."
  print "...."
  sys.exit(0)
 idata = idata + 1


# -------------------------------------------------------------------------------- #

# ------------------------- Set up working directory ----------------------------- # 


# Name of local work directory
subdir=''.join([local_work_dir,'/','ul-run.',str(job_id)])


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
cfstat=''.join([starting_dir,'/ComputeFStatistic_v2'])
mkdata=''.join([starting_dir,'/lalapps_Makefakedata'])
extract_data=''.join([starting_dir,'/ConvertToSFTv2'])

#     ephemeris file
earth=''.join([starting_dir,'/earth03-06.dat'])
sun=''.join([starting_dir,'/sun03-06.dat'])

#      Confidence file
if os.path.exists(res_out):
  shutil.copy(res_out,subdir)

# copy execulables and components to working sub-directory 
shutil.copy(cfstat,subdir)
shutil.copy(mkdata,subdir)
shutil.copy(extract_data,subdir)
shutil.copy(earth,subdir)
shutil.copy(sun,subdir)


# copy over and unzip input polka files

idata=0
while idata < Nfiles:
 # name of output files from ComputeFStatistic_V2
 current_node_id = int(node_id) + int(idata)
 loudest_in_file=''.join([loudest_dir,loudest_files,str(current_node_id)])

 shutil.copy(loudest_in_file,subdir)
 idata=idata+1

#Done with directory setup on the node
print 'Finished Setting up working directory.'
sys.stdout.flush()

# -------------------------------------------------------------------------------- #

# -------------- extract and copy relevant data onto local dir  ------------------- # 



# 1st ifo
datadir1 = './xdata1/'
os.mkdir(datadir1)

narrow_band = ''.join(['./ConvertToSFTv2 --inputSFTs ',data1,' -o ',datadir1,' --fmin ',str(narrow_band_min),' --fmax ',str(narrow_band_max)])
os.system(narrow_band)

if Nifo > 1:
  # 2nd ifo
  datadir2 = './xdata2/'
  os.mkdir(datadir2)

  narrow_band = ''.join(['./ConvertToSFTv2 --inputSFTs ',data2,' -o ',datadir2,' --fmin ',str(narrow_band_min),' --fmax ',str(narrow_band_max)])
  os.system(narrow_band)

if Nifo > 2:
  # 3rd ifo
  datadir3 = './xdata3/'
  os.mkdir(datadir3)

  narrow_band = ''.join(['./ConvertToSFTv2 --inputSFTs ',data3,' -o ',datadir3,' --fmin ',str(narrow_band_min),' --fmax ',str(narrow_band_max)])
  os.system(narrow_band)


# -------------------------------------------------------------------------------- #

# --------------------------  Target 2F Value  -------------------------------- #

idata=0
twoF_loudest = 0
while idata < Nfiles:
 sys.stdout.flush()

 # name of loudest output files from ComputeFStatistic_V2
 current_node_id = int(node_id) + int(idata)
 loudest_in_file=''.join([loudest_dir,loudest_files,str(current_node_id)])

 sys.stdout.flush()
 results_file=open(loudest_in_file,mode='r')
 reading = 0
 sys.stdout.flush()
 while reading < 36:
  sys.stdout.flush()  
  line=results_file.readline()
  reading = reading + 1
 
 words=line.split()
 twoF_loudest1=float(words[2].strip(';'))
 results_file.close()
 if twoF_loudest1 > twoF_loudest: #If new 2F value is higher, set that as the 2F largest value
   twoF_loudest = twoF_loudest1
 idata=idata+1


print 'Target 2F value is ',twoF_loudest
print ' *  *  * '
sys.stdout.flush()

while cont:
  if checkpoint_start: #Are we starting from a checkpointed spot
    checkpoint_start = 0  #If so, don't reset counters
  else: #If we're starting a new cycle, reset
    i = 0
    counter = 0
    last_checkpoint = 0
    iteration = iteration + 1
  
  while i < Ninj:
    
    if fixed_angles == 'True':
      w = 2.0
      while w >= 1.0:  
        x1 = 2.0 * random.uniform(0.0,1.0) - 1.0
        x2 = 2.0 * random.uniform(0.0,1.0) - 1.0
        w = x1 * x1 + x2 * x2
      w = math.sqrt((-2.0 * math.log(w))/ w)
      y1 = x1 * w
      y2 = x2 * w

      iota = iota_mean + y1 * iota_sigma
      cosi = math.cos(iota)
      psi = psi_mean + y2 * psi_sigma
    else:
      cosi = random.uniform(-1.0,1.0)
      psi = random.uniform(0.0,2*pi)

    Ap=(1+cosi**2)/2 * h0
    Ac=cosi*h0
    phi0=random.uniform(0.0,2*pi)
    f0=random.uniform(freq,freq+job_band)
    
    fdot=random.uniform(fdot_start,fdot_end)
    
    fdotdot=random.uniform(fdotdot_start,fdotdot_end)

    # inject the fake signal into 1st ifo data
    os.mkdir('signal1')

    sifo=ifo1
    #run makefakedata
    input_dir = ''.join(['\'',subdir,'/xdata1/*\''])
    input_dir = input_dir.strip(' ')
    makefakedata_args=' '.join(['./lalapps_Makefakedata','-n signal1/ -I ',str(sifo),' -E . -y 03-06 --fmin ',str(inj_fmin),' --Band ',str(inj_band),' --Tsft 1800.0 --Alpha ', str(alpha),' --Delta ', str(delta), ' --aPlus ', str(Ap),' --aCross ',str(Ac),' --psi ',str(psi),' --phi0 ',str(phi0),'-S ',str(ref_time), ' --Freq ',str(f0),' --f1dot ',str(fdot), '--f2dot ',str(fdotdot), ' -D ',input_dir])

    os.system(makefakedata_args)
    
    if Nifo > 1:
    
      # inject the fake signal into 2nd ifo data
      sifo=ifo2
      input_dir = ''.join(['\'',subdir,'/xdata2/*\''])
      input_dir = input_dir.strip(' ')
      # run makefakedata 
      makefakedata_args=' '.join(['./lalapps_Makefakedata','-n signal1/ -I ',str(sifo),' -E . -y 03-06 --fmin ',str(inj_fmin),' --Band ',str(inj_band),' --Tsft 1800.0 --Alpha ', str(alpha),' --Delta ', str(delta), ' --aPlus ', str(Ap),' --aCross ',str(Ac),' --psi ',str(psi),' --phi0 ',str(phi0),'-S ',str(ref_time), ' --Freq ',str(f0),' --f1dot ',str(fdot), '--f2dot ',str(fdotdot), ' -D ',input_dir])
 
      os.system(makefakedata_args)
      

    if Nifo > 2:

      # inject the fake signal into 3rd ifo data
      sifo=ifo3
      input_dir = ''.join(['\'',subdir,'/xdata3/*\''])
      input_dir = input_dir.strip(' ')
      # run makefakedata 
      makefakedata_args=' '.join(['./lalapps_Makefakedata','-n signal1/ -I ',str(sifo),' -E . -y 03-06 --fmin ',str(inj_fmin),' --Band ',str(inj_band),' --Tsft 1800.0 --Alpha ', str(alpha),' --Delta ', str(delta), ' --aPlus ', str(Ap),' --aCross ',str(Ac),' --psi ',str(psi),' --phi0 ',str(phi0),'-S ',str(ref_time), ' --Freq ',str(f0),' --f1dot ',str(fdot), '--f2dot ',str(fdotdot), ' -D ',input_dir])
      
      os.system(makefakedata_args)
 

    # mismatch the template variables

    half_dFreq = float(dFreq/2)
    half_dfdot = float(dfdot/2)
    half_dfdotdot = float(dfdotdot/2)
    Freq_vec = random.uniform(-half_dFreq,half_dFreq)
    fdot_vec = random.uniform(-half_dfdot,half_dfdot)
    fdotdot_vec = random.uniform(-half_dfdotdot,half_dfdotdot)
    Freq_first = f0 + Freq_vec - freq_window
    Freq_band = 2*freq_window
    fdot_first = fdot + fdot_vec - fdot_window
    fdot_band = 2*fdot_window
    fdotdot_first = fdotdot + fdotdot_vec - fdotdot_window
    fdotdot_band = 2*fdotdot_window

    #Search for the injected signal with ComputeFStatistic_v2

   
    cfstat_args=' '.join(['./ComputeFStatistic_v2','-a ',str(alpha),'-d ',str(delta),' -f ',str(Freq_first),' -b ', str(Freq_band),' -s ',str(fdot_first),' -m ',str(fdot_band), '--f2dot ',str(fdotdot_first),  '--f2dotBand ',str(fdotdot_band),' --dFreq ',str(dFreq),' --df1dot ',str(dfdot),' --df2dot ',str(dfdotdot),' --DataFiles=./signal1/* -E . -y 03-06 -F 0.0 --outputFstat Fstat1  --outputLoudest loudestFound'])

    os.system(cfstat_args)    
  
    #Open file containing loudest event found    
    fstats1_file=open('loudestFound',mode='r')
    reading = 0
    sys.stdout.flush()
    while reading < 36:
     sys.stdout.flush()  
     line=fstats1_file.readline()
     reading = reading + 1
 
    words=line.split()
    sF1=float(words[2].strip(';'))
    fstats1_file.close()


    F1=float(sF1)

    #If the largest value found in the injection is less than or equal to our largest value found in the actual search, we don't count it   
    if F1 <= twoF_loudest:
      os.system('rm -rf signal1')            
      sys.stdout.flush()

    #If the largest value found in the injection is greater than our largest value found in the actual search, we count it
    if F1 > twoF_loudest:

      counter=counter+1

      sys.stdout.flush()
      os.system('rm -rf signal1')            

    #increment i (tracking the number of injections) 
    i=i+1
    
    #If we've done injections_per_checkpoint number of injections, add current status of injections to checkpoint file
    if (i/injections_per_checkpoint) > last_checkpoint:
      last_checkpoint = last_checkpoint + 1
      checkpoint_out=''.join(['Checkpoint.data-',str(job_id)])
      outCheckpointData_file = open(checkpoint_out,mode='a')
      print >>outCheckpointData_file,i,counter,last_checkpoint,iteration,Ninj,tol,h0,dh0
      outCheckpointData_file.close()
      shutil.copy(checkpoint_out,starting_dir)
  
  #end i-th injection, next injection
  
  #Save confidence results in file

  confidenceOLD=confidence  
  confidence=float(counter)/float(Ninj) #Calculate confidence

  res_out=''.join(['Confidence.data-',str(job_id)])
  outCdata_file=open(res_out,mode='a')
  print >>outCdata_file,iteration,Ninj,tol,h0,dh0,confidence
  outCdata_file.close()
  shutil.copy(res_out,starting_dir)
  rm_chk = ''.join(['rm ',starting_dir,'/Checkpoint.data-',str(job_id)])
  rm_local_chk = ''.join(['rm ',subdir,'/Checkpoint.data-',str(job_id)])
  os.system(rm_chk)
  os.system(rm_local_chk)    
 
  print 'confidence =',confidence  

  #Determine if we need to do another cycle of injections with a different h0 and/or more injections to improve accuracy
  if abs(confidence-c0) < tol:
    if Ninj < Ninjmax:
      Ninj=Ninj*2
      tol=tol/math.sqrt(2.0)
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

# reopen file to insert the "DONE" flag
res_out=''.join(['Confidence.data-',str(job_id)])
outCdata_file=open(res_out,mode='a')
print >>outCdata_file,'%DONE'
outCdata_file.close()

os.system(rm_subdir)



# -------------------------------------------------------------------------------- #
