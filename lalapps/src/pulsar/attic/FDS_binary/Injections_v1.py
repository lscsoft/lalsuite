#!/usr/bin/env python2
"""
This script has been written to inject signals inito a given dataset
and extract them using the search pipeline.  This version injects
signals from sources in eccentric orbits.
"""

# import required modules
import sys
import os
import getopt
import re
import shutil
import string
import tempfile
import math
import time 
import ConfigParser

# append the lalapps python path
sys.path.append('@PYTHONLIBDIR@')

# program usage
def usage():
  msg = """\
Usage: GenerateFreqMeshFile [options]

  -h, --help               display this message
  -c, --configfile         the name of the input configuration file
  -f, --fmin               the minimum injection frequency
  -F, --fmax               the maximum injection frequency
  -H, --gwamplitude        the gw amplitude
  -e, --eccentricity       the orbital eccentricity
  -p, --pdatadir           the primary detector dataset for injecting into
  -s, --sdatadir           the secondary detector dataset for injecting into
  -T, --ptemplatefile      the name of the primary orbital template file
  -t, --stemplatefile      the name of the secondary orbital template file
  -n, --ntrials            the number of injections required
  -O, --outdir             the output directory for all injection results 
  -x, --id                 an ID number to help identify this set of injections
  """
  print >> sys.stderr, msg

# parse the command line options to figure out what we should do
shortop = "hc:f:F:H:e:p:s:T:t:n:O:x:"
longop = [
"help",
"configfile=",
"fmin=",
"fmax=",
"gwamplitude=",
"eccentricity=",
"pdatadir=",
"sdatadir=",
"ptemplatefile=",
"stemplatefile=",
"ntrials=",
"outdir=",
"id=",
]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

# default options
configfile = None
fmin = None
fmax = None
h0 = None
ecc = None
pdatadir = None
sdatadir = None
ptemplatefile = None
stemplatefile = None
ntrials = None
outdir = None
id = None

# process options
for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-c", "--configfile"):
    configfile = a
  elif o in ("-f", "--fmin"):
    fmin = a
  elif o in ("-F", "--fmax"):
    fmax = a
  elif o in ("-H", "--gwamplitude"):
    h0 = a
  elif o in ("-e", "--eccentricity"):
    ecc = (float)(a)
  elif o in ("-p", "--pdatadir"):
    pdatadir = a
  elif o in ("-s", "--sdatadir"):
    sdatadir = a
  elif o in ("-T", "--ptemplatefile"):
    ptemplatefile = a
  elif o in ("-t", "--stemplatefile"):
    stemplatefile = a
  elif o in ("-n", "--ntrials"):
    ntrials = a
  elif o in ("-O", "--outdir"):
    outdir = a
  elif o in ("-x", "--id"):
    id = a
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

if not configfile:
  print >> sys.stderr, "No configfile file specified."
  print >> sys.stderr, "Use --configfile FILE to specify location."
  sys.exit(1)

if not fmin:
  print >> sys.stderr, "No minimum frequency specified."
  print >> sys.stderr, "Use --fmin REAL8 to specify a value."
  sys.exit(1)

if not fmax:
  print >> sys.stderr, "No maximum frequency specified."
  print >> sys.stderr, "Use --fmax REAL8 to specify a value."
  sys.exit(1)

if not h0:
  print >> sys.stderr, "No h0 value specified."
  print >> sys.stderr, "Use --h0 REAL8 to specify a value."
  sys.exit(1)

if not ecc:
  print >> sys.stderr, "No eccentricity specified."
  print >> sys.stderr, "Setting ecc = 0.0."
  ecc = 0.0

if not pdatadir:
  print >> sys.stderr, "No primary data directory specified."
  print >> sys.stderr, "Use --pdatadir PATH to specify a location."
  sys.exit(1)

if not sdatadir:
  print >> sys.stderr, "No secondary data directory specified."
  print >> sys.stderr, "Use --sdatadir PATH to specify a location."
  sys.exit(1)

if not ptemplatefile:
  print >> sys.stderr, "No primary template file specified."
  print >> sys.stderr, "Use --ptemplatefile FILE to specify a location."
  sys.exit(1)

if not stemplatefile:
  print >> sys.stderr, "No secondary template file specified."
  print >> sys.stderr, "Use --stemplatefile FILE to specify a location."
  sys.exit(1)

if not ntrials:
  print >> sys.stderr, "No number of injections specified."
  print >> sys.stderr, "Use --ntrials INT4 to specify a value."
  sys.exit(1)

if not outdir:
  print >> sys.stderr, "No output directory specified."
  print >> sys.stderr, "Use --outdir PATH to specify a location."
  sys.exit(1)

if not id:
  print >> sys.stderr, "No ID number specified."
  print >> sys.stderr, "Use --id INT4 to specify a value."
  sys.exit(1)

######################################################################
# read input config file

# create the config parser object and read in the ini file
cp = ConfigParser.ConfigParser()
cp.read(configfile)

sensitivity_code = cp.get('condor','sensitivity')
makemesh_code = cp.get('condor','makemesh')
search_code = cp.get('condor','search')
coincidence_code = cp.get('condor','coincidence')
injections_code = cp.get('condor','injections')
upperlimit_code = cp.get('condor','upperlimit')
makefakedata_code = cp.get('condor','makefakedata')
makesubmesh_code = cp.get('condor','makesubmesh')
binaryinput_code = cp.get('condor','binaryinput')
binaryrandominput_code = cp.get('condor','binaryrandominput')
getdataparams_code = cp.get('condor','getdataparams')
selectdata_code = cp.get('condor','selectdata')
freqmeshfile_code = cp.get('condor','freqmeshfile')
periapseshift_ecc_code = cp.get('condor','periapseshift_ecc')
ephdir = cp.get('ephemeris','ephdir')
yr = cp.get('ephemeris','yr')
fulldata_primary = cp.get('data','fulldata_primary')
fulldata_secondary = cp.get('data','fulldata_secondary')
tsft = (int)(cp.get('data','tsft'))
resultsdir = cp.get('results','resultsdir')
det_primary = cp.get('detectors','det_primary')
det_secondary = cp.get('detectors','det_secondary')
tspan = cp.get('sensitivity','tspan')
tstep = cp.get('sensitivity','tstep')
sensitivity_nodes = (int)(cp.get('nodes','sensitivity_nodes'))
search_nodes = (int)(cp.get('nodes','search_nodes'))
coincidence_nodes = (int)(cp.get('nodes','coincidence_nodes'))
injections_nodes = (int)(cp.get('nodes','injections_nodes'))
sourcefile = cp.get('source','sourcefile')
source = cp.get('source','source')
dterms = cp.get('search','dterms')
windowsize = cp.get('search','windowsize')
dopplermax = cp.get('search','dopplermax')
overres = cp.get('search','overres')
binary = cp.get('search','binary')
thresh_primary = cp.get('thresholds','thresh_primary')
thresh_secondary = cp.get('thresholds','thresh_secondary')
mismatch = cp.get('makemesh','mismatch')
meshband = cp.get('makemesh','band')
co_nbins = cp.get('coincidence','nbins')
h0min = (float)(cp.get('injections','h0min'))
h0max = (float)(cp.get('injections','h0max'))
h0step = (float)(cp.get('injections','h0step'))
injectionband = cp.get('injections','injectionband')
tempworkarea = cp.get('injections','tempworkarea')
fband = cp.get('injections','fband')
confidence = cp.get('upperlimit','confidence')
minres = cp.get('general','minres')

# define arg-periapse min and max values
argpmin = 1.0*math.pi
argpmax = -1.0*math.pi
#argpmin = 1.5
#argpmax = 1.5

######################################################################
# determine where we are working (BE MOST TOTALLY CAREFUL)

# define current directory
pwd = os.getcwd()
#print "present working directory is %s" %(pwd)

# testing only
#os.system('ls -l /scratch/tmp')
#temp1 = "/scratch/tmp/cm"
#temp2 = "/scratch/tmp/chrism"
#if os.path.exists(temp1):
  #os.system('ls -l /scratch/tmp/cm')
#if os.path.exists(temp2):
  #os.system('ls -l /scratch/tmp/chrism')

# check existence of and make temporary directory
if not os.path.exists(tempworkarea):
  try: os.mkdir(tempworkarea)
  except:
    print "ERROR : failed to create temporary work area [%s], exiting." %(tempworkarea)
    sys.exit(1)

# defining specific temporary work area directory name
tempworkarea = "%s/injections_%.3f-%.3f_%s_%.3e_%s" %(tempworkarea,float(fmin),float(fmax),h0,ecc,id)

#print tempworkarea

# checking if the temporary work area exists
if os.path.exists(tempworkarea):
  print "WARNING : temporary work area %s already exists, removing it !!!" %(tempworkarea)
  print "WARNING : deleting contents of %s directory !!!" %(tempworkarea)
  try: shutil.rmtree(tempworkarea)
  except:
    print "ERROR : failed to delete contents of %s directory !!!" %(tempworkarea)
    sys.exit(1)

# making the temporary work area
try: os.mkdir(tempworkarea)
except:
  print "ERROR : failed to make temporary work area %s !!!" %(tempworkarea)
  sys.exit(1)

# move to that directory
try: os.chdir(tempworkarea)
except:
  print "ERROR : failed to change directory to %s !!!" %(tempworkarea)
  sys.exit(1)

######################################################################
# make all of the directories required in the work area

# define some directory names
p_datasetdir = tempworkarea + '/p_datasetdir'
s_datasetdir = tempworkarea + '/s_datasetdir'
p_searchresultsdir = tempworkarea + '/p_searchresultsdir'
s_searchresultsdir = tempworkarea + '/s_searchresultsdir'
co_resultsdir = tempworkarea + '/co_resultsdir'

# make these directories
try: os.mkdir(p_datasetdir)
except:
  print "ERROR : failed to create directory %s !!!" %(p_datasetdir)
  sys.exit(1)

try: os.mkdir(s_datasetdir)
except:
  print "ERROR : failed to create directory %s !!!" %(s_datasetdir)
  sys.exit(1)  

try: os.mkdir(p_searchresultsdir)
except:
  print "ERROR : failed to create directory %s !!!" %(p_searchresultsdir)
  sys.exit(1)

try: os.mkdir(s_searchresultsdir)
except:
  print "ERROR : failed to create directory %s !!!" %(s_searchresultsdir)
  sys.exit(1)

try: os.mkdir(co_resultsdir)
except:
  print "ERROR : failed to create directory %s !!!" %(co_resultsdir)
  sys.exit(1)

######################################################################
# copy the required files to the working location

binaryinput_temp = tempworkarea + '/binaryinput'
binaryrandominput_temp = tempworkarea + '/binaryrandominput'
makemesh_temp = tempworkarea + '/makemesh'
makesubmesh_temp = tempworkarea + '/makesubmesh'
search_temp = tempworkarea + '/search'
makefakedata_temp = tempworkarea + '/makefakedata'
coincidence_temp = tempworkarea + '/coincidence'
getdataparams_temp = tempworkarea + '/getdataparams'
freqmeshfile_temp = tempworkarea + '/freqmeshfile'
sourcefile_temp = tempworkarea + '/sourcefile'
periapseshift_ecc_temp = tempworkarea + '/periapseshift'

try: shutil.copyfile(binaryinput_code,binaryinput_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(binaryinput_code,binaryyinput_temp)
  sys.exit(1)

try: shutil.copyfile(binaryrandominput_code,binaryrandominput_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(binaryrandominput_code,binaryrandominput_temp)
  sys.exit(1)

try: shutil.copyfile(makemesh_code,makemesh_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(makemesh_code,makemesh_temp)
  sys.exit(1)
  
try: shutil.copyfile(makesubmesh_code,makesubmesh_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(makesubmesh_code,makesubmesh_temp)
  sys.exit(1)
  
try: shutil.copyfile(search_code,search_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(search_code,search_temp)
  sys.exit(1)

try: shutil.copyfile(makefakedata_code,makefakedata_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(makefakedata_code,makefakedata_temp)
  sys.exit(1)

try: shutil.copyfile(coincidence_code,coincidence_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(coincidence_code,coincidence_temp)
  sys.exit(1)

try: shutil.copyfile(getdataparams_code,getdataparams_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(getdataparams_code,getdataparams_temp)
  sys.exit(1)

try: shutil.copyfile(freqmeshfile_code,freqmeshfile_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(freqmeshfile_code,freqmeshile_temp)
  sys.exit(1)

try: shutil.copyfile(sourcefile,sourcefile_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(sourcefile,sourcefile_temp)
  sys.exit(1)

try: shutil.copyfile(periapseshift_ecc_code,periapseshift_ecc_temp)
except:
  print "ERROR : failed to copy %s to %s !!!" %(periapseshift_ecc_code,periapseshift_ecc_temp)
  sys.exit(1)

# make all the copies of the codes executable (its a bit clunky using os.system) 
os.system('chmod +x ' + binaryinput_temp)
os.system('chmod +x ' + binaryrandominput_temp)
os.system('chmod +x ' + makemesh_temp)
os.system('chmod +x ' + makesubmesh_temp)
os.system('chmod +x ' + search_temp)
os.system('chmod +x ' + makefakedata_temp)
os.system('chmod +x ' + coincidence_temp)
os.system('chmod +x ' + getdataparams_temp)
os.system('chmod +x ' + freqmeshfile_temp)
os.system('chmod +x ' + periapseshift_ecc_temp)
  
######################################################################
# define the output file name for this instance of the script

output_file = "%s/injection_%.3f-%.3f_%.3e_%.3e_%d.data" %(outdir,float(fmin),float(fmax),float(h0),float(ecc),(int)(id))

# also define null result string - consists of 30 zeros
null_result = "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n"

######################################################################
# extract some dataset parameters

# define some filenames
p_stampsfile = tempworkarea + '/primary_stamps.data'
s_stampsfile = tempworkarea + '/secondary_stamps.data'
p_dataparamsfile = tempworkarea + '/primary_datasetparams.data'
s_dataparamsfile = tempworkarea + '/secondary_datasetparams.data'

try: os.system("%s --datadir %s --tsft %s --stampsfile %s --outfile %s" %(getdataparams_temp,pdatadir,tsft,p_stampsfile,p_dataparamsfile))
except:
  print "ERROR : failed to extract primary dataset parameters from %s !!!" %(pdatadir)
  sys.exit(1)

try: os.system("%s --datadir %s --tsft %s --stampsfile %s --outfile %s" %(getdataparams_temp,sdatadir,tsft,s_stampsfile,s_dataparamsfile))
except:
  print "ERROR : failed to extract primary dataset parameters from %s !!!" %(sdatadir)
  sys.exit(1)

# interpret the outputs (just need observation times at present)
p_dp = ConfigParser.ConfigParser()
p_dp.read(p_dataparamsfile)
p_tobs = p_dp.get('dataparams','tobs')
p_tspan = p_dp.get('dataparams','tspan')
p_start = p_dp.get('dataparams','start')
s_dp = ConfigParser.ConfigParser()
s_dp.read(s_dataparamsfile)
s_tobs = s_dp.get('dataparams','tobs')
s_tspan = s_dp.get('dataparams','tspan')
s_start = s_dp.get('dataparams','start')

######################################################################
# Make a temporary freqmeshfile

# define the freqmeshfile name
freqmeshfile = tempworkarea + '/temp_freqmeshfile.data'

# open a file for writing
try: fp = open(freqmeshfile, 'w')
except:
  print "ERROR : failed to open freq/mesh file %s" %(filename)
  sys.exit(1)

freqmesh_fband = float(fmax) - float(fmin)
fp.write("%.3f %.3f %.3f %s %s" %(float(fmin),float(fmax),freqmesh_fband,ptemplatefile,stemplatefile))  
  
fp.close()

######################################################################
# setup input template file

# define output filename
makefakedata_templatefile = "mfd_template.input"

# we set eccentricity here equal to ecc
try: os.system("%s --sourcefile %s --source %s --ifo %s --stamps %s --tsft %s --ecc %f --reftime %s --noisedir %s --ephem %s --yr %s --smaxis 1.0 --outfile %s" %(binaryinput_temp,sourcefile,source,det_primary,p_stampsfile,tsft,ecc,p_start,pdatadir,ephdir,yr,makefakedata_templatefile))
except:
  print "ERROR : failed to generate input template file !!!"
  sys.exit(1)

os.system('cat ' + makefakedata_templatefile)

######################################################################
# start the loop over injections

# but first lets define the start time for use later in generating seeds
tzero = float(time.time())
#print "tzero is %6.12f" %(tzero)

i=0
while i < (int)(ntrials):

######################################################################
# need to deal with eccentricity, argument of periapse and time of periapse passage

# as we are recovering signals using a circular orbit template we need
# to deal with signal parameters that describe eccentric orbits.  So we
# need to alter the time of peripase passage if we are selecting a
# random value for the argument of periapse.

######################################################################
# generate random input template file for each detector 

  # define some inputs
  # here we set the seed by the clock, we have chosen 100th second accuracy
  tnow = float(time.time())
  #print "tnow is %6.12f" %(tzero)
  #seed = int(1e2*(tnow -tzero))
  seed = (int)(id)*(int)(ntrials)+(int)(i)
  
  #print "seed is now %d" %(seed)
  p_makefakedata_randomfile = "primary_mfd.random"
  s_makefakedata_randomfile = "secondary_mfd.random"
  p_datasetbase = p_datasetdir + '/SFT'
  s_datasetbase = s_datasetdir + '/SFT'

  try: os.system("%s --ptemplatefile %s --stemplatefile %s --phi --psi --cosiota --h0MIN %s --h0MAX %s --f0MIN %s --f0MAX %s --argpMIN %s --argpMAX %s --seed %s --co --pstampsfile %s --sstampsfile %s --infile %s --poutfile %s --soutfile %s --fband %s --pnoisedir %s --snoisedir %s --psftbase %s --ssftbase %s" %(binaryrandominput_temp,ptemplatefile,stemplatefile,h0,h0,fmin,fmax,argpmin,argpmax,seed,p_stampsfile,s_stampsfile,makefakedata_templatefile,p_makefakedata_randomfile,s_makefakedata_randomfile,fband,pdatadir,sdatadir,p_datasetbase,s_datasetbase))
  except:
    print "ERROR : failed to generate random input files for making fake data !!!"
    sys.exit(1)

  #os.system('cat ' + p_makefakedata_randomfile)
  #os.system('cat ' + s_makefakedata_randomfile)
  
######################################################################
# Extract the random parameters from the data

  p_dp = ConfigParser.ConfigParser()
  p_dp.read("randparams_python.data")
  f_inj = float(p_dp.get('randomparameters','f0'))
  sma_inj = float(p_dp.get('randomparameters','sma'))
  tperi_sec_inj = (int)(p_dp.get('randomparameters','tpsec'))
  tperi_nan_inj = (int)(p_dp.get('randomparameters','tpnano'))
  argp_inj = float(p_dp.get('randomparameters','argp'))
  period_inj = float(p_dp.get('randomparameters','period'))

  #os.system('cat randparams_python.data')
  print "we use the parameters %f %f %d %d %f %d" %(f_inj,sma_inj,tperi_sec_inj,tperi_nan_inj,argp_inj,period_inj)

######################################################################
# now alter the time of periapse passage based on value of argument of periapse

  try: os.system("%s --tperisec %d --tperinano %d --period %f --ecc %f --argp %f" %(periapseshift_ecc_temp,tperi_sec_inj,tperi_nan_inj,period_inj,ecc,argp_inj))
  except:
    print "ERROR : failed to calculate new non-circular orbit periapse passage time !!!"
    sys.exit(1)

  os.system("cat tperi_ecc.data")  

  # read in the parameters from file
  tpn_dp = ConfigParser.ConfigParser()
  tpn_dp.read("tperi_ecc.data")
  tperisec_new_inj = (int)(tpn_dp.get('tperi_ecc_new','tperi_new_sec'))
  tperinano_new_inj = (int)(tpn_dp.get('tperi_ecc_new','tperi_new_nano'))
  tempfile1 = tempworkarea + "/tempfile1"
  tempfile2 = tempworkarea + "/tempfile2"
  #print "tperi_new_inj is %d %d" %(tperisec_new_inj,tperinano_new_inj)

  #os.system("echo 's/^\(orbitTperiSSBsec.*=\) [0-9]+/\1 %d/' %s" %(tperisec_new_inj,p_makefakedata_randomfile))
  os.system("touch %s %s" %(tempfile1,tempfile2))
  #os.system("cat %s" %(p_makefakedata_randomfile))

  # and alter the value of Tp in the primary  mfd input files  
  try: os.system("sed 's/^\\(orbitTperiSSBsec.*=\\) [0-9][0-9]*/\\1 %d/' %s > %s" %(tperisec_new_inj,p_makefakedata_randomfile,tempfile1))
  except:
    print "ERROR : failed to replace primary random periapse time with new eccentric time !!!"
    sys.exit(1)
  try: os.system("sed 's/^\\(orbitTperiSSBns.*=\\) [0-9][0-9]*/\\1 %d/' %s > %s" %(tperinano_new_inj,tempfile1,tempfile2))
  except:
    print "ERROR : failed to replace primary random periapse time with new eccentric time !!!"
    sys.exit(1)
  try: shutil.copyfile(tempfile2,p_makefakedata_randomfile) 
  except:
    print "ERROR : failed to copy temporary file to primary makefakedata input file !!!"
    sys.exit(1)

  #os.system("cat %s" %(tempfile1))
  #os.system("cat %s" %(tempfile2))
  #os.system("cat %s" %(p_makefakedata_randomfile))
  #sys.exit(1)

  # and alter the value of Tp in the primary  mfd input files  
  try: os.system("sed 's/^\\(orbitTperiSSBsec.*=\\) [0-9][0-9]*/\\1 %d/' %s > %s" %(tperisec_new_inj,s_makefakedata_randomfile,tempfile1))
  except:
    print "ERROR : failed to replace secondary random periapse time with new eccentric time !!!"
    sys.exit(1)
  try: os.system("sed 's/^\\(orbitTperiSSBns.*=\\) [0-9][0-9]*/\\1 %d/' %s > %s" %(tperinano_new_inj,tempfile1,tempfile2))
  except:
    print "ERROR : failed to replace secondary random periapse time with new eccentric time !!!"
    sys.exit(1)
  try: shutil.copyfile(tempfile2,s_makefakedata_randomfile) 
  except:
    print "ERROR : failed to copy temporary file to secondary makefakedata input file !!!"
    sys.exit(1)    

  

  #sed 's/^\(orbitT.*=\) [0-9]+/\1 '${var}'/'  FILE > tmp
  #sed 's/^\(orbitT.*=\) [0-9][0-9]*/\1 '${var}'/'  FILE > tmp
  #mv tmp FILE


######################################################################
# make the data for each detector

  try: os.system("%s @%s" %(makefakedata_temp,p_makefakedata_randomfile))
  except:
    print "ERROR : failed to generate fake data for primary detector !!!"
    sys.exit(1)

  try: os.system("%s @%s" %(makefakedata_temp,s_makefakedata_randomfile))
  except:
    print "ERROR : failed to generate fake data for secondary detector !!!"
    sys.exit(1)

######################################################################
# make the submesh for each detector

  # define some file names
  p_subbank = tempworkarea + '/primary_subbank.data'
  s_subbank = tempworkarea + '/secondary_subbank.data'

  try: os.system("%s --sma %s --tpsec %s --tpnano %s --fullbank %s --subbank %s --ephem %s --yr %s --Nsub 1" %(makesubmesh_temp,sma_inj,tperi_sec_inj,tperi_nan_inj,ptemplatefile,p_subbank,ephdir,yr))
  except:
    print "ERROR : failed to generate primary sub-mesh !!!"
    sys.exit(1)
  
  try: os.system("%s --sma %s --tpsec %s --tpnano %s --fullbank %s --subbank %s --ephem %s --yr %s --Nsub 1" %(makesubmesh_temp,sma_inj,tperi_sec_inj,tperi_nan_inj,stemplatefile,s_subbank,ephdir,yr))
  except:
    print "ERROR : failed to generate secondary sub-mesh !!!"
    sys.exit(1)

  #os.system('tail -1 ' + p_subbank)
  #os.system('tail -1 ' + s_subbank)

######################################################################
# determine the frequency bands for searching

  # define primary search frequency and band
  p_searchres = float(1.0/(float(overres)*float(p_tspan)))
  p_fres = float(1.0/float(p_tspan))
  p_fmin = float(p_searchres*math.floor(((float(f_inj)-(float(co_nbins)*float(p_fres)))/float(p_searchres))+0.5))
  p_fband = 2.0*float(co_nbins)*float(p_fres)
  p_fmax = p_fmin + p_fband

  #print "p_searchres is %f p_fres is %f" %(p_searchres,p_fres)
  #print "p_fmin = %6.12f p_fmax = %6.12f p_fband = %6.12f\n" %(p_fmin,p_fmax,p_fband)
  
  # define secondary search frequency and band
  s_searchres = (1.0/(float(overres)*float(s_tspan)))
  s_fres = (1.0/float(s_tspan))
  s_fmin = float(s_searchres*math.floor(((float(f_inj)-(float(co_nbins)*float(s_fres)))/float(s_searchres))+0.5))
  s_fband = 2.0*float(co_nbins)*float(s_fres)
  s_fmax = s_fmin + s_fband

  #print "s_searchres is %f s_fres is %f" %(s_searchres,s_fres)
  #print "s_fmin = %6.12f s_fmax = %6.12f s_fband = %6.12f\n" %(s_fmin,s_fmax,s_fband)

  # define the minimum and maximum frequencies for the coincidence analysis
  co_fmin = s_fmin
  co_fmax = s_fmin + s_fband

  #print "co_fmin is %6.12f co_fmax is %6.12f" %(co_fmin,co_fmax)

######################################################################
# do the searches for each detector

  # define some output variables
  p_outputlabel = "_%s_%.3f-%.3f.data" %(det_primary,p_fmin,p_fmax)
  s_outputlabel = "_%s_%.3f-%.3f.data" %(det_secondary,s_fmin,s_fmax)

  # testing
  #p_Fout = "/raid/1/cm/p_Fout_%d_%d_%e_%f.data" %((int)(i),(int)(id),float(h0),float(argp_inj))
  #s_Fout = "/raid/1/cm/s_Fout_%d_%d_%e_%f.data" %((int)(i),(int)(id),float(h0),float(argp_inj))

  try: os.system("%s --dterms %s --Freq %s --FreqBand %s --overres %s --DataDir %s --EphemDir %s --EphemYear %s --IFO %s --binary --dopplermax %s --Fthreshold %s --windowsize %s --outputLabel %s --binarytemplatefile %s --sourcefile %s --source %s --workingDir %s" %(search_temp,dterms,p_fmin,p_fband,overres,p_datasetdir,ephdir,yr,det_primary,dopplermax,thresh_primary,windowsize,p_outputlabel,p_subbank,sourcefile,source,p_searchresultsdir))
  except:
    print "ERROR : failed to complete search on primary injection !!!"
    sys.exit(1)

  #os.system('cat ' + p_searchresultsdir + '/*')
  #print "done primary search"  

  try: os.system("%s --dterms %s --Freq %s --FreqBand %s --overres %s --DataDir %s --EphemDir %s --EphemYear %s --IFO %s --binary --dopplermax %s --Fthreshold %s --windowsize %s --outputLabel %s --binarytemplatefile %s --sourcefile %s --source %s --workingDir %s" %(search_temp,dterms,s_fmin,s_fband,overres,s_datasetdir,ephdir,yr,det_secondary,dopplermax,thresh_secondary,windowsize,s_outputlabel,s_subbank,sourcefile,source,s_searchresultsdir))
  except:
    print "ERROR : failed to complete search on secondary injection !!!"
    sys.exit(1)

  #os.system('cat ' + s_searchresultsdir + '/*')
  #print "done secondary search"

######################################################################
# do the coincidence analysis

  try: os.system("%s --fmin %s --fmax %s --nbins %s --sdataparamsfile %s --presultsdir %s --sresultsdir %s --coresultsdir %s --freqmeshfile %s --ephdir %s --yr %s" %(coincidence_temp,fmin,fmax,co_nbins,s_dataparamsfile,p_searchresultsdir,s_searchresultsdir,co_resultsdir,freqmeshfile,ephdir,yr))
  except:
    print "ERROR : failed to do coincidence analysis !!!"
    sys.exit(1)

  #print "done coincidence"

######################################################################
# store the results

  # the results should be stored in the coincidence directory in a single file
  # there should be a single result but there could be two so we need to deal with this
  # so we select the loudest one

  # extract the filelist
  filelist = os.listdir(co_resultsdir)

  # open the main results file
  try: fpo = open(output_file, 'a')
  except:
    print "ERROR : failed to open injections results file %s" %(output_file)
    sys.exit(1)

  #print "length of filelist is %d" %(len(filelist))   
 
  # store the result only if we have one
  if (len(filelist)>0):

    #print "found the co result"
    filename=co_resultsdir + '/' + filelist[0]

    #print "opening file %s" %(filename)
    
    # open the coincidence results file
    try: fp = open(filename, 'r')
    except:
      print "ERROR : strange, can't open coincidence results file %s" %(filename)
      sys.exit(1)

    fp.flush() 
    #os.system("cat %s" %(filename))  
   
    # loop over each line in this file
    loudest = 999999999.0
    flag=0
    for line in fp:
      temp = line.rstrip().split()
      co_sig = float(temp[22])
      #print "loudest is %f co_sig is %f" %(float(loudest),float(co_sig))
      if (co_sig < loudest):
        loudest = co_sig
        loudest_line = line
        flag=1

    fp.close
          
    if (flag==1): 
      fpo.write((str)(f_inj) + ' ' + loudest_line)
      #print "RESULT %d = %f %s" %(i,f_inj,loudest_line)
    else:
      fpo.write((str)(f_inj) + null_result)
      #print "RESULT %d = %f %s" %(i,f_inj,null_result)

  else:
    fpo.write((str)(f_inj) + null_result)
    #print "RESULT %d = %f %s" %(i,f_inj,null_result)

  fpo.close

######################################################################
# now we need to tidy up after ourselves within the loop (BEING CAREFUL)

  if os.path.exists(p_datasetdir):
    try: shutil.rmtree(p_datasetdir)
    except:
      print "ERROR : failed to delete primary fake data set in injections !!!"
      sys.exit(1)
    try: os.mkdir(p_datasetdir)
    except:
      print "ERROR : failed to create primary fake data set directory in injections !!!"
      sys.exit(1)
    
  if os.path.exists(s_datasetdir):
    try: shutil.rmtree(s_datasetdir)
    except:
      print "ERROR : failed to delete secondary fake data set in injections !!!"
      sys.exit(1)
    try: os.mkdir(s_datasetdir)
    except:
      print "ERROR : failed to create secondary fake data directory set in injections !!!"
      sys.exit(1)  

  if os.path.exists(p_searchresultsdir):
    try: shutil.rmtree(p_searchresultsdir)
    except:
      print "ERROR : failed to delete primary results in injections !!!"
      sys.exit(1)
    try: os.mkdir(p_searchresultsdir)
    except:
      print "ERROR : failed to create primary results directory set in injections !!!"
      sys.exit(1)    

  if os.path.exists(s_searchresultsdir):
    try: shutil.rmtree(s_searchresultsdir)
    except:
      print "ERROR : failed to delete secondary results in injections !!!"
      sys.exit(1)
    try: os.mkdir(s_searchresultsdir)
    except:
      print "ERROR : failed to create secondary results directory in injections !!!"
      sys.exit(1) 

  if os.path.exists(co_resultsdir):
    try: shutil.rmtree(co_resultsdir)
    except:
      print "ERROR : failed to delete coincidence results in injections !!!"
      sys.exit(1)
    try: os.mkdir(co_resultsdir)
    except:
      print "ERROR : failed to create coincidence results directory in injections !!!"
      sys.exit(1)  

  if os.path.exists(p_makefakedata_randomfile):
    try: os.remove(p_makefakedata_randomfile)
    except:
      print "ERROR : failed to delete primary random input file in injections !!!"
      sys.exit(1)

  if os.path.exists(s_makefakedata_randomfile):
    try: os.remove(s_makefakedata_randomfile)
    except:
      print "ERROR : failed to delete secondary random input file in injections !!!"
      sys.exit(1)

  if os.path.exists(p_subbank):
    try: os.remove(p_subbank)
    except:
      print "ERROR : failed to delete primary submesh file in injections !!!"
      sys.exit(1)

  if os.path.exists(s_subbank):
    try: os.remove(s_subbank)
    except:
      print "ERROR : failed to delete secondary submesh file in injections !!!"
      sys.exit(1)

  if os.path.exists("randparams_python.data"):
    try: os.remove("randparams_python.data")
    except:
      print "ERROR : failed to delete python random parameter config file in injections !!!"
      sys.exit(1)

  if os.path.exists(tempfile1):
    try: os.remove(tempfile1)
    except:
      print "ERROR : failed to delete temporary file in injections !!!"
      sys.exit(1)

  if os.path.exists(tempfile2):
    try: os.remove(tempfile2)
    except:
      print "ERROR : failed to delete temporary file in injections !!!"
      sys.exit(1)

  if os.path.exists("tperi_ecc.data"):
    try: os.remove("tperi_ecc.data")
    except:
      print "ERROR : failed to delete periapse shift outfile in injections !!!"
      sys.exit(1)      

######################################################################
# increment the loop and end it
  i=i+1
  #print "i is %d and ntrials is %d" %(i,(int)(ntrials))

######################################################################
# now we need to tidy up after ourselves (BEING CAREFULL)

# change back to original directory
try: os.chdir(pwd)
except:
  print "ERROR : failed to change back to original working directory %s !!!" %(pwd)
  sys.exit(1)

# now delete the temporary working directory if it exists 
if os.path.exists(tempworkarea):
  try: shutil.rmtree(tempworkarea)
  except:
    print "ERROR : failed to delete temporary working directory %s !!!" %(tempworkarea)
    sys.exit(1)
    
######################################################################
# finished

sys.exit(0)

