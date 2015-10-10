from optparse import OptionParser
import argparse
import io
import os
import ast
import sys
import getpass
import ConfigParser
from lalapps import inspiralutils

user=getpass.getuser()

## 2013 Salvatore Vitale, Michalis Agathos. Python to setup both injection and init files for TIGER runs. 
## Will create the injection, the configs, the folders, and call the pipeline exec file.
## 2014 User can pass pre-existing xml table.
## 2014 Will look for science and veto segments and generate times for unvetoed injections.
## 2014 All options moved to configuration file.

################################################################################
#
#   DEFINE CONSTANTS
#
################################################################################

usage='''%prog [options] config.ini
Setup a batch of runs for TIGER and invoke lalinference_multi_pipe based on the config.ini file provided as command-line argument.

A pre-existing injection file may be used with the -I option, otherwise one will be generated automatically.
The base directory for standard output and the output directory for post-processing should also be provided using the -O and -P options respectively.

Additionally, a local directory for log output may be specified with the -L option and a scratch directory with the -S option.
'''

'''
MGparams_approx_dic = {
"TaylorF2Test":["dphi0","dphi1", "dphi2", "dphi3", "dphi4", "dphi5", "dphi5l", "dphi6", "dphi6l", "dphi7"],
"TaylorT4Test":["dchi0","dchi1", "dchi2", "dchi3", "dchi4", "dchi5", "dchi5l", "dchi6", "dchi6l", "dchi7"],
"SpinTaylorT4Test":["dchi0","dchi1", "dchi2", "dchi3", "dchi4", "dchi5", "dchi5l", "dchi6", "dchi6l", "dchi7"]
}
'''

################################################################################
#
#   DEFINE FUNCTIONS
#
################################################################################

def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = range(r)
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)

def createCombinationsFromList(list):
        outputlist = []
        outputlist.append("GR")
        for L in xrange(len(list)):
                for subset in combinations(list, L+1):
                        """
                        temp = ""
                        for subsetitem in subset:
                                temp += subsetitem
                        """
                        outputlist.append(subset)

        return outputlist

def write_pipe_init(dic, cp):
    #cp_temp = cp_ConfigParser(cp)
    cp_temp = cp #FIXME: deep copy
    #FIXME: cp_temp.read_dict(dic) instead!
    for sec in dic.keys():
        for opt in dic[sec].keys():
            val = dic[sec][opt]
            if not type(val)==str:
                val=str(val)
            cp_temp.set(sec, opt, val)

    with open("pipeline.ini", "w") as ofile:
        cp_temp.write(ofile)

def ensure_dir(f):
    if not os.path.isdir(f):
        os.makedirs(f)

# Deep copy ConfigParser instance
def cp_ConfigParser(config):
    config_string = io.StringIO()
    config.write(config_string)
    # We must reset the buffer ready for reading.
    config_string.seek(0) 
    new_config = ConfigParser.RawConfigParser()
    new_config.read_file(config_string)
    return new_config

################################################################################
#
#   PARSE COMMAND LINE OPTIONS
#
################################################################################

#parser = OptionParser()
parser = argparse.ArgumentParser(usage)
parser.add_argument('config', metavar='CONFIG_FILE', type=str, nargs='+', help='A (list of) configuration file(s) containing sections and options for the entire pipeline')
parser.add_argument("-O", dest='basedir', type=str, help="Path to base directory")
parser.add_argument("-P", dest='postproc', type=str, help="Path to post-processing output (optional)", default=None)
parser.add_argument("-I", dest='injfile', type=str, help="Path to a pre-existing injection .xml file (optional)", default=None)
parser.add_argument("-L", dest='logdir', type=str, help="Path to log directory (optional)", default=None)
parser.add_argument("-S", dest='scratchdir', type=str, help="Path to scratch directory (optional)", default=None)
parser.add_argument("-g",'--gid', dest='gid', help="run from a graceID id (optional)", default=None)
parser.add_argument("--condor-submit",action="store_true",default=False,help="Automatically submit the condor dag")
args = parser.parse_args()

config_file = args.config

if args.basedir is not None:
    basefolder = os.path.abspath(args.basedir)
else:
    print 'Error: No base directory provided! Exiting...'
    sys.exit(1)

if args.postproc is not None:
    webdir = os.path.abspath(args.postproc)
else:
    webdir = None

if args.injfile is not None:
    injfile = os.path.abspath(args.injfile)
    if os.path.exists(injfile):
        print 'TIGER: Reading injections from file'
    else:
        print 'Error: Cannot find xml file for injections: ' + injfile
        sys.exit(1)
else:
    injfile = None

if args.logdir is not None:
    logdir = os.path.abspath(args.logdir)
else:
    logdir = None
    
if args.scratchdir is not None:
    scratchdir = os.path.abspath(args.scratchdir)
else:
    scratchdir = None

################################################################################
#
#   READ TIGER VARIABLES FROM CONFIG FILE
#
################################################################################

cp = ConfigParser.RawConfigParser()
cp.optionxform = str
# FIXME: bring all defaults here
cp.read(config_file)

if not cp.has_section('tiger'):
  print 'Invalid configuration file! No "tiger" section found.'
  sys.exit(1)


dic_engine = {}

# Your inspinj seed. The inspnest dataseed will be created from this, adding three zeros at the end (e.g. inspinj 7001 --> inspnest 7001000)
inspinj_seed = cp.getint('tiger', 'seed') 
dataseed = 1000*inspinj_seed
if cp.has_option('analysis', 'dataseed'):
    dataseed = cp.getint('analysis', 'dataseed')
dic_analysis={"dataseed":dataseed,}

# A descriptive name of what is being run. Will be the name of the postprocessing result folder and a part of the injection file, for the records
tiger_tag = cp.get('tiger', 'tiger-tag') 

# Setup output directory names
basefolder = os.path.join(basefolder, tiger_tag, str(inspinj_seed))
ensure_dir(basefolder)

if webdir is None:
    webdir = cp.get('paths', 'webdir')
webdir = os.path.join(webdir, tiger_tag, str(inspinj_seed)) 
baseurl = os.path.join(cp.get('paths', 'baseurl'), tiger_tag, str(inspinj_seed))

# This has to be either GR or MG for modified gravity, or NO for analyzing pure data
type_inj = cp.get('tiger','type-inj') 
gpstimefile = None
gid=None
# The injection approximant string, and PN order, e.g. inj_approx=TaylorF2  inj_pnorder=threePointFivePN
inj_approx = cp.get('tiger', 'inj-approx') 
inj_PNorder = cp.get('tiger', 'inj-pnorder')

if cp.has_option('tiger', 'dense-inj'):
    whereinj='often'
    inj_every=cp.getint('tiger', 'inj-every')
else:
    whereinj='middle'

#remote_script='svitale@login.nikhef.nl:/project/gravwav/safe_append.sh' ## This is the remote file which appends to the database
#remote_database='TF2Cut_db.txt'   ## Succesful runs are appended to remote_database. Failed runs are appended to 'remote_database'_failed

# This is the number of signals created in the xml file. Inspnest will analize all of them.
num_events = cp.getint('tiger', 'num-events') 
if type_inj!='NO':
  # GPS start time before the 1st injection
  sta_time = cp.getint('input','gps-start-time') 
  # GPS end time after the last injection
  end_time = cp.getint('input','gps-end-time') 

timeslides=False
# Check if timeslides should be used for injections
if cp.has_option('input', 'timeslides'):
  timeslides = cp.getboolean('input', 'timeslides')

# List of IFO identifiers
ifos = ast.literal_eval(cp.get('analysis','ifos'))


# Parse options for modified gravity injections

if type_inj == 'MG':
  # Read list of modified gravity parameters.
  mgpars = ast.literal_eval(cp.get('tiger','mg-params'))

  # Read modified parameter shifts. If type_inj is 'GR' this will be ignored.
  mgshifts = ast.literal_eval(cp.get('tiger','mg-shifts'))      

  # Read distribution list for the shifts. (default: constant for all)
  # Available options: 'const' for constant shift, 'uniform' for uniform from 0 to mg-shifts or 'gauss' for normal distributed as N(mg-shifts, mg-sigmas)
  if cp.has_option('tiger', 'mg-distr'):
      mgdistr = cp.get('tiger', 'mg-distr')
  if len(mgpars) is not len(mgshifts):
    print 'Error: Number of MG shift values does not match number of MG parameters'
    sys.exit(1)

  # Currently inspinj only supports constant modGR parameters.
  if mgdistr != 'const':
      print "Error: Only mg-distr='const' is available at the moment"
      sys.exit(1)

  # Read stdevs for normal-distributed deviations. (default: 0 for all)
  if mgdistr=='gauss':
    if cp.has_option('tiger', 'mg-sigmas'):
      mgsigmas = ast.literal_eval(cp.get('tiger', 'mg-sigmas'))
      if len(mgpars) is not len(mgsigmas):
          print 'Error: Number of MG stdev values does not match number of MG parameters'
          sys.exit(1)
    else:
      print 'Error: Gaussian distribution requested but no mg-sigmas provided'
      sys.exit(1)
elif type_inj == 'NO':
    if not cp.has_option('input','gps-time-file') and not args.gid:
        print "Error: TIGER called without injections but no gps-time-file provided"
        sys.exit(1)
    if cp.has_option('input','gps-time-file'):
      gpstimefile = cp.get('input','gps-time-file')
    elif args.gid is not None:
      gid=args.gid 
# FIXME: What calibration options are used and where? ([calibration]?)
#CALIB_SEED = cp.getint('calibration', 'calib-seed')

hypotheses = ast.literal_eval(cp.get('tiger', 'test-params'))


#Define defaults
dic_inj={
#"min_snr":8,
#"max_snr":25, # NetSNR
#"start_freq":30, # f_low for SNR calculation
#"coinc_snr":5.5,
#"min_m":1,
#"max_m":3,
"seed":inspinj_seed,
"waveform":inj_approx+inj_PNorder,
}

#Read rest of dict from config file

distrange = ast.literal_eval(cp.get('tiger', 'inj-dist-range'))

mass1range = ast.literal_eval(cp.get('tiger', 'inj-m1-range'))
mass2range = ast.literal_eval(cp.get('tiger', 'inj-m2-range'))
dic_inj.update({
"min_m1":mass1range[0], "max_m1":mass1range[1],
"min_m2":mass2range[0], "max_m2":mass2range[1],
"min_mtot":mass1range[0]+mass2range[0], "max_mtot":mass1range[1]+mass2range[1],
"min-distance":distrange[0], "max-distance":distrange[1],
})

if cp.has_option('tiger', 'inj-spins'):
    inject_spins = cp.get('tiger', 'inj-spins')
    spin1range = ast.literal_eval(cp.get('tiger', 'inj-a1-range'))
    spin2range = ast.literal_eval(cp.get('tiger', 'inj-a2-range'))
    dic_inj.update({"enable-spin":"", "min-spin1":spin1range[0], "max-spin1":spin1range[1], "min-spin2":spin2range[0], "max-spin2":spin2range[1]})
    if cp.has_option('tiger','inj-a-gaussian'):
        spin1mean = cp.get('tiger','inj-a1-mean')
        spin2mean = cp.get('tiger','inj-a2-mean')
        spin1stdev = cp.get('tiger','inj-a1-stdev')
        spin2stdev = cp.get('tiger','inj-a2-stdev')
        dic_inj.update({"mean-spin1":spin1mean, "stdev-spin1":spin1stdev, "mean-spin2":spin2mean, "stdev-spin2":spin2stdev})
    if inject_spins=='aligned':
        dic_inj.update({"aligned":"",})
else:
    dic_inj.update({"disable-spin":"",})


if type_inj=="MG":
    dic_inj.update({'enable-dchi':''})
    for mgp,mgs in zip(mgpars,mgshifts):
        dic_inj.update({mgp:mgs})
        print 'modGR at ', mgp, ' by ', mgs

#uname=os.uname()
#if  any(["atlas" in i for i in uname]):
#    WWWdir="WWW/LSC"
#else:
#    WWWdir="public_html"


################################################################################
#
#   REAL DATA: FETCHING THE SEGMENT FILES & GENERATING INJECTION TIMES
#
################################################################################

seglen = 300
psdlen = 1024
if cp.has_option('engine','seglen'):
    seglen = cp.getint('engine','seglen')
if cp.has_option('input','max-psd-length'):
    psdlen = cp.getint('input','max-psd-length')

if injfile is None and gpstimefile is None and gid is None and not cp.has_option('lalinference', 'fake-cache'):
  from lalinference.tiger import make_injtimes
  print 'TIGER: Generating science and veto segment files for real data'
  if not (cp.has_option('input','gps-start-time') and cp.has_option('input','gps-end-time')):
    print "make_injtimes needs both gps start and end time"
    sys.exit(1)
    
  segfolder = os.path.join(basefolder, 'segments')
  ensure_dir(segfolder)
  curdir = os.getcwd()
  os.chdir(segfolder)

  vetoCats = ast.literal_eval(cp.get('segments', 'veto-categories'))

  vetoDefFile = inspiralutils.downloadVetoDefFile(cp, True)
  inspiralutils.generate_veto_cat_files(cp, vetoDefFile, True)

  # Populate dictionary with veto filenames for each IFO and generate the veto segment files
  veto_dic = {}
  seg_dic = {}


  # Download science segment files for each IFO
  for IFO in ifos:
    (seg_dic[IFO], veto_dic[IFO]) = inspiralutils.findSegmentsToAnalyze(cp, IFO, vetoCats, generate_segments=True, use_available_data=False, data_quality_vetoes=True)


  # Feed science and veto segment files to make_injtimes and generate injection times and timeslides
  print 'TIGER: Generating GPS times for unvetoed injections'
  ensure_dir(os.path.join(basefolder, 'injtimes'))
  IFOdict={}
  for ifo in ifos: 
    IFOdict[ifo]=make_injtimes.IFO(ifo, seg_dic[ifo], veto_dic[ifo][4], minlen=max(psdlen, seglen))
  os.chdir(curdir)


  # Generate injtimes file
  if timeslides:
    print 'TIGER: Injection times will be generated with timeslides'
    timesdict = {}
    # timeslidefile = os.path.join(basefolder, 'injtimes', 'slides_%s_%s_%s.dat'%(str(sta_time), str(end_time)))
    # Generate timeslides file
    for i in IFOdict.keys():
    # Generate triggers for each IFO (for timeslides only)
      timesdict[i] = IFOdict[i].getTrigTimes(whereInj=whereinj, interval=inj_every, lmargin=seglen)

    # Generate timeslide output
    injtimesfile, slidefile = make_injtimes.generateTimeslides(timesdict, num_events, ref=ifos[0], outfolder=os.path.join(basefolder,'injtimes'))
    #label=''.join(ifos) + '_' + str(num_events)
    #injtimesfile = os.path.join(basefolder,'injtimes','injtimes_%s.dat'%(label))
    # slidefile = os.path.join(basefolder,'injtimes','timeslides_%s.dat'%(label))
    cp.set('input', 'timeslide-file', slidefile)
  else:
    print 'TIGER: Injection times will be generated on coincident unvetoed time (no timeslides)'
    # Combine segments from differnt IFOs to get multi-IFO unvetoed injections
    compIFO = IFOdict.values()[0]
    for ifo in IFOdict.values()[1:]:
        compIFO = make_injtimes.getDoubles(compIFO, ifo)
    injtimesfile = os.path.join(basefolder, 'injtimes', 'injtimes_%s_%s.dat'%(compIFO._name, str(num_events)))
    compIFO.getTrigTimes(whereInj=whereinj, interval=inj_every, lmargin=seglen, n=num_events, outfile=injtimesfile)
  dic_inj.update({"t-distr":"file", "time-file":injtimesfile})

elif gpstimefile is None and gid is None:
  time_step=((end_time-2-sta_time-seglen)/(num_events-1))
  dic_inj.update({"time-step":time_step,})


################################################################################
#
#   CREATING THE INJECTION FILE
#
################################################################################

if injfile is None and gpstimefile is None and gid is None:

  print "TIGER: Creating the xml file\n"
  print tiger_tag
      
  #inspinjname=os.path.join(basefolder,'injections_%s_%s_SNR_%s_%s.xml'%(tiger_tag,dic_inj['seed'],dic_inj['min_snr'],dic_inj['max_snr']))
  inspinjname=os.path.join(basefolder,'injections_%s_%s.xml'%(tiger_tag,dic_inj['seed']))


  dic_inj.update({
  "lalapps_inspinj":cp.get('tiger', 'lalapps_inspinj'),
  "gps-start-time":sta_time+seglen,
  "gps-end-time":end_time-2,
  "output": inspinjname,
#  "min_tot_m":(dic_inj["min_m"]*2),
#  "max_tot_m":(dic_inj["max_m"]*2)
  })


  #--min-snr MIN_SNR --max-snr MAX_SNR --snr-distr volume  --ligo-fake-psd LALSimAdLIGO --virgo-fake-psd LALSimAdVirgo --ligo-start-freq START_FREQ --virgo-start-freq START_FREQ --ifos H1,L1,V1"
  string="LALAPPS_INSPINJ --f-lower 10.0 --gps-start-time GPS-START-TIME --gps-end-time GPS-END-TIME --seed SEED --waveform WAVEFORM --d-distr volume --l-distr random --i-distr uniform --min-mass2 MIN_M2 --max-mass2 MAX_M2 --min-mass1 MIN_M1 --max-mass1 MAX_M1 --m-distr componentMass --min-mtotal MIN_MTOT --max-mtotal MAX_MTOT --amp-order 0 --time-step TIME-STEP --output OUTPUT "

  for p in dic_inj.keys():
      if not p.upper() in string:
          if type(dic_inj[p])==str:
              string=string+ " --"+p +" "+dic_inj[p]
          else:
              string=string+ " --"+p +" "+repr(dic_inj[p])
      else:
          string=string.replace(p.upper()+" ", "%s "%(dic_inj[p]))
  print string+"\n"
  os.system(string)

elif injfile is not None:
  dic_inj.update({"output":injfile})




################################################################################
#
#   CREATE ALL POSSIBLE COMBINATIONS FROM A LIST
#
################################################################################

curdir = os.getcwd()
foldernames=""
parser_paths=""
allcombinations = createCombinationsFromList(hypotheses)

for run in allcombinations:

    # DETERMINE FOLDER NAME
    foldername = ""
    subhy=""
    if type(run) is not tuple:
        foldername=run
        subhy=run
    else:
        for item in run:
            foldername += item
            subhy+= item +","
        subhy=subhy[:-1]
    ensure_dir(os.path.join(basefolder,foldername))
    dic_path={"webdir":os.path.join(webdir, foldername),"baseurl":os.path.join(baseurl, foldername)}
    os.chdir(os.path.join(basefolder,foldername))
    foldernames+=foldername+' '
    parser_paths+=str(os.path.join(basefolder,foldername,"pipeline.ini"))+" "
    if subhy!='GR':
        dic_engine.update({'grtest-parameters':subhy})
    else:
        dic_engine.update({'grtest-parameters':''})
    dic = {'engine':dic_engine, 'analysis':dic_analysis, 'paths':dic_path}
    write_pipe_init(dic, cp)
    

#foldernames=foldernames[:-1]
if type_inj == "NO":
    pipestring="%s "%cp.get('tiger','lalinference_multi_pipe')
    if gid is not None:
      pipestring+="--gid %s"%gid
    condor_s=""
    if args.condor_submit:
      condor_s=" --condor-submit "
    pipestring+=condor_s

    pipestring+=" -F %s -r %s" %(foldernames, basefolder)
 
else:
    pipestring="%s -I %s -F %s -r %s" %(cp.get('tiger','lalinference_multi_pipe'), dic_inj['output'], foldernames, basefolder)

if logdir is not None:
    pipestring+= " -p %s "%logdir
#if scratchdir is not None:
#    pipestring+= " -l %s "%scratchdir
pipestring+=" %s "%parser_paths

print pipestring
os.system(pipestring)

# RETURN TO CURRENT WORKING DIRECTORY
os.chdir(curdir)
