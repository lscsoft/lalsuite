"""
Running the omega pipeline with condor for an arbitrary amount of sources
"""

__author__ = "Jeroen Meidam"
__credits__ = ["Jeroen Meidam"]
__maintainer__ = "Jeroen Meidam"
__email__ = "jeroen.meidam@ligo.org"
__status__ = ""

usage="""  omegascans_dag.py config.ini [options]
  This script creates a DAG to run the omega pipeline for multiple detectors
  and multiple sources on a cluster.
  Framefiles are automatically looked up with gw_data_find

  It requires a config file and a sourcefile containing trigtimes
  and detector timeslides.

  Run with --example to create an example config and sourcefile.

  It is possible to include auxiliary channels in the scans only for H1 and L1.
  For Virgo this is only possible on one of the Virgo clusters
"""

###############################################################################
#
# LOAD LIBRARIES
#
###############################################################################

from lal import pipeline
from optparse import OptionParser
import uuid
import os
import sys
import ast
from subprocess import Popen,PIPE

from configparser import ConfigParser

###############################################################################
#
# DEFINITIONS
#
###############################################################################

DefaultFrameTypes = {'L1':'L1_RDS_R_L1',#Includes auxiliary channels
                     'H1':'H1_RDS_R_L1',#Includes auxiliary channels
                     'V1':'T1300121_V1_EARLY_RECOLORED_V2'}

DefaultConfigFiles = {'L1':'/archive/home/omega/configurations/S6/S6b/L0L1-RDS_R_L1-selected.txt',
                      'H1':'/archive/home/omega/configurations/S6/S6b/H0H1-RDS_R_L1-selected.txt',
                      'V1':''}

RecoloredFrameTypes = {'L1':'T1200307_V4_EARLY_RECOLORED_V2',
                       'H1':'T1200307_V4_EARLY_RECOLORED_V2',
                       'V1':'T1300121_V1_EARLY_RECOLORED_V2'}

ConfigRecoloredV = \
"""
# QScan configuration file

# Scans H1 and L1 gravitational-wave channel data
# from the S5 version 1 calibrated data set

# Shourov K. Chatterji
# shourov@ligo.caltech.edu
# 2006-02-11

[Context,Context]

[Parameters,Parameter Estimation]

[Notes,Notes]

[Gravitational,Gravitational wave data]

{
  channelName:                 'V1:h_16384Hz'
  frameType:                   'T1300121_V1_EARLY_RECOLORED_V2'
  sampleFrequency:             4096
  searchTimeRange:             128
  searchFrequencyRange:        [32 Inf]
  searchQRange:                [3.32 141]
  searchMaximumEnergyLoss:     0.2
  whiteNoiseFalseRate:         1e-3
  alwaysPlotFlag:               1
  searchWindowDuration:        0.5
  plotTimeRanges:              [1 8 100]
  plotFrequencyRange:          []
  plotNormalizedEnergyRange:   [0 25.5]
}
"""

ConfigRecoloredL = \
"""
# QScan configuration file

# Scans H1 and L1 gravitational-wave channel data
# from the S5 version 1 calibrated data set

# Shourov K. Chatterji
# shourov@ligo.caltech.edu
# 2006-02-11

[Context,Context]

[Parameters,Parameter Estimation]

[Notes,Notes]

[Gravitational,Gravitational wave data]

{
  channelName:                 'L1:LDAS-STRAIN'
  frameType:                   'T1200307_V4_EARLY_RECOLORED_V2'
  sampleFrequency:             4096
  searchTimeRange:             128
  searchFrequencyRange:        [32 Inf]
  searchQRange:                [3.32 141]
  searchMaximumEnergyLoss:     0.2
  whiteNoiseFalseRate:         1e-3
  alwaysPlotFlag:               1
  searchWindowDuration:        0.5
  plotTimeRanges:              [1 8 100]
  plotFrequencyRange:          []
  plotNormalizedEnergyRange:   [0 25.5]
}
"""

ConfigRecoloredH = \
"""
# QScan configuration file

# Scans H1 and L1 gravitational-wave channel data
# from the S5 version 1 calibrated data set

# Shourov K. Chatterji
# shourov@ligo.caltech.edu
# 2006-02-11

[Context,Context]

[Parameters,Parameter Estimation]

[Notes,Notes]

[Gravitational,Gravitational wave data]

{
  channelName:                 'H1:LDAS-STRAIN'
  frameType:                   'T1200307_V4_EARLY_RECOLORED_V2'
  sampleFrequency:             4096
  searchTimeRange:             128
  searchFrequencyRange:        [32 Inf]
  searchQRange:                [3.32 141]
  searchMaximumEnergyLoss:     0.2
  whiteNoiseFalseRate:         1e-3
  alwaysPlotFlag:               1
  searchWindowDuration:        0.5
  plotTimeRanges:              [1 8 100]
  plotFrequencyRange:          []
  plotNormalizedEnergyRange:   [0 25.5]
}
"""

ConfigsRecolored = {'L1':ConfigRecoloredL,
                    'H1':ConfigRecoloredH,
                    'V1':ConfigRecoloredV}

ExampleConfig = \
"""
[paths]
#The basedir is where the dag and sub files end up.
basedir=/home/jmeidam/tiger_runs/omegascans/example

#Within the outdir you will have one folder for each trigtime
#  and within that, a resultsfolder (the resultspage) for each ifo
out=/home/jmeidam/public_html/omegascans/example

#A file that contains a list of sources, each of which has at least
#  a trigtime, and timeslide for each detector (can be 0.0)
sourcefile=/home/jmeidam/tiger_runs/omegascans/example/omegascanslist.txt



[omegapipe]
executable=/home/omega/opt/omega/bin/wpipeline



[analysis]
ifos=['L1','H1','V1']

frametypes={'L1':'L1_RDS_R_L1','H1':'H1_RDS_R_L1','V1':'T1300121_V1_EARLY_RECOLORED_V2'}

accounting_group=ligo.dev.o1.cbc.testgr.tiger

#oddslimit can be empty if all sources in sourcefile are to be scanned
#  in this example, only odds > 20 are scanned
oddslimit={'low':20.0}

#comment out "basicscan" to include auxiliary channels as well
basicscan=

#leave out to use defaults, also ignored when basicscan is set
#conigfiles={'L1':'conig_L1.txt','H1':'conig_H1.txt','V1':'conig_V1.txt'}



[cbcwikitable]
#A table will be written out to the output folder that can be pasted into
#a CBC wiki page. The table includes links to all the scans.

#When outlier is given a value, the CBCwiki_table will only include entries
#  where LogOdds > outlier. Not needed if already set with "oddslimit"
outlier=20.0

#the webdir is used to create links to the omega scan pages.
#  Make sure it corresponds to the "out" directory in the [paths] section
webdir=https://ldas-jobs.ligo.caltech.edu/~jmeidam/omegascans/example
"""

ExampleSourceFile = \
"""trigtime H1 L1 V1 logodds
966944343 0 1869348 593939 107.23
970687764 0 -1963504 -903565 156.14
969100904 0 -109089 -2545560 96.04
"""

###############################################################################
#
# FUNCTIONS
#
###############################################################################

def mkdirs(path):
  """
  Helper function. Make the given directory, creating intermediate
  dirs if necessary, and don't complain about it already existing.
  """
  if os.access(path,os.W_OK) and os.path.isdir(path): return
  else: os.makedirs(path)

def system_call(command):
    p = Popen([command], stdout=PIPE,stderr=PIPE, shell=True)
    out, err = p.communicate()
    return out, err#p.stdout.read(),p.stderr.read()

def get_framedir(trigtime,timeslide,ifo,frametype):
  """
  Use gw_data_find to find the frame directory name which is used as
  an argument to wpipeline.
  """
  scantime = trigtime+timeslide

  #get framefile location with datafind:
  command="/usr/bin/gw_data_find --observatory=%s --url-type=file --type=%s --gps-start-time=%.3f --gps-end-time=%.3f"%(ifo[0],frametype,scantime,scantime)
  datafind_stdout,datafind_stderr = system_call(command)
  if not datafind_stdout:
    print(datafind_stderr)
    exit("gw_data_find failed, exiting...")
  gwf_file = datafind_stdout.replace("file://localhost","").strip()
  gwf_dir = os.path.dirname(gwf_file).strip()

  #An extra check to be sure the output was not gibberish
  if not os.path.isfile(gwf_file):
    exit("%s is not a file or does not exist, exiting..."%gwf_file)

  return gwf_dir


def fix_subfile(subfile):
  """
  The omega pipeline executable is followed by "scan", which
  needs to be in front of all the other arguments and options.
  This function corrects the final subfile to this end.
    - It removes "scan" from the executable
    - Places it as the first item in "arguments"
  """
  with open(subfile, 'r') as f:
    lines = f.readlines()
  i=0
  i_exec = 0
  i_args = 0
  for l in lines:
    if l.split()[0] == "executable":
      i_exec=i
    if l.split()[0] == "arguments":
      i_args=i
    i+=1

  #remove last item of the executable line, which should be "scan".
  if lines[i_exec].strip().split()[-1] == "scan":
    execitem = lines[i_exec].strip().split()[:-1]
    lines[i_exec] = ' '.join(execitem)
    lines[i_exec] += '\n'

  #add scan as the first argument
  x = lines[i_args]
  lines[i_args] = x.replace('"','" scan',1)

  #rewrite subfile
  with open(subfile, 'w') as f:
    for line in lines:
      f.write(line)


def fix_scriptfile(dagpath,dagfilename,executable):
  """
  If the scan argument is not included in the executable,
  it will not show up in the sh file.
  This function corrects that.
  """
  scriptname = ".".join(dagfilename.split(".")[:-1]) + ".sh"
  with open(os.path.join(dagpath,scriptname), 'r') as f:
    lines = f.readlines()

  if len(executable.split())>1:
    if executable.split()[1] == 'scan':
      print("scan found")
      return

  #do the following if scan is not present
  i=0
  for l in lines:
    spl = l.strip().split()
    if len(spl) > 1:
      if len(spl[0].split('/')) > 1: #the line we need
        x = spl[0]
        x += " scan"
        for item in spl[1:]:
          x+=" %s"%item
        x+='\n'
        lines[i]=x
    i+=1

  #rewrite
  with open(os.path.join(dagpath,scriptname), 'w') as f:
    for line in lines:
      f.write(line)

###############################################################################
#
# CLASSES
#
###############################################################################

class OmegaScansDAG(pipeline.CondorDAG):
  def __init__(self,cp):
    self.config=cp

    ## Setup some paths and filenames
    self.basepath = self.config.get('paths','basedir')
    self.sourcefile = self.config.get('paths','sourcefile')
    self.outpath = self.config.get('paths','out')
    mkdirs(self.outpath)
    self.daglogfile=os.path.join(self.basepath,'omegascans-'+str(uuid.uuid1())+'.log')
    pipeline.CondorDAG.__init__(self,self.daglogfile)
    if cp.has_option('paths','logdir'):
      self.logpath=cp.get('paths','logdir')
    else:
      self.logpath=os.path.join(self.basepath,'log')
    mkdirs(self.logpath)
    self.dagfilename="omegascans"
    self.submitFile = os.path.join(self.basepath,'omegascans.sub')

    #get ifos
    if cp.has_option('analysis','ifos'):
      self.ifos=ast.literal_eval(cp.get('analysis','ifos'))
    else:
      self.ifos=['H1','L1','V1']

    #Set the DAG's filename
    self.set_dag_file(os.path.join(self.basepath,self.dagfilename))

    #Set the config files and frametypes
    self.frametypes=ast.literal_eval(cp.get('analysis','frametypes'))
    if cp.has_option('analysis','configfiles'):
      self.configfiles=ast.literal_eval(cp.get('analysis','configfiles'))
    else:
      if cp.has_option('analysis','basicscan'):
        self.configfiles={}
        self.frametypes=RecoloredFrameTypes
        for ifo in self.ifos:
          self.configfiles[ifo]=os.path.join(self.basepath,"basic_config_%s.txt"%ifo)
          with open(os.path.join(self.basepath,self.configfiles[ifo]),'w') as f:
            f.write(ConfigsRecolored[ifo])
      else:
        if 'V1' in self.ifos:
          #For V1 there is no other file available on non-virgo clusters
          DefaultConfigFiles['V1']=os.path.join(self.basepath,"basic_config_V1.txt")
          with open(os.path.join(self.basepath,DefaultConfigFiles['V1']),'w') as f:
            f.write(ConfigsRecolored['V1'])
        self.configfiles=DefaultConfigFiles


    #get info (e.g. trigtimes and timeslides) from data file
    info = self.ReadInfoFromFile(self.ifos)
    Nsources = len(info['trigtime'])

    #create the job and add the nodes
    job = OmegaScanJob(self.config,self.submitFile,self.logpath)

    #set limits for analysis
    self.oddslimit_set = False
    if self.config.has_option('analysis','oddslimit'):
      self.oddslimit=ast.literal_eval(self.config.get('analysis','oddslimit'))
      if self.oddslimit.has_key('low') or self.oddslimit.has_key('high'):
        self.oddslimit_set = True
      if not self.oddslimit.has_key('low'):
        self.oddslimit['low'] = float('-inf')
      if not self.oddslimit.has_key('high'):
        self.oddslimit['high'] = float('inf')
    else:
      self.oddslimit={'low':float('-inf'),'high':float('inf')}



    #This is where we store data for a wiki table
    self.table_entries=[]

    print("calling gw_data_find for each node...")
    for n in range(Nsources):
      timeslides={}
      if self.oddslimit_set:
        if info['logodds'][n] > self.oddslimit['low'] and info['logodds'][n] < self.oddslimit['high']:
          for ifo in self.ifos:
            node = OmegaScanNode(info['trigtime'][n],info[ifo][n],ifo,self.frametypes[ifo],self.configfiles[ifo],self.outpath,job)
            timeslides[ifo] = info[ifo][n]
            self.add_node(node)
          self.add_table_entry(info['trigtime'][n],timeslides,self.ifos,info['logodds'][n])
      else:
        for ifo in self.ifos:
          node = OmegaScanNode(info['trigtime'][n],info[ifo][n],ifo,self.frametypes[ifo],self.configfiles[ifo],self.outpath,job)
          timeslides[ifo] = info[ifo][n]
          self.add_node(node)
        self.add_table_entry(info['trigtime'][n],timeslides,self.ifos,0.0)
    print("success!")

    self.write_table()

  def add_table_entry(self,trigtime,timeslides,ifos,logodds):
    table_entry = {}
    for ifo in ifos:
      table_entry[ifo]=timeslides[ifo]
    table_entry['trigtime']=trigtime
    table_entry['logodds']=logodds

    self.table_entries.append(table_entry)

  def write_table(self):
    """
    writes a cbc wiki type table with some information and links
    """

    filename = os.path.join(self.outpath,"CBCwiki_table.txt")
    print("writing CBC wiki table to \"%s\""%filename)

    outlier=None
    if self.config.has_option('cbcwikitable','outlier'):
      outlier=float(self.config.get('cbcwikitable','outlier'))

    fp = open(filename,'w')

    webdir=None
    if self.config.has_option('cbcwikitable','webdir'):
      webdir=self.config.get('cbcwikitable','webdir')

    header="||'''trigtime'''"
    for ifo in self.ifos:
      header+="||'''timeslide %s'''"%ifo
    for ifo in self.ifos:
      header+="||'''injtime %s'''"%ifo
    header+="||'''log odds'''||"

    fp.write(header+"\n")

    entries_to_write=[]
    if outlier:
      for e in self.table_entries:
        if float(e['logodds']) > outlier:
          entries_to_write.append(e)
    else:
      entries_to_write = self.table_entries

    for e in entries_to_write:
      string="||%d"%e['trigtime']
      #timeslide entries
      for ifo in self.ifos:
        string+="||%d"%e[ifo]
      #injection time entries
      for ifo in self.ifos:
        time = float(e['trigtime'])+float(e[ifo])
        if webdir:
          link = os.path.join(webdir,"scan_%d"%e['trigtime'],"%s_%.2f"%(ifo,time))
          string+="||[[%s|%d]]"%(link,time)
        else:
          string+="||%d"%time
      string+="||%.3f||"%e['logodds']
      fp.write(string+"\n")

    fp.close()


  def ReadInfoFromFile(self,ifos):
    """
    Reads in a file containing columns for at least trigtime, timeslide ifo1, timeslide ifo2 and timeslide ifo3 resp.
    """
    headerdict={}
    info = {}

    oddslimit_set = False
    if self.config.has_option('analysis','oddslimit'):
      oddstest = ast.literal_eval(self.config.get('analysis','oddslimit'))
      if oddstest.has_key('low') or oddstest.has_key('high'):
        oddslimit_set = True



    with open(self.sourcefile,'r') as f:

      data=f.readlines()

      #Read the header to get the IFO names
      header = data[0].strip().split()
      Nsources = len(data)-1

      #redefine data without the header
      data = data[1:]

      #which header label corresponds to which collum
      for i in range(len(header)):
        headerdict[header[i]] = i

      for ifo in ifos:
        if not headerdict.has_key("timeslide_"+ifo) and not headerdict.has_key(ifo):
          sys.exit("ERROR: Not all ifos from config file are present in \""+self.sourcefile+"\"")

      if not headerdict.has_key('trigtime') and not headerdict.has_key('injtime'):
        sys.exit("ERROR: The \""+self.sourcefile+"\" does not contain \"trigtime\" column.")

      if oddslimit_set:
        #values for logodds are required if a limit was set
        if not headerdict.has_key('logodds') and not headerdict.has_key('logO'):
          sys.exit("ERROR: An odds limit was set, but \""+self.sourcefile+"\" does not contain \"logodds\" column.")

      #change header dict to use preferred format
      for ifo in ifos:
        if headerdict.has_key("timeslide_"+ifo):
          val = headerdict["timeslide_"+ifo]
          headerdict[ifo]=val
      if headerdict.has_key("injtime"):
        val = headerdict["injtime"]
        headerdict["trigtime"]=val
      if oddslimit_set:
        if headerdict.has_key("logO"):
          val = headerdict["logO"]
          headerdict["logodds"]=val

      #read in the remaining lines
      for ifo in ifos:
        timeslist=[]
        for n in range(Nsources):
          linesplit = data[n].strip().split()
          col = headerdict[ifo]
          timeslist.append(float(linesplit[col]))
        info[ifo]=timeslist

      #read the trigtimes
      trigtimeslist=[]
      if oddslimit_set: logoddslist=[]
      for n in range(Nsources):
        linesplit = data[n].strip().split()
        coltrig = headerdict['trigtime']
        trigtimeslist.append(float(linesplit[coltrig]))
        if oddslimit_set:
          colodds = headerdict['logodds']
          logoddslist.append(float(linesplit[colodds]))

      info['trigtime'] = trigtimeslist
      if oddslimit_set:
        info['logodds'] = logoddslist

    return info

class OmegaScanJob(pipeline.CondorDAGJob):
  def __init__(self,cp,submitFile,logdir):
    self.basepath=cp.get('paths','basedir')

    pipeline.CondorDAGJob.__init__(self,'vanilla',cp.get('omegapipe','executable'))

    self.set_sub_file(os.path.abspath(submitFile))
    self.machine_count=str(1)
    self.machine_memory=str(1024) # default value
    self.add_condor_cmd('accounting_group',cp.get('analysis','accounting_group'))
    self.add_condor_cmd('RequestMemory',str(2000))
    self.add_arg("--report")
    self.set_stdout_file(os.path.join(logdir,'omegascans-$(cluster)-$(process)-$(node).out'))
    self.set_stderr_file(os.path.join(logdir,'omegascans-$(cluster)-$(process)-$(node).err'))

class OmegaScanNode(pipeline.CondorDAGNode):

  def __init__(self,trigtime,timeslide,ifo,frametype,configfile,outdir,job,logodds=None):
    self.framedir=get_framedir(trigtime,timeslide,ifo,frametype)
    pipeline.CondorDAGNode.__init__(self,job)
    self.__finalized=False
    scantime = trigtime+timeslide

    #sometimes there are unwanted lockfiles lying around.
    #they will be removed, but the user is warned
    lockfile = os.path.join(outdir,"scan_%.0f"%trigtime,"%s_%.2f"%(ifo,scantime),"lock.txt")
    if os.path.isfile(lockfile):
      print("WARNING: lock file found in output directory.\n         Deleting it now, but check that you did not forget to stop some ongoing run.")
      system_call('rm %s'%lockfile)

    self.add_var_arg('%.3f'%scantime)
    self.add_var_opt("outdir",os.path.join(outdir,"scan_%.0f"%trigtime,"%s_%.2f"%(ifo,scantime)))
    self.add_var_opt("configuration",configfile)
    self.add_var_opt("framecache",self.framedir)

  def finalize(self):
    if self.__finalized:
      return
    self.__finalized=True





def main():

  #############################################################################
  #
  # ARGUMENT PARSING
  #
  #############################################################################

  parser=OptionParser(usage)
  parser.add_option("-e","--example",default=False,dest="example",action="store_true",help="Create example config.ini and an example sourcefile")
  (opts,args) = parser.parse_args()

  if opts.example:
    with open("omega_config.ini","w") as f:
      f.write(ExampleConfig)
    with open("omegascanslist.txt","w") as f:
      f.write(ExampleSourceFile)

    print("Example files \"omega_config.ini\" and \"omegascanslist.txt\" are created")
    sys.exit(0)

  if len(args) != 1:
    parser.print_help()
    sys.exit("ERROR: Must provide one config.ini")

  cp=ConfigParser()
  cp.optionxform = str
  try:
    cp.read_file(open(args[0]))
  except AttributeError:
    cp.readfp(open(args[0]))

  dag=OmegaScansDAG(cp)

  dag.write_sub_files()
  dag.write_dag()
  dag.write_script()

  #fix the sub and sh files
  #This is required because pipeline.py does not yet have the ability to add
  #a specific argument before all other arguments and options ('scan' in this case)
  fix_subfile(dag.submitFile)
  fix_scriptfile(cp.get('paths','basedir'),dag.get_dag_file(),cp.get('omegapipe','executable'))

  print('Successfully created DAG file.')
  fulldagpath=os.path.join(cp.get('paths','basedir'),dag.get_dag_file())
  print('Now run condor_submit_dag %s\n'%(fulldagpath))




###############################################################################
#
# START THE MAIN FUNCTION
#
###############################################################################

if __name__ == "__main__":
	# START THE MAIN FUNCTION IF RUN AS A SCRIPT. OTHERWISE, JUST LOADS THE CLASS
	# AND FUNCTION DEFINITIONS
	exit(main())
