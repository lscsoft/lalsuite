"""
plot_inspiral_hipe.in - generates plots from ihope output

This script generates a condor DAG to make plots that summarize the
results of an inspiral search
"""

__author__ = 'Patrick Brady <patrick@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

##############################################################################
# import standard modules
import sys, os, copy, math
import socket, time
import re, string
from optparse import *
import tempfile
from glue.pipeline import DeepCopyableConfigParser as dcConfigParser
import urlparse
from pylal import CoincInspiralUtils
from glue import iterutils
##############################################################################
# import the modules we need to build the pipeline
from glue import pipeline

import inspiral
from inspiralutils import determine_sieve_patterns
import inspiralutils
##############################
# Convert ifo list to string
def combo2str(combo):
  ifo_list = ""
  for ifo in combo:
    ifo_list += ifo
  return ifo_list

def set_common_options(cp, this_section):
  """
  """
  common_options = ["cache-file", "gps-start-time", "gps-end-time", "enable-output", "output-path"]
  for common_option in common_options:
    try:
      cp.set(this_section, common_option,  string.strip(cp.get('common',common_option)))
    except:
      print "warning:", common_option, "is strongly recommended in the [common] section of your ini file"
      pass

  return cp

def addVetoFileForPlotinspiral( job, veto_file ):
   """
   Add another veto_file to the option "veto-file" for plotinspiral
   to veto the on-source segment in the case of the external trigger search.
   Might insert the following function into pipeline.py to make
   this more sophisticated:
   """

   # see if there is already an value set for this option
   value = job.get_opt( "veto-file" )

   # then add this veto_file or set a new option
   if value:
     job.add_opt( "veto-file", value+','+veto_file )
   else:
     job.add_opt( "veto-file", veto_file )


#############################################################################
#
#  MAIN PROGRAM
#
##############################################################################
usage = """usage: %prog [options] 
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")
    
parser.add_option("-u", "--user-tag",action="store",type="string",\
    default=None,metavar=" USERTAG",\
    help="tag the jobs with USERTAG (overrides value in ini file)")

parser.add_option("-G", "--g1-data",action="store_true",default=False,\
    help="analyze g1 data")

parser.add_option("-H", "--h1-data",action="store_true",default=False,\
    help="analyze h1 data")

parser.add_option("-B", "--h2-data",action="store_true",default=False,\
    help="analyze h2 data")

parser.add_option("-L", "--l1-data",action="store_true",default=False,\
    help="analyze l1 data")

parser.add_option("-V", "--v1-data",action="store_true",default=False,\
    help="analyze v1 data")

parser.add_option("-S", "--one-ifo",action="store_true",default=False,\
    help="analyze single ifo data (not usable for GEO)")

parser.add_option("-D", "--two-ifo",action="store_true",default=False,\
    help="analyze two interferometer data")

parser.add_option("-T", "--three-ifo",action="store_true",default=False,\
    help="analyze three interferometer data")

parser.add_option("-Q", "--four-ifo",action="store_true",default=False,\
    help="analyze four interferometer data")
  
parser.add_option("-R", "--five-ifo",action="store_true",default=False,\
    help="analyze five interferometer data")

parser.add_option("-a", "--analyze-all",action="store_true",default=False,\
    help="analyze all combinations of IFOs (overrides --n-ifo options)")

parser.add_option("-i", "--plotinspiral" ,action="store_true",default=False,\
    help="run plotinspiral to summarize filtering output")

parser.add_option("-t", "--plotthinca" ,action="store_true",default=False,\
    help="run plotthinca to summarize coincidence")

parser.add_option("-n", "--plotnumtemplates" ,action="store_true",default=False,\
    help="run plotnumtemplates to plot trigbanks and templtbanks")

parser.add_option("-e", "--plotethinca" ,action="store_true",default=False,\
    help="run plotethinca parameter for a pair of ifo combos")

parser.add_option("-M", "--plotinspmissed" ,action="store_true",default=False,\
    help="run plotinspmissed to plot missed triggers")

parser.add_option("-j", "--plotinspinj" ,action="store_true",default=False,\
    help="run plotinspinj to plot inspiral injections.")

parser.add_option("-s", "--plotsnrchi" ,action="store_true",default=False,\
    help="run plotsnrchi to plot  snr vs chisq for a glob of triggers.")

parser.add_option("-z", "--plotinspiralrange" ,action="store_true",default=False,\
    help="run plotinspiralrange to plot range of inspiral (Mpc) vs end time (in days).")

parser.add_option("", "--ploteffdistcut" ,action="store_true",default=False,\
    help="run ploteffdistcut to plot the effective distance cut.")

parser.add_option("", "--plotgrbtimeslidestats" ,action="store_true",default=False,\
    help="run pylal_grbtimeslide_stats to plot statistics of the GRB background.")

parser.add_option("-f", "--config-file",action="store",type="string",\
    metavar=" FILE",help="use configuration file FILE")

parser.add_option("-p", "--log-path",action="store",type="string",\
    metavar=" PATH",help="directory to write condor log file")
parser.add_option("", "--first-stage",action="store_true", default=False,\
    metavar="FIRSTSTAGE",help="to make plots for output of the first inspiral stage")
parser.add_option("", "--second-stage",action="store_true", default=False,\
    metavar="SECONDSTAGE",help="to make plots for output of the second inspiral stage")
parser.add_option("-P", "--priority",action="store",type="int",\
    metavar=" PRIO",help="run jobs with condor priority PRIO")

#################################

parser.add_option("-w", "--write-script", action="store_true", default=False,
      help="write the workflow to a locally executable script")



command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#################################
#ALPHABETS USED IN THIS CODE.
# a, e, f, i, j, n, p, s, t, u, v, w, z; B, D, G, H, L, M, P, Q, R, S, T, V

#################################
# if --version flagged
if opts.version:
  print "$Id$"
  sys.exit(0)

#################################
# Sanity check of input arguments
################################

# Checks for config file
if not opts.config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use --config-file FILE to specify location."
  sys.exit(1)

# Checks for log path
if not opts.log_path:
  print >> sys.stderr, "No log file path specified."
  print >> sys.stderr, "Use --log-path PATH to specify a location."
  sys.exit(1)

# Checks for at least one ifo is specified  
if not opts.g1_data and not opts.h1_data and not opts.h2_data and \
    not opts.l1_data and not opts.v1_data:
  print >> sys.stderr, "No ifos specified.  Please specify at least one of"
  print >> sys.stderr, "--g1-data, --h1-data, --h2-data, --l1-data, --v1-data"
  print >> sys.stderr, "or use --analyze-all to analyze all ifos"
  sys.exit(1)

# Checks for H1 and H2 data when running ploteffdistcut
if opts.ploteffdistcut and not (opts.h1_data and opts.h2_data):
  print >> sys.stderr, "How can I plot effdistcut when H1 and H2 are not"
  print >> sys.stderr, "being analysed??? I will not run ploteffdistcut."
  opts.ploteffdistcut = False

# Checks for plots which require more than two ifos
if ((opts.plotthinca or opts.plotethinca or opts.plotinspmissed or \
    opts.plotinspinj or opts.plotsnrchi or opts.ploteffdistcut) and \
    not (opts.two_ifo or opts.three_ifo or opts.four_ifo or opts.analyze_all)):
  print >> sys.stderr, "No number of ifos given. Please specify at least one of"
  print >> sys.stderr, "--two-ifo, --three-ifo, --four-ifo"
  print >> sys.stderr, "or use --analyze-all to analyze all ifos"
  sys.exit(1)

##############################################################################
# create the config parser object and read in the ini file
cp = dcConfigParser()
cp.read(opts.config_file)
 
# Checks if the directory as given by output path exists or not. If it exists, it skips else it creates the directory.
this_section="common"
config = set_common_options(cp, this_section)
dir = config.get('common','output-path') 

if not os.path.isdir(dir):
  os.mkdir(dir)

# Checks if input-user-tag is specified.
tag = cp.get('pipeline','input-user-tag')
if tag=="":
  input_user_tag = None
else: 
  input_user_tag = tag + '_'

# Copies the number of slides to the codes that use it
num_slides = cp.get('pipeline','num-slides')
if num_slides:
  for section in ['plotthinca', 'plotethinca', 'ploteffdistcut']:
    cp.set(section,'num-slides',num_slides)

###############################
# Construct list based on ifos supplied at command line

ifolist = [ifo for ifo in ('G1','H1', 'H2', 'L1', 'V1') \
            if getattr(opts, "%s_data" % ifo.lower())]

if opts.two_ifo:
   ifo_combo=list(iterutils.choices(ifolist,2))
if opts.three_ifo:
   ifo_combo=list(iterutils.choices(ifolist,2)) + list(iterutils.choices(ifolist,3))
if opts.four_ifo:
   ifo_combo=list(iterutils.choices(ifolist,2)) + list(iterutils.choices(ifolist,3)) + \
             list(iterutils.choices(ifolist,4))
if opts.analyze_all: 
   ifo_combo=CoincInspiralUtils.get_ifo_combos(ifolist)


##############################################################################
# try to make a directory to store the cache files and job logs
try: os.mkdir('logs')
except: pass

##############################################################################
# create the config parser object and read in the ini file
#cp = ConfigParser.ConfigParser()
#cp.read(opts.config_file)

##############################################################################
# if a usertag has been specified, override the config file
if opts.user_tag is not None:
  usertag = opts.user_tag
  cp.set('pipeline','user-tag', usertag)
elif cp.has_option("pipeline", "user-tag"):
  usertag = string.strip(cp.get('pipeline','user-tag'))
else:
  usertag = None

##############################################################################
 # if exttrig flag is specified create trigger xml filename
if cp.has_option("pipeline", "exttrig"):

  # get the path to the GRB analysis directory
  exttrig_path = string.strip(cp.get('pipeline','exttrig'))

  # extract the usertag: its the last name in the path
  usertag =  exttrig_path.split('/')[-1]

  # create the external trigger filename and the onsource file
  exttrig_file = exttrig_path+'/trigger'+usertag+'.xml'
  exttrig_onsource = exttrig_path+'/onSourceSeg.txt'
else:
  exttrig_file = None
  exttrig_onsource = None
  
##############################################################################
# create a log file that the Condor jobs will write to
basename = re.sub(r'\.ini',r'',opts.config_file)
tempfile.tempdir = opts.log_path
if usertag is not None:
  tempfile.template = basename + '.' + usertag + '.dag.log.'
else:
  tempfile.template = basename + '.dag.log.'
logfile = tempfile.mktemp()
fh = open( logfile, "w" )
fh.close()

##############################################################################
# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile)
if usertag is not None:
  dag.set_dag_file(basename + '.' + usertag )
else:
  dag.set_dag_file(basename )

# set better submit file names than the default
if usertag is not None:
  subsuffix = '.' + usertag + '.sub'
else:
  subsuffix = '.sub'

##############################################################################
# create the Condor jobs that will be used in the DAG
all_jobs = []



# inspiral:
if opts.plotinspiral:
  
  """ Plotinspiral is looped over all single ifos. Plotinspiral 
      has now an option called opts.description, which sieves the 
      cache file by description "SIRE" and an ifo by ifo-type namely 
      "H1", "H2" or "L1". For a particular ifo file corresponding 
      all times are sieved.
  """
  plotinspiral_jobs = {}
  this_section = "plotinspiral"
  cp = set_common_options(cp, this_section)   
  cp_second = copy.deepcopy(cp)
  for ifo in ifolist:
    if opts.first_stage:
      if cp.has_option('plotinspiral','snr-chisq'):cp.remove_option('plotinspiral','snr-chisq')
      if cp.has_option('plotinspiral','log-snr-chisq'):cp.remove_option('plotinspiral','log-snr-chisq')
      if cp.has_option('plotinspiral','hist-chisq'):cp.remove_option('plotinspiral','hist-chisq')
      if cp.has_option('plotinspiral','cum-hist-snr-chi'):cp.remove_option('plotinspiral','cum-hist-snr-chi')
      if cp.has_option('plotinspiral','hist-snr-chi'):cp.remove_option('plotinspiral','hist-snr-chi')    
      job = inspiral.PlotInspiralJob(cp)
      job.add_opt("ifo-times", ifo)
      job.add_opt("ifo-tag", "FIRST_" + ifo)
      job.set_sub_file( basename + '.plotinspiral_FIRST_' + ifo + subsuffix )
      if exttrig_onsource: addVetoFileForPlotinspiral( job, exttrig_onsource)
      pattern_dict = inspiralutils.determine_sieve_patterns(cp, "plotinspiral",
        "FIRST", input_user_tag)
      for opt, val in pattern_dict.iteritems():
        val = val.replace("CAT_4","")
        val = val.replace("CAT_2","")
        val = val.replace("CAT_3","")
        job.add_opt(opt, val)
      plotinspiral_jobs[ifo + "FIRST"] = job
      all_jobs.extend(plotinspiral_jobs.values())
    if opts.second_stage:
      job = inspiral.PlotInspiralJob(cp_second)
      job.add_opt("ifo-times", ifo)
      job.add_opt("ifo-tag", "SECOND_" + ifo)
      job.set_sub_file( basename + '.plotinspiral_SECOND_' + ifo + subsuffix)
      if exttrig_onsource: addVetoFileForPlotinspiral( job, exttrig_onsource)
      pattern_dict = inspiralutils.determine_sieve_patterns(cp_second, "plotinspiral", "SECOND", input_user_tag)
      for opt, val in pattern_dict.iteritems():
        job.add_opt(opt, val)
      plotinspiral_jobs[ifo + "SECOND"] = job
      all_jobs.extend(plotinspiral_jobs.values())
 


#plotthinca
if opts.plotthinca:
  plotthinca_jobs = {}
  this_section = "plotthinca"
  cp = set_common_options(cp, this_section)
  for ifos in  ifo_combo:
    combostring = combo2str(ifos)
    if opts.first_stage:
      pattern_dict = determine_sieve_patterns(cp, "plotthinca", "FIRST",
        input_user_tag)
      plotthinca_jobs[combostring + "FIRST"] = inspiral.PlotThincaJob(cp)
      for ifo in ifos:
        plotthinca_jobs[combostring + "FIRST"].add_opt(ifo.lower() + "-triggers","  ")
        plotthinca_jobs[combostring + "FIRST"].add_opt("ifo-times", combostring)
        plotthinca_jobs[combostring + "FIRST"].add_opt("ifo-tag", "FIRST_" + combostring)
        plotthinca_jobs[combostring + "FIRST"].add_opt("statistic", "snr")
        for opt, val in pattern_dict.iteritems():
          val = val.replace("_CAT_4","")
          val = val.replace("_CAT_2","")
          val = val.replace("_CAT_3","")
          plotthinca_jobs[combostring + "FIRST"].add_opt(opt, val)
        plotthinca_jobs[combostring + "FIRST"].set_sub_file( basename + '.plotthinca_FIRST_' + combostring + subsuffix )
    if opts.second_stage:
      pattern_dict = determine_sieve_patterns(cp, "plotthinca",
        "SECOND_" + combostring, input_user_tag)
      plotthinca_jobs[combostring + "SECOND"] = inspiral.PlotThincaJob(cp)
      for ifo in ifos:
        plotthinca_jobs[combostring + "SECOND"].add_opt(ifo.lower() + "-triggers","  ")
        plotthinca_jobs[combostring + "SECOND"].add_opt("ifo-times", combostring)
        plotthinca_jobs[combostring + "SECOND"].add_opt("ifo-tag", "SECOND_" + combostring)
        plotthinca_jobs[combostring + "SECOND"].add_opt("statistic", "new_snr")
        for opt, val in pattern_dict.iteritems():
          plotthinca_jobs[combostring + "SECOND"].add_opt(opt, val)
        plotthinca_jobs[combostring + "SECOND"].set_sub_file( basename + '.plotthinca_SECOND_' + combostring + subsuffix )


  all_jobs.extend(plotthinca_jobs.values())  

#plotnumtemplates
if opts.plotnumtemplates:
  """ Plotnumtemplates needs no looping. This code now separates cache files
      either by description = "TMPLTBANK" or "TRIGBANK". In this way TMPLTBANKS 
      and TRIGBANKS for all ifos namely "H1", "H2" and "L1" are sieved together.
  """
  plotnumtemplates_jobs = {}
  this_section = "plotnumtemplates"
  cp = set_common_options(cp, this_section)
  plotnumtemplates_jobs = inspiral.PlotNumtemplatesJob(cp)
  plotnumtemplates_jobs.set_sub_file( basename + '.plotnumtemplates_' + subsuffix)
  pattern_dict = determine_sieve_patterns(cp, this_section, "", input_user_tag)
  for opt, val in pattern_dict.iteritems():
    print >> sys.stderr, opt, val
    if opt=="bank-pattern": val = "TMPLTBANK"
    plotnumtemplates_jobs.add_opt(opt, val)
  

  all_jobs.append(plotnumtemplates_jobs)


# plotethinca
if opts.plotethinca:

  plotethinca_jobs = {}
  this_section = "plotethinca"
  cp = set_common_options(cp, this_section)

  for ifos in ifo_combo:
    combostring = combo2str(ifos)
    ifosoption=ifos[0]
    for i in range(1,len(ifos)):
      ifosoption+=' --ifo '+ ifos[i]
      if opts.first_stage:
        plotethinca_jobs[combostring + "FIRST"] = inspiral.PlotEthincaJob(cp)
        plotethinca_jobs[combostring + "FIRST"].add_opt("ifo-times", combostring)
        plotethinca_jobs[combostring + "FIRST"].add_opt("ifo-tag", "FIRST_" + combostring)
        pattern_dict = determine_sieve_patterns(cp, "plotethinca",
                "FIRST", input_user_tag)
        for opt, val in pattern_dict.iteritems():
          val = val.replace("_CAT_4","")
          val = val.replace("_CAT_2","")
          val = val.replace("_CAT_3","")
          plotethinca_jobs[combostring + "FIRST"].add_opt(opt, val)
        plotethinca_jobs[combostring + "FIRST"].set_sub_file( basename + '.plotethinca_FIRST_' + combostring + subsuffix )
      if opts.second_stage:
        plotethinca_jobs[combostring + "SECOND"] = inspiral.PlotEthincaJob(cp)
        plotethinca_jobs[combostring + "SECOND"].add_opt("ifo-times", combostring)
        plotethinca_jobs[combostring + "SECOND"].add_opt("ifo-tag", "SECOND_" + combostring)
        pattern_dict = determine_sieve_patterns(cp, "plotethinca",
                "SECOND_" + combostring, input_user_tag)
        for opt, val in pattern_dict.iteritems():
          plotethinca_jobs[combostring + "SECOND"].add_opt(opt, val)
        plotethinca_jobs[combostring + "SECOND"].set_sub_file( basename + '.plotethinca_SECOND_' + combostring + subsuffix )

  all_jobs.extend(plotethinca_jobs.values())


#plotinspmissed

if opts.plotinspmissed:

  plotinspmissed_jobs = {}
  this_section = "plotinspmissed"
  cp = set_common_options(cp, this_section)
  for ifos in ifo_combo:
    combostring = combo2str(ifos)
    if opts.first_stage:
      if cp.has_option('plotinspmissed','do-followup'):
        cp.remove_option ('plotinspmissed','do-followup')
      plotinspmissed_jobs[combostring +"FIRST"] = inspiral.PlotInspmissedJob(cp)
      pattern_dict = determine_sieve_patterns(cp, this_section,
        "FIRST", input_user_tag)
      for opt, val in pattern_dict.iteritems():
        val = val.replace("CAT_4","")
        val = val.replace("CAT_2","")
        val = val.replace("CAT_3","")
        plotinspmissed_jobs[combostring + "FIRST"].add_opt(opt, val)
      plotinspmissed_jobs[combostring + "FIRST"].add_opt("ifo-times", combostring)
      plotinspmissed_jobs[combostring + "FIRST"].add_opt("ifo-tag", "FIRST_" + combostring)
      plotinspmissed_jobs[combostring +"FIRST"].set_sub_file( basename + '.plotinspmissed_FIRST_' + combostring + subsuffix )
    if opts.second_stage:
      plotinspmissed_jobs[combostring +"SECOND"] = inspiral.PlotInspmissedJob(cp)
      plotinspmissed_jobs[combostring + "SECOND"].add_opt("ifo-times", combostring)
      pattern_dict = determine_sieve_patterns(cp, this_section,
        "SECOND", input_user_tag)
      for opt, val in pattern_dict.iteritems():
        plotinspmissed_jobs[combostring + "SECOND"].add_opt(opt, val)
      plotinspmissed_jobs[combostring +"SECOND"].add_opt("ifo-tag", "SECOND_" + combostring)
      plotinspmissed_jobs[combostring +"SECOND"].set_sub_file( basename + '.plotinspmissed_SECOND_' + combostring + subsuffix )

  all_jobs.extend(plotinspmissed_jobs.values())


#ploteffdistcut
if opts.ploteffdistcut:

  ploteffdistcut_jobs = {}
  this_section = "ploteffdistcut"
  cp = set_common_options(cp, this_section)

  ifos = ('H1','H2')
  combostring = combo2str(ifos)
  if opts.first_stage:
    ploteffdistcut_jobs[combostring +"FIRST"] = inspiral.PlotEffdistcutJob(cp)
    pattern_dict = determine_sieve_patterns(cp, this_section,
      "FIRST", input_user_tag)
    for opt, val in pattern_dict.iteritems():
      val = val.replace("CAT_4","")
      val = val.replace("CAT_2","")
      val = val.replace("CAT_3","")
      ploteffdistcut_jobs[combostring + "FIRST"].add_opt(opt, val)
    ploteffdistcut_jobs[combostring + "FIRST"].add_opt( \
          "ifo-times", combostring)
    ploteffdistcut_jobs[combostring + "FIRST"].add_opt("ifoa", "H1")
    ploteffdistcut_jobs[combostring + "FIRST"].add_opt("ifob", "H2")
    ploteffdistcut_jobs[combostring + "FIRST"].add_opt("gps-start-time", 
          cp.get("common","gps-start-time"))
    ploteffdistcut_jobs[combostring + "FIRST"].add_opt("gps-end-time", 
          cp.get("common","gps-end-time"))
    ploteffdistcut_jobs[combostring + "FIRST"].add_opt("ifo-tag", 
          "FIRST_" + combostring)
    ploteffdistcut_jobs[combostring +"FIRST"].set_sub_file( basename + \
          '.ploteffdistcut_FIRST_' + combostring + subsuffix )
  if opts.second_stage:
    ploteffdistcut_jobs[combostring +"SECOND"] = inspiral.PlotEffdistcutJob(cp)
    ploteffdistcut_jobs[combostring + "SECOND"].add_opt("ifo-times", 
          combostring)
    pattern_dict = determine_sieve_patterns(cp, this_section,
      "SECOND", input_user_tag)
    for opt, val in pattern_dict.iteritems():
      ploteffdistcut_jobs[combostring + "SECOND"].add_opt(opt, val)
    ploteffdistcut_jobs[combostring + "SECOND"].add_opt("ifoa", "H1")
    ploteffdistcut_jobs[combostring + "SECOND"].add_opt("ifob", "H2")
    ploteffdistcut_jobs[combostring + "SECOND"].add_opt("gps-start-time",
          cp.get("common","gps-start-time"))
    ploteffdistcut_jobs[combostring + "SECOND"].add_opt("gps-end-time",
          cp.get("common","gps-end-time"))
    ploteffdistcut_jobs[combostring +"SECOND"].add_opt("ifo-tag", 
          "SECOND_" + combostring)
    ploteffdistcut_jobs[combostring +"SECOND"].set_sub_file( basename + '.ploteffdistcut_SECOND_' + combostring + subsuffix )

  all_jobs.extend(ploteffdistcut_jobs.values())


#plotinspinj
if opts.plotinspinj:
  plotinspinj_jobs = {}
  this_section = "plotinspinj"
  cp = set_common_options(cp, this_section)

  for ifo in ifolist:
    if opts.first_stage:
      plotinspinj_jobs[ifo + "FIRST"] = inspiral.PlotInspinjJob(cp)
      pattern_dict = determine_sieve_patterns(cp, this_section, "FIRST", 
          input_user_tag)
      for opt, val in pattern_dict.iteritems():
        val = val.replace("CAT_4","")
        val = val.replace("CAT_2","")
        val = val.replace("CAT_3","")
        plotinspinj_jobs[ifo + "FIRST"].add_opt(opt, val)
      plotinspinj_jobs[ifo + "FIRST"].add_opt("ifo-times", ifo)
      plotinspinj_jobs[ifo + "FIRST"].add_opt("ifo-tag", "FIRST_" + ifo)
      plotinspinj_jobs[ifo + "FIRST"].add_opt("title-text",ifo)
      plotinspinj_jobs[ifo + "FIRST"].set_sub_file( basename + '.plotinspinj_FIRST_' + ifo + subsuffix )
    if opts.second_stage:
      plotinspinj_jobs[ifo + "SECOND"] = inspiral.PlotInspinjJob(cp)
      pattern_dict = determine_sieve_patterns(cp, this_section, "SECOND", 
          input_user_tag)
      for opt, val in pattern_dict.iteritems():
        plotinspinj_jobs[ifo + "SECOND"].add_opt(opt, val)
      plotinspinj_jobs[ifo + "SECOND"].add_opt("ifo-times", ifo)
      plotinspinj_jobs[ifo + "SECOND"].add_opt("ifo-tag", "SECOND_" + ifo)
      plotinspinj_jobs[ifo + "SECOND"].add_opt("title-text",ifo)
      plotinspinj_jobs[ifo + "SECOND"].set_sub_file( basename + '.plotinspinj_SECOND_' + ifo + subsuffix )

  all_jobs.extend(plotinspinj_jobs.values())

#plotsnrchi
if opts.plotsnrchi:
  plotsnrchi_jobs = {}
  this_section = "plotsnrchi"
  cp = set_common_options(cp, this_section)

  if opts.first_stage:
    print >>sys.stderr, "warning: plotsnrchi can only run with second stage inputs.  Skipping..."
  if opts.second_stage:
    for ifo in ifolist:   
      plotsnrchi_jobs[ifo] = inspiral.PlotSnrchiJob(cp)
      plotsnrchi_jobs[ifo].add_opt("ifo-times", ifo)
      pattern_dict = determine_sieve_patterns(cp, this_section, "", 
          input_user_tag)
      for opt, val in pattern_dict.iteritems():
        plotsnrchi_jobs[ifo].add_opt(opt, val)
      plotsnrchi_jobs[ifo].set_sub_file( basename + '.plotsnrchi_SECOND_' + ifo + subsuffix )

  all_jobs.extend(plotsnrchi_jobs.values())

#plotinspiralrange
if opts.plotinspiralrange:
  """ Plotnumtemplates needs no looping. This code now separates cache files
      either by description = "TMPLTBANK" or "INSPIRAL". In this way TMPLTBANKS
      and TRIGBANKS for all ifos namely "H1", "H2" and "L1" are sieved together.
  """
  plotinspiralrange_jobs = {}
  this_section = "plotinspiralrange"
  cp = set_common_options(cp, this_section)
  
  pattern_dict = determine_sieve_patterns(cp, this_section, "", input_user_tag)
  plotinspiralrange_jobs["FIRST"] = inspiral.PlotInspiralrangeJob(cp)
  plotinspiralrange_jobs["FIRST"].set_sub_file( basename + '.plotinspiralrange_' + subsuffix)
  for opt, val in pattern_dict.iteritems():
    if opt=="bank-pattern":
      val = "TMPLTBANK"
    else:
      val = val.replace("_*_","_")
    plotinspiralrange_jobs["FIRST"].add_opt(opt, val)
  all_jobs.extend(plotinspiralrange_jobs.values())

# create the command for pylal_grbtimeslide_stats
if opts.plotgrbtimeslidestats:
  this_section = "grbtimeslidestats"
  cp = set_common_options(cp, this_section)
  pattern_dict = determine_sieve_patterns(cp, this_section, "", input_user_tag)

  if opts.second_stage:
    plotgrbtimeslidestats_job = inspiral.PlotGRBtimeslideStatsJob(cp)   
    for opt, val in pattern_dict.iteritems():
      plotgrbtimeslidestats_job.  add_opt(opt,val)
    all_jobs.append(plotgrbtimeslidestats_job)

#############################################################################
# set the usertag in the jobs
if usertag is not None:
  for job in all_jobs:
    job.add_opt('user-tag',usertag)

# Setting priority of plotting jobs.
for job in all_jobs:
  job.add_condor_cmd('priority',str(opts.priority))
 
##############################################################################
#  The meat of the DAG generation comes below
#
##############################################################################

if opts.plotinspiral:
  if opts.first_stage:
    for ifo in ifolist:
    # add an plotinspiral job
      plotinspiral_node=inspiral.PlotInspiralNode(plotinspiral_jobs[ifo + "FIRST"])
      dag.add_node(plotinspiral_node)
  if opts.second_stage:
    for ifo in ifolist:
    # add an plotinspiral job
      plotinspiral_node=inspiral.PlotInspiralNode(plotinspiral_jobs[ifo + "SECOND"])
      dag.add_node(plotinspiral_node)

if opts.plotnumtemplates:
  # add an plotnumtemplates job
  plotnumtemplates_node=inspiral.PlotNumtemplatesNode(plotnumtemplates_jobs)
  dag.add_node(plotnumtemplates_node)

if opts.plotthinca:
  # add an plotthinca job
  if opts.first_stage:
    for ifos in ifo_combo:
      combostring = combo2str(ifos)
      plotthinca_node=inspiral.PlotThincaNode(plotthinca_jobs[combostring + "FIRST"])
      dag.add_node(plotthinca_node)
  if opts.second_stage:
    for ifos in ifo_combo:
      combostring = combo2str(ifos)
      plotthinca_node=inspiral.PlotThincaNode(plotthinca_jobs[combostring + "SECOND"])
      dag.add_node(plotthinca_node)

if opts.plotethinca:
  if opts.first_stage:
    for ifos in ifo_combo:
      combostring = combo2str(ifos)
      plotethinca_node =inspiral.PlotEthincaNode(plotethinca_jobs[combostring + "FIRST"])
      dag.add_node(plotethinca_node)
  if opts.second_stage:
    for ifos in ifo_combo:
      combostring = combo2str(ifos)
      plotethinca_node =inspiral.PlotEthincaNode(plotethinca_jobs[combostring + "SECOND"])
      dag.add_node(plotethinca_node)

if opts.plotinspmissed:
  if opts.first_stage:
    for ifos in ifo_combo:
      combostring = combo2str(ifos)
      plotinspmissed_node =inspiral.PlotInspmissedNode(plotinspmissed_jobs[combostring + "FIRST"])
      dag.add_node(plotinspmissed_node)
  if opts.second_stage:
    for ifos in ifo_combo:
      combostring = combo2str(ifos)
      plotinspmissed_node =inspiral.PlotInspmissedNode(plotinspmissed_jobs[combostring + "SECOND"])
      dag.add_node(plotinspmissed_node)

if opts.ploteffdistcut:
  if opts.first_stage:
    ifos = "H1H2"
    combostring = combo2str(ifos)
    ploteffdistcut_node =inspiral.PlotEffdistcutNode(ploteffdistcut_jobs[combostring + "FIRST"])
    dag.add_node(ploteffdistcut_node)
  if opts.second_stage:
    ifos = "H1H2"
    combostring = combo2str(ifos)
    ploteffdistcut_node =inspiral.PlotEffdistcutNode(ploteffdistcut_jobs[combostring + "SECOND"])
    dag.add_node(ploteffdistcut_node)

if opts.plotinspinj:
  if opts.first_stage:
    for ifo in ifolist:
    # add an plotinspinj job
      plotinspinj_node =inspiral.PlotInspinjNode(plotinspinj_jobs[ifo + "FIRST"])
      dag.add_node(plotinspinj_node)
  if opts.second_stage:
    for ifo in ifolist:
    # add an plotinspinj job
      plotinspinj_node =inspiral.PlotInspinjNode(plotinspinj_jobs[ifo + "SECOND"])
      dag.add_node(plotinspinj_node)

if opts.plotsnrchi:
   for ifo in ifolist:
  # add an plotinspinj job
    plotsnrchi_node =inspiral.PlotSnrchiNode(plotsnrchi_jobs[ifo])
    dag.add_node(plotsnrchi_node)

if opts.plotinspiralrange:
  # add an plotnumtemplates job
  # if opts.first_stage:
  plotinspiralrange_node=inspiral.PlotInspiralrangeNode(plotinspiralrange_jobs["FIRST"])
  dag.add_node(plotinspiralrange_node)
  #if opts.second_stage:
  #  plotinspiralrange_node=inspiral.PlotInspiralrangeNode(plotinspiralrange_jobs["SECOND"])
  #  dag.add_node(plotinspiralrange_node)

if opts.plotgrbtimeslidestats:
  # add an pylal_grbtimeslide_stats job
  plotgrbtimeslides_node=inspiral.PlotGRBtimeslideStatsNode(plotgrbtimeslidestats_job)
  dag.add_node(plotgrbtimeslides_node)


# Add number of retries to plotting jobs as specified by ihope.ini
if cp.has_option("pipeline", "retry-plot-jobs"):
  num_retries = cp.getint("pipeline", "retry-plot-jobs")
  for node in dag.get_nodes():
    node.set_retry(num_retries)

##############################################################################
# Step 10: Write out the DAG, help message and log file
dag.write_sub_files()
dag.write_dag()

if opts.write_script:
  dag.write_script()

##############################################################################  
# write a message telling the user that the DAG has been written
print "\nCreated a DAG file which can be submitted by executing"
print "\n   condor_submit_dag", dag.get_dag_file()
print "\nfrom a condor submit machine (e.g. hydra.phys.uwm.edu)"

##############################################################################
# write out a log file for this script
if usertag:
  log_fh = open(basename + '.plotter.' + usertag + '.log', 'w')
else:
  log_fh = open(basename + '.plotter.log', 'w')
  
# FIXME: the following code uses obsolete CVS ID tags.
# It should be modified to use git version information.
log_fh.write( "$Id$" + "\n" )
log_fh.write( "$Name$" + "\n\n" )
log_fh.write( "Invoked with arguments:" )
for arg in command_line:
  if arg[0] == '-':
    log_fh.write( "\n" )
  log_fh.write( arg + ' ')

log_fh.write( "\n" )
log_fh.write( "Config file has CVS strings:\n" )
#log_fh.write( cp.get('pipeline','version') + "\n" )
#log_fh.write( cp.get('pipeline','cvs-tag') + "\n\n" )

log_fh.close()

sys.exit(0)

