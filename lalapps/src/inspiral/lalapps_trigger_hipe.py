# coding=utf-8
"""
trigger_hipe.in - externally triggered hierarchical inspiral pipeline driver
script
"""
from __future__ import division

__author__ = 'Patrick Brady <patrick@gravity.phys.uwm.edu>, Nickolas Fotopoulos <foton@caltech.edu>, Ian Harry <ian.harry@astro.cf.ac.uk>'

##############################################################################
import copy
import glob
import math
import sha
import os
import sys
from glue.pipeline import DeepCopyableConfigParser as dcConfigParser
import optparse
import shutil
import subprocess

import numpy

from glue.ligolw import ligolw
from glue.ligolw import utils
from glue.ligolw import table
from glue.ligolw import lsctables
from glue import lal
from glue import segments
from glue import segmentsUtils
from glue import pipeline
from pylal import grbsummary
from pylal import exttrig_dataquery
from lalapps import inspiralutils

##############################################################################
# define a few utility functions that make this job easier

def mkdir( directory, overwrite ):
  try:
    os.mkdir(directory)
  except OSError:
    if overwrite:
      print "Warning: Overwriting contents in existing directory %s" % directory
    else:
      raise

def hash_n_bits(input_string, n_bits, hash_lib=sha):
  """
  Return the first n_bits (as an integer) of a hash on input_string.  Defaults
  to SHA-1.

  n_bits=31 is useful for, say, random seeds that must fit into a LAL INT4.
  """
  # Some test examples:
  #
  # In [30]: hex(hash_n_bits("test", 28) << (32 - 28))
  # Out[30]: '0xa94a8fe0L'
  #
  # In [31]: hex(hash_n_bits("test", 29) << (32 - 29))
  # Out[31]: '0xa94a8fe0L'
  #
  # In [32]: hex(hash_n_bits("test", 30) << (32 - 30))
  # Out[32]: '0xa94a8fe4L'
  #
  # In [33]: hex(hash_n_bits("test", 31) << (32 - 31))
  # Out[33]: '0xa94a8fe4L'
  #
  # In [34]: hex(hash_n_bits("test", 32))
  # Out[34]: '0xa94a8fe5L'

  hashed_string = hash_lib.new(input_string).digest()
  n_bytes = int(math.ceil(n_bits / 8))

  def bytes2int(byte_string):
    result = 0
    for byte in byte_string:
      result <<= 8
      result += ord(byte)
    return result

  return bytes2int(hashed_string[:n_bytes]) >> (8*n_bytes - n_bits)

def createInjectionFile(hipe_dir, cp, cpinj, injrun, injection_segment,
  source_file, ipn_gps=None, usertag=None, verbose=False):
  """
  Creates an master injection file containing all injections for this run.
  Also reads the file and returns its contents
  """
  cpinj = copy.deepcopy(cpinj)

  # get the number of injections to be made
  for opt in ['exttrig-inj-start','exttrig-inj-stop']:
    value = int(cpinj.get(injrun,opt))
    cpinj.remove_option(injrun,opt)
    if 'start' in opt:
      injStart = value
    else:
      injEnd = value
  seed = hash_n_bits(hipe_dir, 31)
  numberInjections = injEnd - injStart + 1  # e.g., 1 through 5000 inclusive

      
  # get the jitter parameters
  if cpinj.has_option(injrun, "jitter-skyloc"):
    jitter_sigma_deg = cpinj.getfloat(injrun, "jitter-skyloc")
    cpinj.remove_option(injrun, "jitter-skyloc")
  else:
    jitter_sigma_deg = None


  # check if the specific Fermi systematic error needs to 
  # be added to the location jittering
  if cpinj.has_option(injrun, "jitter-skyloc-fermi"):
    jitter_skyloc_fermi = cpinj.getboolean(injrun, "jitter-skyloc-fermi")
    cpinj.remove_option(injrun, "jitter-skyloc-fermi")
  else:
    jitter_skyloc_fermi = False

  # check if we should align the total angular momentum
  if cpinj.has_option(injrun, "align-total-spin"):
    align_total_spin = cpinj.getboolean(injrun, "align-total-spin")
    cpinj.remove_option(injrun, "align-total-spin")
  else:
    align_total_spin = False

  # set all the arguments
  argument = []
  for (opt,value) in cpinj.items(injrun):
    argument.append("--%s %s" % (opt, value) )

  # add arguments on times and time-intervals
  interval = abs(injection_segment)
  injInterval = interval / numberInjections
  argument.append(" --gps-start-time %d" % injection_segment[0] )
  argument.append(" --gps-end-time %d" % injection_segment[1] )
  argument.append(" --time-interval %f" % injInterval )
  argument.append(" --time-step %f" % injInterval )
  argument.append(" --seed %d" % seed )
  argument.append(" --user-tag %s" % usertag)

  # set output file as exttrig-file or IPN file with IPN GPS time
  if ipn_gps:
    argument.append(" --ipn-gps-time %d" % ipn_gps )
  else:
    argument.append(" --exttrig-file %s" % source_file )


  # execute the command
  executable = cp.get("condor", "inspinj")
  arguments = " ".join(argument)
  inspiralutils.make_external_call(executable + " " + arguments,
    show_command=verbose)

  # recreate the output filename
  injFile = "HL-INJECTIONS_" + str(seed)
  if usertag is not None:
    injFile += "_" + usertag
  injFile += "-%d-%d.xml" % (injection_segment[0], abs(injection_segment))

  # move it into the GRB directory to avoid clutter
  new_injFile = hipe_dir + "/" + injFile
  os.rename(injFile, new_injFile)

  # jitter the sky locations of the injections
  if jitter_sigma_deg is not None:
    # rename the original, then have ligolw_cbc_jitter_skyloc create a new one
    os.rename(new_injFile, new_injFile + ".prejitter")
    cmd = ["ligolw_cbc_jitter_skyloc"]
    if jitter_skyloc_fermi:
      cmd.append("--apply-fermi-error")
    cmd.extend(["--jitter-sigma-deg",
      str(jitter_sigma_deg), "--output-file", new_injFile,
      new_injFile + ".prejitter"])
    if verbose:
      print " ".join(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
      raise subprocess.CalledProcessError(p.returncode, "%s: %s" % (" ".join(cmd), err))

  # rotate the binary so that total angular momentum has the current inclination
  if align_total_spin:
    # rename the original then have ligolw_cbc_align_total_spin create a new one
    os.rename(new_injFile, new_injFile + ".prealign")
    cmd = ["ligolw_cbc_align_total_spin", "--output-file", new_injFile,
      new_injFile + ".prealign"]
    if verbose:
      print " ".join(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
      raise subprocess.CalledProcessError(p.returncode, "%s: %s" % (" ".join(cmd), err))

  # read in the file and the tables
  doc = utils.load_filename(new_injFile)
  sims = table.get_table(doc, lsctables.SimInspiralTable.tableName)

  return sims, injInterval, numberInjections, new_injFile

##############################################################################
# inspiral_hipe wrapper

class hipe_run(object):
  """
  This class is intended to represent single run of lalapps_inspiral_hipe.
  """
  def __init__(self, hipe_dir, base_cp, ifos, log_path, source_file,
               dag_base_name, span,
               usertag=None, verbose=False, only_run=None, dont_run=None,
               do_coh_PTF=False):
    self._hipe_dir = hipe_dir
    self._log_path = log_path
    self._cp = copy.deepcopy(base_cp)
    self._ifos = ifos
    self._verbose = verbose
    self._span = span

    self.dag_base_name = dag_base_name
    self.set_span(span)

    # create directory here
    mkdir(self._hipe_dir, overwrite=True)

    ## Determine arguments to inspiral_hipe
    self._hipe_args = ["--output-segs", "--log-path=%s" % log_path,
      "--config-file=%s", "--write-script"]

    # decide what stages to run
    self.determine_stages(only_run, dont_run)

    # determine how many ifos to analyze
    n = len(self._ifos)
    if n < 1 or n > 4:
      raise ValueError, "cannot handle less than one or more than four IFOs"
    number_words = {1: "one", 2: "two", 3: "three", 4: "four"}
    for i in range(n):
      self._hipe_args.append("--%s-ifo" % number_words[i+1])

    # set individual ifos to analyze
    for ifo in self._ifos:
      self._cp.set("input", "%s-segments" % ifo.lower(), "../offSourceSeg.txt")
      self._hipe_args.append("--%s-data" % ifo.lower())

    # set the source file for thinca
    if not do_coh_PTF:
      self._cp.set('thinca','exttrig','../../'+source_file)

    # set the usertag
    self.set_usertag(usertag)

  def get_usertag(self):
    return self.usertag

  def set_usertag(self, usertag):
    self.usertag = usertag
    if usertag is not None:
      self._cp.set("pipeline", "user-tag", usertag)
    else:
      self._cp.remove_option("pipeline", "user-tag")

  def determine_stages(self, only_run=None, dont_run=None):
    """
    Determine what stages should inspiral_hipe should construct.  Store them
    in self._stages.  If only_run is not None, use those stages only.  If
    dont_run is specified, the given stages will not be run.
    """
    allowed_stages = set(["datafind", "template-bank", "inspiral",
      "coincidence", "trigbank", "inspiral-veto", "second-coinc",
      "sire-inspiral", "sire-inspiral-veto", "coire-coinc",
      "coire-second-coinc", "coherent-bank", "coherent-inspiral",
      "cohire", "summary-coherent-inspiral-triggers" ])

    if only_run is not None:
      if not (set(only_run) <= allowed_stages):
        raise ValueError, "stage not in allowed stages"
      self._stages = only_run
    else:
      if dont_run is not None:
        allowed_stages -= set(dont_run)
      self._stages = allowed_stages

  def set_span(self, span):
    """
    Set start and end times, required for cache to be produced.
    """
    self._cp.set("input", "gps-start-time", str(span[0]))
    self._cp.set("input", "gps-end-time", str(span[1]))
    self._span = span

  def set_numslides(self, numslides):
    if numslides == 0: numslides = ""
    self._cp.set('input', 'num-slides', str(numslides))

  def remove_longslides(self):
    if self._cp.has_option('input','do-long-slides'):
      self._cp.remove_option('input','do-long-slides')

  def set_longslides(self):
    if not self._cp.has_option('input','do-long-slides'):
      self._cp.set('input','do-long-slides','')

  def remove_shortslides(self):
    if self._cp.has_option('coh_PTF_inspiral','do-short-slides'):
      self._cp.remove_option('coh_PTF_inspiral','do-short-slides')

  def set_injections(self, injrun, numberInjFiles):
    """
    Turn this analysis into an injection run, using the injrun section from
    the config file cpinj.
    """
    # set the start and stop seeds from the injection config file
    self._cp.set('pipeline', 'exttrig-inj-start', '1')
    self._cp.set('pipeline', 'exttrig-inj-stop', str(numberInjFiles))
    self._cp.set("condor", "inspinj", "/bin/true")
    self._hipe_args.append("--noop-inspinj")

  def write_dag(self):
    """
    Create directory and run lalapps_inspiral_hipe to create the DAG.
    """
    # write out an appropriate ini file
    ini_file_name = self.dag_base_name + ".ini"
    os.chdir(self._hipe_dir)
    if not os.path.isdir(self._log_path):  # in case of relative path
      os.mkdir(self._log_path)
    cp_file = open(ini_file_name, 'w')
    self._cp.write(cp_file)
    cp_file.close()

    arguments = " ".join(self._hipe_args)
    arguments += " " + " ".join(["--" + stage for stage in self._stages])

    # run inspiral_hipe
    inspiralutils.make_external_call(\
      " ".join((self._cp.get("condor", "inspiral_hipe"),
               arguments % ini_file_name)),
      show_stdout=False, show_command=self._verbose)

    # create link to datafind directory's cache directory
    if not "datafind" in self._stages:
      if os.path.isdir("cache") and not os.path.islink("cache"):
        os.rmdir("cache")
      if not os.path.exists("cache"):
        os.symlink("../datafind/cache", "cache")

    # set inspinj jobs to run in the local universe
    inspinj_sub_files = glob.glob("*inspinj*.sub")
    for sub_file in inspinj_sub_files:
      contents = open(sub_file).read().replace("universe = standard", "universe = local")
      open(sub_file, "w").write(contents)

    os.chdir("../../")

  def get_cache_name(self):
    return self._hipe_dir + "/" + \
      inspiralutils.hipe_cache(self._ifos, self.usertag, self._span[0],
                               self._span[1])

  def get_dag_name(self):
    dag_name = self.dag_base_name
    if self.usertag is not None:
      dag_name = dag_name + "." + self.usertag
    dag_name = dag_name + ".dag"
    return dag_name

  def get_dag_path(self):
    """
    Determine the full path of the DAG
    """
    return self._hipe_dir + "/" + self.get_dag_name()

  def create_dag_node(self):
    """
    Return a CondorDAGNode that represents this entire DAG.
    """
    dir, fname = os.path.split(self.get_dag_path())
    job = pipeline.CondorDAGManJob(fname, dir)
    node = pipeline.CondorDAGNode(job)
    return node

  def finalize(self, uberdag=None, parent=None):
    """
    Write the DAG and its sub_file.  Add node to uberdag and add dependency
    to parent, if these kwargs are not None.  Return the CondorDAGNode
    representing this DAG.
    """
    self.write_dag()
    node = self.create_dag_node()
    if parent is not None:
      node.add_parent(parent)
    if uberdag is not None:
      uberdag.add_node(node)
    return node

##############################################################################
# coh_PTF hipe wrapper

class ptf_run(hipe_run):
  """
  This class is intended to represent a single run of lalapps_cohPTF_hipe
  """
  def __init__(self, hipe_dir, base_cp, ifos, log_path, source_file,
               dag_base_name, span,usertag=None, verbose=False,
               only_run=None, dont_run=None):
    hipe_run.__init__(self, hipe_dir, base_cp, ifos, log_path, source_file,
               dag_base_name, span,
               usertag=usertag, verbose=verbose, only_run=only_run,
               dont_run=dont_run,do_coh_PTF=True)

    ## Determine arguments to cohPTF_hipe
    temp_args=[]
    for arg in self._hipe_args:
      if arg != "--noop-inspinj" and arg != "--ringdown":
        temp_args.append(arg)
    self._hipe_args = temp_args

    # decide what stages to run
    self.determine_PTF_stages(only_run, dont_run)

  def determine_PTF_stages(self, only_run=None, dont_run=None):
    """
    Determine what stages should inspiral_hipe should construct.  Store them
    in self._stages.  If only_run is not None, use those stages only.  If
    dont_run is specified, the given stages will not be run.
    """
    allowed_stages = set(["inspiral","splitbank" ])

    if only_run is not None:
      if not (set(only_run) <= allowed_stages):
        raise ValueError, "stage not in allowed stages"
      self._stages = only_run
    else:
      if dont_run is not None:
        allowed_stages -= set(dont_run)
      self._stages = allowed_stages
      
  def write_dag(self):
    """
    Create directory and run lalapps_cohPTF_hipe to create the DAG.
    """
    # write out an appropriate ini file
    ini_file_name = self.dag_base_name + ".ini"
    os.chdir(self._hipe_dir)
    if not os.path.isdir(self._log_path):  # in case of relative path
      os.mkdir(self._log_path)
    cp_file = open(ini_file_name, 'w')
    self._cp.write(cp_file)
    cp_file.close()

    arguments = " ".join(self._hipe_args)
    arguments += " " + " ".join(["--" + stage for stage in self._stages])

    # run inspiral_hipe
    inspiralutils.make_external_call(\
      " ".join((self._cp.get("condor","cohPTF_hipe"),
               arguments % ini_file_name)),
      show_stdout=False, show_command=self._verbose)

    if os.path.isdir("cache") and not os.path.islink("cache"):
      os.rmdir("cache")
    if not os.path.exists("cache"):
      os.symlink("../datafind/cache", "cache")

    if cp.has_option('splitbank-meta','bank-file'):
      shutil.copy("../../" + cp.get('splitbank-meta','bank-file'),'.')
    if cp.has_option('coh_PTF_inspiral','bank-veto-templates'):
      shutil.copy("../../" + cp.get('coh_PTF_inspiral','bank-veto-templates'),'.')

    os.chdir("../../")

##############################################################################
# convenience functions

def concatenate_files(input_names, output_name):
  """
  Concatenate several files to a file named output_cache.
  """
  inspiralutils.make_external_call("cat %s > %s" % \
                                   (" ".join(input_names), output_name))

def make_ligolw_add(orig_cache, pattern, outfile, logdir, cp, ligolw_job=None):
  """
  Take a cache, sieve it for pattern, create a cache-file for ligolw_add,
  and create a LigolwAddNode that will create a new
  output cache with the matching files replaced by the single ligolw_added
  file.
  """
  # determine the files to ligolw_add
  sub_cache = orig_cache.sieve(description=pattern)
  if len(sub_cache) == 0:
    print >>sys.stderr, "warning: no files on which to run ligolw_add"
    return None

  # create the cache-file
  cachefile = os.path.basename( outfile )[:-3]+'cache'
  sub_cache.tofile(open(cachefile,'w'))

  if ligolw_job is None:
    ligolw_job = pipeline.LigolwAddJob(logdir, cp)
    ligolw_job.set_universe("local")

  node = pipeline.LigolwAddNode(ligolw_job)
  node.add_output_file(outfile)
  node.add_var_opt("input", cachefile)
  node.add_var_opt("output", outfile)

  # return the cache-file without the ones
  # just extracted above
  new_cache = lal.Cache([entry for entry in orig_cache if entry not in sub_cache])
  new_cache.extend(lal.Cache.from_urls([outfile]))
  return node, new_cache

##############################################################################
# define usage and command line options and arguments - parse
usage = """usage: %prog ...

Set up a some Condor DAGs to generate CBC candidates associated with an
external trigger, plus the off-source and time-slide triggers required for
background estimation and software injections for efficiency measurements.

Directory structure:

GRBXXXXXX/
  datafind/
  injectionXX/
  injectionYY/
  ...
  onoff/
"""
parser = optparse.OptionParser(usage, version="$Id$")

parser.add_option("-v", "--verbose", action="store_true",default=False,\
  help="make things verbose" )

# read GRB information
parser.add_option("", "--name", action="append",
  default=[], help="Name of the GRB (e.g. 100122A).")
parser.add_option("", "--time", action="store", type="int",
  default=None, help="GPS time of the GRB trigger.")
parser.add_option("", "--ra", action="store", type="float",
  default=None, help="Right ascension of the GRB in degrees.")
parser.add_option("", "--dec", action="store", type="float",
  default=None, help="Declination of the GRB in degrees.")
parser.add_option("", "--grb-file", action="store", type="string",
  default=None, metavar=" FILE", help="Full path to GRB xml file.")

# what parameters?
parser.add_option("--offset",action="store",type="int",\
    default=2500, help="Length of time to query for segments in both directions (seconds).")
parser.add_option("-c","--number-buffer-left",action="store",type="int",\
    default=0, metavar=" BUFFERLEFT",\
    help="Specifies the number of buffer trials to the left of the on-source "\
    "segment (default: 0)")
parser.add_option("-d","--number-buffer-right",action="store",type="int",\
    default=0, metavar=" BUFFERRIGHT",\
    help="Specifies the number of buffer trials to the right of the on-source "\
    "segment (default: 0)")
parser.add_option("-t", "--padding-time", type="int", default=0,
  help="padding time on each side of an analysis segment in seconds")
parser.add_option("-n", "--num-trials", type="int",
  help="the off-source segment's full extent should contain NUM_TRIALS "\
  "segments of the length of the on-source segment; this includes the sum "\
  "of on-source, buffer, and off-source trials.")
parser.add_option("-s", "--symmetric", action="store_true", default=False,
  help="restrict the off-source trials to be arranged symmetrically about "\
       "the on-source segment (default: asymmetrical distribution allowed)")
parser.add_option("-o", "--overwrite-dir", action="store_true",default=False,\
  help="overwrite contents in already existing directories" )

# read in the config file
parser.add_option("-f","--config-file",action="store",type="string",\
  default=None, metavar=" FILE", help="use configuration file FILE" )
parser.add_option("-g","--injection-config",action="store",type="string",\
  default=None, metavar=" FILE", help="use configuration file FILE" )
parser.add_option("-p", "--log-path",action="store",type="string",\
    metavar=" PATH",help="directory to write condor log file")
parser.add_option("-u","--user-tag",action="store",type="string",\
    default=None, help="user defined tag for naming the DAG file" )

# Add some extra capabilities
parser.add_option("", "--skip-datafind", action="store_false", default=True,
  dest="do_datafind", help="skip a query for data file locations (place "\
  "the file caches in ./GRBxxxxxx/datafind/cache).")
parser.add_option("", "--skip-onoff", action="store_false", default=True,
  dest="do_onoff", help="skip the onoff DAG")
parser.add_option("", "--do-slides", action="store_true", default=False,
  dest="do_slides", help="do time slides")
parser.add_option("", "--do-coh-PTF", action="store_true", default=False,
  dest="do_coh_PTF", help="Use this flag if doing a coh PTF analysis")
parser.add_option("", "--ipn", action="store_true", default=False,\
  dest="ipn", help="Use this flag if analysing a GRB detected by the IPN, default: %default")
parser.add_option("", "--extend", action="store_true", default=False,
  help="")
parser.add_option("", "--useold", action="store_true", default=False,
  help="")
parser.add_option("", "--make-plots", action="store_true", default=False,
  help="")
parser.add_option("", "--make-xml", action="store_true", default=False,
  help="")

( opts , args ) = parser.parse_args()

required_opts = ["name", "time", "padding_time", "config_file", "injection_config", "log_path"]
for opt in required_opts:
    if getattr(opts, opt) is None:
        raise ValueError, "--%s is a required option" % opt
if not opts.grb_file and (not opts.time or not opts.name):
    raise ValueError, "Either a valid GRB xml file must be specified or the GPS time and name of the GRB!"

##############################################################################
# find available data

if opts.grb_file:
  xmldoc    = utils.load_filename(opts.grb_file, gz=opts.grb_file.endswith('.gz'))
  ext_table = table.get_table(xmldoc,lsctables.ExtTriggersTable.tableName)
  grb_time  = ext_table[0].start_time
  grb_name  = os.path.basename(opts.grb_file)[3:-4]
  grb_ra    = ext_table[0].event_ra
  grb_dec   = ext_table[0].event_dec
else:
  grb_name = opts.name[0]
  grb_time = int(opts.time)
  grb_ra   = float(opts.ra)
  grb_dec  = float(opts.dec)

exttrig_config_file, grb_ifolist, onSourceSegment, offSourceSegment = exttrig_dataquery.exttrig_dataquery(grb_name, grb_time, grb_ra, grb_dec, opts.offset, opts.config_file, opts.extend, opts.useold, opts.make_plots, opts.make_xml)

##############################################################################
# create the config parser object and exttrig_dataquery ini file
cp = dcConfigParser()
cp.read(exttrig_config_file)

ext_trigs = grbsummary.load_external_triggers('grb%s.xml'%opts.name[0])

if ext_trigs is None:
  print >>sys.stderr, "No external triggers found.  Nothing to do."
  sys.exit()

if len(opts.name) > 0:
  temp_list = filter(lambda row: row.event_number_grb in opts.name, ext_trigs)
  ext_trigs = table.new_from_template(ext_trigs)
  ext_trigs.extend(temp_list)

  if len(ext_trigs) != len(opts.name):
    missing = set(opts.name) - set([row.event_number_grb for row in ext_trigs])
    raise ValueError, "could not find the following requested GRBs: " \
      + "".join(missing)

if opts.verbose:
  print "Will construct workflows to analyze:", \
    ", ".join(["GRB" + trig.event_number_grb for trig in ext_trigs])

if opts.do_coh_PTF:
  for program in ['coh_PTF_inspiral','coh_PTF_spin_checker']:
    if not opts.ipn:
      cp.set(program,'right-ascension',str(ext_trigs[0].event_ra))
      cp.set(program,'declination',str(ext_trigs[0].event_dec))
    cp.set(program,'trigger-time',str(ext_trigs[0].start_time))
    cp.set(program,'trigger-time-ns',str(ext_trigs[0].start_time_ns))

############################################################################
# set up the über dag for all intervals and all injections

tag = opts.config_file.rstrip(".ini")
if opts.user_tag:
  tag += '_'+opts.user_tag
uberdag = pipeline.CondorDAG("%s/%s_uberdag.log" % (opts.log_path, tag))
uberdag.set_dag_file("%s_uberdag" % tag)

##############################################################################
# loop over the GRBs and construct their required sub-DAGs

grb_caches = []

for grb in ext_trigs:

  # name and the directory
  idirectory = "GRB" + str(grb.event_number_grb)
  if opts.verbose:
    print "* Constructing workflow for", idirectory
  mkdir(idirectory, opts.overwrite_dir)

  ##############################################################################
  # create the source-file which contains data about the position and time
  # and place it into the corresponding directory, also set the IPN file for the
  # very same purpose when analysing IPN GRBs
  source_file = idirectory+"/triggerGRB"+str(grb.event_number_grb)+".xml"
  grbsummary.write_rows([grb], lsctables.ExtTriggersTable, source_file)

  ##############################################################################
  # set up the segment including the off-source segment

  grb_ifolist.sort()
  ifo_times = "".join(grb_ifolist)

  if offSourceSegment is None:
    print >>sys.stderr, "Warning: insufficient multi-IFO data to construct an off-source segment for GRB %s; skipping" % grb.event_number_grb
    continue
  elif opts.verbose:
    print "Sufficient off-source data has been found in", ifo_times, "time."

  # write out the segment list to a segwizard file
  offsource_segfile = idirectory + "/offSourceSeg.txt"
  segmentsUtils.tosegwizard(open(offsource_segfile, "w"),
                            segments.segmentlist([offSourceSegment]))
  onsource_segfile = idirectory + "/onSourceSeg.txt"
  segmentsUtils.tosegwizard(file(onsource_segfile, "w"),
                            segments.segmentlist([onSourceSegment]))
  segLen = abs( onSourceSegment )
  bufferSegment = segments.segment( onSourceSegment[0]-opts.number_buffer_left*segLen,\
                                    onSourceSegment[1]+opts.number_buffer_right*segLen)
  buffer_segfile = idirectory + "/bufferSeg.txt"
  segmentsUtils.tosegwizard(file(buffer_segfile, "w"),
                            segments.segmentlist([bufferSegment]))

  if opts.verbose:
    print "on-source segment: ", onSourceSegment
    print "off-source segment: ", offSourceSegment

  ############################################################################
  # set up the analysis dag for this interval
  #
  # In doing this, we simply drop the configuration file into the
  # sub-directory, modify as needed, and then run the appropriate DAG
  # generation script.  In slightly more detail, steps are:
  #
  # 1. datafind
  # 2. onoff-source analysis
  # 3. time-slide analysis
  # 4. injection runs
  #    a. generate injections
  #    b. split injections into files so they are safely separated
  #    c. generate DAG

  hipe_caches = []

  if opts.do_datafind:
    if opts.verbose: print "Creating datafind DAG"
    datafind_run = ["datafind"]
    if opts.do_coh_PTF:
      datafind_run.append("template-bank")
    datafind_dag = hipe_run(idirectory + "/datafind", cp, grb_ifolist,
                            opts.log_path, source_file, "datafind",
                            offSourceSegment,
                            usertag=idirectory + "_DATAFIND",
                            verbose=opts.verbose, only_run=datafind_run,
                            do_coh_PTF=opts.do_coh_PTF)
    datafind_node = datafind_dag.finalize(uberdag=uberdag, parent=None)
  else:
    # the other nodes will try to use datafind_node as a parent; set to None.
    datafind_node = None

  if opts.do_onoff:
    if opts.verbose: print "Creating DAG for zero lag analysis"
    if opts.do_coh_PTF:
      onoff_analysis = ptf_run(idirectory + "/onoff", cp, grb_ifolist,
                        opts.log_path, source_file, "zero_lag",
                        offSourceSegment,
                        usertag=idirectory + "_ZERO_LAG_CATEGORY_1",
                        verbose=opts.verbose)
      onoff_analysis.remove_longslides()
    else:
      onoff_analysis = hipe_run(idirectory + "/onoff", cp, grb_ifolist,
                        opts.log_path, source_file, "zero_lag",
                        offSourceSegment,
                        usertag=idirectory + "_ZERO_LAG_CATEGORY_1",
                        verbose=opts.verbose,
                        dont_run=["datafind"])
      onoff_analysis.set_numslides(0)
    onoff_node = onoff_analysis.finalize(uberdag=uberdag, parent=datafind_node)
    hipe_caches.append(onoff_analysis.get_cache_name())

  if opts.do_slides:
    if opts.verbose: print "Creating DAG for time slide analysis"
    if opts.do_coh_PTF:
      slide_analysis = ptf_run(idirectory + "/timeslides", cp, grb_ifolist,
                        opts.log_path, source_file, "slides",
                        offSourceSegment,
                        usertag=idirectory + "_TIME_SLIDES_CATEGORY_1",
                        verbose=opts.verbose)
      slide_analysis.set_longslides()
    else:
      slide_analysis = hipe_run(idirectory + "/timeslides", cp, grb_ifolist,
                        opts.log_path, source_file, "slides",
                        offSourceSegment,
                        usertag=idirectory + "_TIME_SLIDES_CATEGORY_1",
                        verbose=opts.verbose,
                        dont_run=["datafind"])
    slide_node = slide_analysis.finalize(uberdag=uberdag, parent=datafind_node)
    hipe_caches.append(slide_analysis.get_cache_name()) 
  ############################################################################
  # create the config parser object and read in the ini file
  if opts.injection_config:
    cpinj = dcConfigParser()
    cpinj.read(opts.injection_config)

    ############################################################################
    # set up the injection dag for this interval
    for injrun in cpinj.sections():
      if opts.verbose: print "Creating DAG for", injrun
      hipe_dir = os.path.join(idirectory, injrun)
      usertag = idirectory + "_" + injrun

      if opts.do_coh_PTF:
        injection_analysis = ptf_run(hipe_dir, cp, grb_ifolist,
          opts.log_path, source_file, injrun, offSourceSegment,
          usertag=usertag, verbose=opts.verbose)
      else:
        injection_analysis = hipe_run(hipe_dir, cp, grb_ifolist,
          opts.log_path, source_file, injrun, offSourceSegment,
          usertag=usertag, verbose=opts.verbose,
          dont_run=["datafind", "sire-inspiral", "sire-inspiral-veto",
                    "coire-coincidence"])

      # create the master injection file and get its content
      if opts.verbose: print "  Creating injections..."
      if opts.ipn:
        ipnstart = ext_trigs[0].start_time
      else:
        ipnstart = None
      sims, injInterval, numberInjections, injFile =\
          createInjectionFile(hipe_dir, cp, cpinj, injrun,\
                              offSourceSegment.contract(opts.padding_time),\
                              source_file, ipnstart, usertag, opts.verbose)

      new_injfiles = [injFile]

      # split the injections
      if opts.verbose: print "  Splitting injections..."
      safetyInterval = 800 # safety time in seconds between two injections
      deltaIndex = int(safetyInterval // injInterval) + 1
      for i in range(deltaIndex):
          inj_file_name = "%s/HL-INJECTION_%s_%d-%d-%d.xml" % \
              (hipe_dir, usertag, i+1, offSourceSegment[0],
               abs(offSourceSegment))
          sims_chunk = sims[i:numberInjections:deltaIndex]
          grbsummary.write_rows(sims_chunk, lsctables.SimInspiralTable, inj_file_name)
          new_injfiles.append(inj_file_name)

      # fixing the right parameters
      if opts.verbose: print "  Writing DAG..."
      if not opts.do_coh_PTF:
        injection_analysis.set_numslides(0)
      else:
        injection_analysis.remove_longslides()
        injection_analysis.remove_shortslides()
      injection_analysis.set_injections(injrun, deltaIndex)
      injection_node = injection_analysis.finalize(uberdag=uberdag,
                                                   parent=datafind_node)

      if not opts.do_coh_PTF:
        # ligolw_add FOUND and MISSED files
        injection_cache_filename = injection_analysis.get_cache_name()
        injection_cache = lal.Cache.fromfile(open(injection_cache_filename))
        ligolw_added_cache_filename = injection_cache_filename.replace("HIPE",
          "HIPE_LIGOLW_ADDED")

        # add a ligolw_add command to concatenate all the found injections
        found_filename = "%s/%s-COIRE_INJECTION_FOUND_SECOND_%s_%s-%d-%d.xml" %\
          (hipe_dir, "".join(grb_ifolist), "".join(grb_ifolist), usertag,
           offSourceSegment[0], abs(offSourceSegment))
        ligolw_add_1, new_cache = make_ligolw_add(injection_cache,
                                       "COIRE_*_FOUND_SECOND", found_filename,
                                       opts.log_path, cp)
        ligolw_add_1.add_parent(injection_node)
        uberdag.add_node(ligolw_add_1)

        # add a ligolw_add command to concatenate all the missed injections
        missed_filename = found_filename.replace("FOUND", "MISSED")
        ligolw_add_2, new_cache = make_ligolw_add(new_cache,
                                  "COIRE_*_MISSED_SECOND", missed_filename,
                                  opts.log_path, cp)
        ligolw_add_2.add_parent(injection_node)

        # Create a temporay cache-file containing two files (F/M)
        # which will be created by the ligolw_add step
        new_cache.extend(lal.Cache.from_urls(new_injfiles, coltype=int))
        new_cache.tofile(open(ligolw_added_cache_filename, "w"))
        uberdag.add_node(ligolw_add_2)

        # add these two files to the regular main cache
        hipe_caches.append(ligolw_added_cache_filename)

  if opts.user_tag:
    grb_cache_name = idirectory + "/" + idirectory +"_"+ opts.user_tag+".cache"
  else:
    grb_cache_name = idirectory + "/" + idirectory + ".cache"
  concatenate_files(hipe_caches, grb_cache_name)
  grb_caches.append(grb_cache_name)

if opts.verbose: print "Writing sub files and uber-DAG..."
uberdag.write_sub_files()
uberdag.write_dag()

concatenate_files(grb_caches, tag + ".cache")

print """
If you run the uber-dag, you probably have to tell Condor not to
complain that your DAG logs are on NFS volumes:

bash users:
export _CONDOR_DAGMAN_LOG_ON_NFS_IS_ERROR=FALSE

tcsh users:
setenv _CONDOR_DAGMAN_LOG_ON_NFS_IS_ERROR FALSE
"""
