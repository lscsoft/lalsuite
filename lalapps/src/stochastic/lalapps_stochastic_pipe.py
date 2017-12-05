"""
stochastic_pipe.in - SGWB Standalone Analysis Pipeline
                   - Pipeline DAG Driver Script

Copyright (C) 2004-2006 Adam Mercer

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
"""

## \file
#
# <dl>
# <dt>Name</dt><dd>
# <tt>lalapps_stochastic_pipe</tt> --- python script to generate Condor DAGs
# to run the stochastic pipeline.</dd>
#
# <dt>Synopsis</dt><dd>
# <tt>lalapps_stochastic_pipe</tt>
# [<tt>--help</tt>]
# [<tt>--version</tt>]
# [<tt>--user-tag</tt> <i>TAG</i>]
# [<tt>--datafind</tt>]
# [<tt>--stochastic</tt>]
# [<tt>--stopp</tt>]
# [<tt>--priority</tt>]
# <tt>--config-file</tt> <i>FILE</i>
# <tt>--log-path</tt> <i>PATH</i></dd>
#
# <dt>Description</dt><dd>
# <tt>lalapps_stochastic_pipe</tt> builds a stochastic search DAG suitable
# for running at the various LSC Data Grid sites. The configuration file
# specifies the parameters needed to run the analysis jobs contained
# within the pipeline. It is specified with the <tt>--config-file</tt>
# option. Examples of typical .ini files can be found in the directory
# <tt>lalapps/src/stochastic/example</tt>.
#
# The .ini file contains several sections. The <tt>[condor]</tt> section
# contains the names of executables which will run the various stages of
# the pipeline. The <tt>[pipeline]</tt> section gives the CVS details of the
# pipeline. The <tt>[input]</tt> section specifies the segment list and the
# minimum and maximum segment duration for the jobs. The <tt>[datafind]</tt>
# section specifies the frame types for the two input streams for passing
# onto LSCdataFind. The <tt>[detectors]</tt> section specifies the detector
# pair to cross correlate. The <tt>[calibration]</tt> section specifies the
# location and name of the various calibration cache files. The
# <tt>[stochastic]</tt> section specifies static options to pass to
# <tt>lalapps_stochastic</tt>, i.e. options that are not automatically
# generated. Finally the <tt>[stopp]</tt> section specifies the static
# options to pass onto <tt>lalapps_stopp</tt>.</dd>
#
# <dt>Options</dt><dd>
# <dl>
# <dt><tt>--help</tt></dt><dd> Display usage information</dd>
# <dt><tt>--version</tt></dt><dd> Display version information</dd>
# <dt><tt>--user-tag</tt> <i>TAG</i></dt><dd> The tag for the job</dd>
# <dt><tt>--datafind</tt></dt><dd> Run LSCdataFind as part of the DAG to create
# the cache files for each science segment</dd>
# <dt><tt>--stochastic</tt></dt><dd> Run <tt>lalapps_stochastic</tt> on the data</dd>
# <dt><tt>--stopp</tt></dt><dd> Run <tt>lalapps_stopp</tt> on the data</dd>
# <dt><tt>--priority</tt> <i>PRIO</i></dt><dd> Run jobs with condor priority <i>PRIO</i></dd>
# <dt><tt>--config-file</tt> <i>FILE</i></dt><dd> Use configuration file <i>FILE</i></dd>
# <dt><tt>--log-path</tt> <i>PATH</i></dt><dd> Directory to write condor log file</dd>
# </dl></dd>
#
# <dt>Configuration File Options</dt><dd>
# The configuration file is a standard python configuration that can be
# easily parsed using the standard ConfigParser module from python. The
# configuration file consists of sections, led by a "[section]" header
# and followed by "name: value" entries. The optional values can contain
# format strings which refer to other values in the same section. Lines
# beginning with "\#" or ";" are ignored and may be used to provide
# comments.
#
# The first section required is "[condor]", this specfies all the
# parameters associated with condor.
# <dl>
# <dt>universe</dt><dd>
# Specfies the condor universe under which to run, this should be set to
# "standard" if LALApps has been compiled with the
# <tt>--enable-condor</tt>, otherwise it should be set to "vanilla".</dd>
# <dt>datafind</dt><dd>
# Specifies the location of the <tt>LSCdataFind</tt> script.</dd>
# <dt>stochastic</dt><dd>
# Specifies the location of the <tt>lalapps_stochastic</tt> binary.</dd>
# <dt>stopp</dt><dd>
# Specifies the location of the <tt>lalapps_stopp</tt> binary.</dd>
# </dl>
#
# The next section is "[pipeline]", this specifies version information
# of the pipeline, currently only the version of configuration file is
# specified.
# <dl>
# <dt>version</dt><dd>
# A CVS <tt>\\f$Id\\f$</tt> keyword specifying the version of the
# configuration file.</dd>
# </dl>
#
# The next section is "[input]" which specifies parameters regarding
# the input data.
# <dl>
# <dt>segments</dt><dd>
# Specifies the location of the segment list, outputed from
# <tt>SegWizard</tt> listing the Science Segments to analyse.</dd>
# <dt>min_length</dt><dd>
# Specifies the minimum length of Science Segment to analyse</dd>
# <dt>max_length</dt><dd>
# Specifies the maximum length at which the split the Science Segments
# into.</dd>
# <dt>channel</dt><dd>
# Currently the pipeline infrastructure requires that the "[input]"
# section contains the option "channel", this can be set to anything as
# it not used within the \c LSCdataFind class.</dd>
# </dl>
#
# The next section is "[datafind]" which specifies parameters for
# <tt>LSCdataFind</tt>.
# <dl>
# <dt>type-one</dt><dd>
# Specifies the frame type of the first data stream for <tt>LSCdataFind</tt>
# to find.</dd>
# <dt>type-two</dt><dd>
# Specifies the frame type of the second data stream for
# <tt>LSCdataFind</tt> to find.</dd>
# <dt>server</dt><dd>
# Specifies which LDRdataFindServer to use.</dd>
# </dl>
#
# The next section is "[detectors]" which specifies the detectors to use
# for the search.
# <dl>
# <dt>detector-one</dt><dd>
# Specifies the detector for the first data stream.</dd>
# <dt>detector-two</dt><dd>
# Specifies the detector for the second data stream.</dd>
# </dl>
#
# The next section is "[calibration]" which specifies the location of
# the calibration files.
# <dl>
# <dt>path</dt><dd>
# Specifies the path to the calibration cache files.</dd>
# <dt>L1</dt><dd>
# Specifies the cache file for the Livingston 4 km detector.</dd>
# <dt>H1</dt><dd>
# Specifies the cache file for the Hanford 4 km detector.</dd>
# <dt>H2</dt><dd>
# Specifies the cache file for the Hanford 2 km detector.</dd>
# </dl>
#
# The next section is "[stochastic]", this is used to specify any static
# options to pass onto <tt>lalapps_stochastic</tt>, i.e. options that are
# not automatically generated. See
# \ref stochastic.c for the accepted options for
# <tt>lalapps_stochastic</tt>. The options that are automatically generated
# are, <tt>--gps-start-time</tt>, <tt>--gps-end-time</tt>,
# <tt>--ifo-one</tt>, <tt>--ifo-two</tt>, <tt>--frame-cache-one</tt>,
# <tt>--frame-cache-two</tt>, <tt>--calibration-cache-one</tt>, and\\
# <tt>--calibration-cache-two</tt>.
#
# The next section is "[stopp]", this is used to specfiy any static
# options to pass onto <tt>lalapps_stopp</tt>, i.e. options that are not
# automatically generated. See \ref stopp.c for the
# accepted options for <tt>lalapps_stopp</tt>. The XML output files from
# <tt>lalapps_stochastic</tt> are automatically added.</dd>
#
# <dt>Example</dt><dd>
# Generate a DAG to run a stochastic search on a pair of interferometers
# specified in the configuration file. The generated DAG is then submitted
# with \c condor_submit_dag
#
# \code
# > lalapps_stochastic_pipe --log-path /home/ram/dag_logs \
#      --datafind --stochastic --stopp --config-file stochastic_H1L1.ini
# > condor_submit_dag stochastic_H1L1.dag
# \endcode</dd>
#
# <dt>Author</dt><dd>
# Adam Mercer</dd>
# </dl>

__author__ = 'Adam Mercer <ram@star.sr.bham.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'

# import required modules
import sys
import os
import getopt
import re
import string
import tempfile
import ConfigParser

# import the lalapps pipeline modules
from glue import pipeline
import stochastic

# program usage
def usage():
  msg = """\
Usage: lalapps_stochastic_pipe [options]

  -h, --help               display this message
  -v, --version            print version information and exit
  -u, --user-tag TAG       tag the jobs with TAG

  -d, --datafind           run LSCdataFind to create frame cache files
  -s, --stochastic         run lalapps_stochastic
  -t, --stopp              run lalapps_stopp

  -P, --priority PRIO      run jobs with condor priority PRIO

  -f, --config-file FILE   use configuration file FILE
  -l, --log-path PATH      directory to write condor log file
  """
  print >> sys.stderr, msg

# parse the command line options to figure out what we should do
shortop = "hvu:dstP:f:l:"
longop = [
  "help",
  "version",
  "user-tag=",
  "datafind",
  "stochastic",
  "stopp",
  "priority=",
  "config-file=",
  "log-path="
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

# default options
usertag = None
do_datafind = None
do_stochastic = None
do_stopp = None
condor_prio = None
config_file = None
log_path = None

# process options
for o, a in opts:
  if o in ("-v", "--version"):
    print "lalapps_stochastic_pipe version", __version__
    print "Built on top of stochastic.py version", stochastic.version()
    sys.exit(0)
  elif o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-u", "--user-tag"):
    usertag = a
  elif o in ("-d", "--datafind"):
    do_datafind = 1
  elif o in ("-s", "--stochastic"):
    do_stochastic = 1
  elif o in ("-t", "--stopp"):
    do_stopp = 1
  elif o in ("-P", "--priority"):
    condor_prio = a
  elif o in ("-f", "--config-file"):
    config_file = a
  elif o in ("-l", "--log-path"):
    log_path = a
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

if not config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use --config-file FILE to specify location."
  sys.exit(1)

if not log_path:
  print >> sys.stderr, "No log file path specified."
  print >> sys.stderr, "Use --log-path PATH to specify a location."
  sys.exit(1)

# try and make a directory to store the cache files and job logs
try: os.mkdir('cache')
except: pass
try: os.mkdir('logs')
except: pass

# create the config parser object and read in the ini file
cp = ConfigParser.ConfigParser()
cp.read(config_file)

# create a log file that the Condor jobs will write to
basename = re.sub(r'\.ini',r'',config_file)
tempfile.tempdir = log_path
if usertag:
  tempfile.template = basename + '.' + usertag + '.dag.log'
else:
  tempfile.template = basename + '.dag.log.'
logfile = tempfile.mktemp()
fh = open( logfile, "w" )
fh.close()

# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile)
if usertag:
  dag.set_dag_file(basename + '.' + usertag)
else:
  dag.set_dag_file(basename)

# create the Condor jobs that will be used in the DAG
df_job = stochastic.LSCDataFindJob('cache','logs',cp)
stoch_job = stochastic.StochasticJob(cp)
stopp_job = stochastic.StoppJob(cp)

# set file names
if usertag:
  subsuffix = '.' + usertag + '.sub'
else:
  subsuffix = '.sub'
df_job.set_sub_file(basename + '.datafind' + subsuffix)
stoch_job.set_sub_file(basename + '.stochastic' + subsuffix)
stopp_job.set_sub_file(basename + '.stopp' + subsuffix)

# set the usertag in the jobs
if usertag:
  stoch_job.add_opt('user-tag',usertag)

# set the condor job priority
if condor_prio:
  df_job.add_condor_cmd('priority',condor_prio)
  stoch_job.add_condor_cmd('priority',condor_prio)
  stopp_job.add_condor_cmd('priority',condor_prio)

# get the pad and chunk lengths from the values in the ini file
min_length = int(cp.get('input','min_length'))
max_length = int(cp.get('input','max_length'))

# read science segs that are greater or equal to a chunk from the input file
data = pipeline.ScienceData()
data.read(cp.get('input','segments'),min_length)

# create the chunks from the science segments
data.make_chunks(max_length)
data.make_short_chunks_from_unused(min_length)

# get the ifos
ifo1 = cp.get('detectors','detector-one')
ifo2 = cp.get('detectors','detector-two')

# get the frame types
type1 = cp.get('datafind','type-one')
type2 = cp.get('datafind','type-two')

# get the LSCdataFind server
server = cp.get('datafind','server')

# set up stopp job
stopp = stochastic.StoppNode(stopp_job)

# create all the LSCdataFind jobs to run in sequence
prev_df = None

# loop through the analysis data
for seg in data:
  # find the data with LSCdataFind
  df1 = stochastic.LSCDataFindNode(df_job)
  df1.set_start(seg.start())
  df1.set_end(seg.end())
  df1.set_observatory(ifo1[0])
  df1.set_type(type1)
  df1.set_server(server)
  if prev_df:
    df1.add_parent(prev_df)
  prev_df = df1

  # find data with LSCdataFind for different detector
  if ifo1[0] != ifo2[0]:
    df2 = stochastic.LSCDataFindNode(df_job)
    df2.set_start(seg.start())
    df2.set_end(seg.end())
    df2.set_observatory(ifo2[0])
    df2.set_type(type2)
    df2.set_server(server)
    df2.add_parent(prev_df)
    prev_df = df2

  # add LSCdataFind node(s) if required
  if do_datafind:
    dag.add_node(df1)
    if ifo1[0] != ifo2[0]:
      dag.add_node(df2)

  # if same site, only one LSCdataFind job is run
  if ifo1[0] == ifo2[0]:
    df2 = df1

  # set up stochastic job
  for chunk in seg:
    stoch = stochastic.StochasticNode(stoch_job)
    stoch.set_start(chunk.start())
    stoch.set_end(chunk.end())
    stoch.set_ifo_one(ifo1)
    stoch.set_ifo_two(ifo2)
    stoch.set_cache_one(df1.get_output())
    stoch.set_cache_two(df2.get_output())
    stoch.set_calibration_one(ifo1,chunk.start())
    stoch.set_calibration_two(ifo2,chunk.start())

    # add stochastic output file to stopp job
    stopp.add_var_arg(stoch.get_output())

    # add dependancy for stopp jobs on stochastic jobs
    if do_stochastic:
      stopp.add_parent(stoch)

    # add dependancy on LSCdataFind jobs, if required
    if do_datafind:
      stoch.add_parent(df1)
      if ifo1[0] != ifo2[0]:
        stoch.add_parent(df2)

    # add stochastic node
    if do_stochastic:
      dag.add_node(stoch)

# add stopp node
if do_stopp:
  dag.add_node(stopp)

# write out the DAG
dag.write_sub_files()
dag.write_dag()

# write a message telling the user that the DAG has been written
print "Created DAG file", dag.get_dag_file()
if do_datafind:
  print """\nAs you are running LSCdataFind jobs, do not forget to initialise your
grid proxy certificate on the condor submit machine by running the
commands

  $ unset X509_USER_PROXY
  $ grid-proxy-init -hours 72

Enter your pass phrase when promted. The proxy will be valid for 72 hours.
If you expect the LSCdataFind jobs to take longer to complete, increase
the time specified in the -hours option to grid-proxy-init. You can check
that the grid proxy has been sucessfully created by executing the command:

  $ grid-cert-info -all -file /tmp/x509up_u`id -u`

This will also give the expiry time of the proxy."""

print """\nThe DAG can now be summitted by executing the following command
on a condor submit machine

  $ condor_submit_dag""", dag.get_dag_file()

print """\nIt is currently recommended that you pass condor_submit_dag the option
-maxjobs 30, to limit the maximum of jobs. This is due to condor being
designed to run jobs that take several hours to complete, whereas the
stochastic jobs are complete within a few minutes. Running the stochastic
pipeline on a cluster without the -maxjobs option will essentially bring
the cluster to a standstill."""

# write out a logfile
if usertag:
  log_fh = open(basename + '.pipeline.' + usertag + '.log', 'w')
else:
  log_fh = open(basename + '.pipeline.log', 'w')

# FIXME: the following code uses obsolete CVS ID tags.
# It should be modified to use git version information.
log_fh.write("$Id$" + "\n\n")
log_fh.write("Invoked with arguments:\n")
for o, a in opts:
  log_fh.write(o + ' ' + a + '\n')

log_fh.write("\nConfig file has the CVS strings:\n")
log_fh.write(cp.get('pipeline','version') + "\n\n")

log_fh.write("Parsed " + str(len(data)) + " science segments\n")
total_data = 0
for seg in data:
  for chunk in seg:
    total_data += len(chunk)
print >> log_fh, "total data =", total_data

print >> log_fh, "\n", data
for seg in data:
  print >> log_fh, seg
  for chunk in seg:
    print >> log_fh, chunk

sys.exit(0)

# vim: et syntax=python
