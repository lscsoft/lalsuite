"""
noise_comp_pipe - standalone noise comparison pipeline driver script

This script produces the condor submit and dag files to run
the noise comparison between h(t) and calibrated DARM_ERR
"""

__author__ = 'Xavier Siemens<siemens@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'


# import standard modules
import sys, os, shutil
import optparse
from optparse import *
import getopt, re, string
import tempfile
import ConfigParser
from pylab import *

# import the modules we need to build the pipeline
from glue import pipeline, lal
import strain


# MAIN PROGRAM
########## Option Parser ######################################################

usage = """usage: %prog [options]"""
parser = OptionParser( usage )

parser.add_option("-l", "--log-path",action="store",type="string",\
    metavar=" PATH",help="directory to write condor log file")

parser.add_option("-i", "--config-file",action="store",type="string",\
    metavar=" FILE",help="ini file")

parser.add_option("-f", "--basename",action="store",type="string",\
    metavar=" STR",help="basename for .dag file (excluding the .dag)")

parser.add_option("-S", "--segment-filename",action="store",type="string",\
    metavar=" FILE",help="segments file")

parser.add_option("-w", "--write-script",action="store_true",default=False,\
    help="write out a shell script workflow DEFAULT: FALSE")

parser.add_option("-d", "--data-find",action="store_true",default=False,\
    help="run data find portions of the pipeline DEFAULT: FALSE")

parser.add_option("-n", "--noise-comp",action="store_true",default=False,\
    help="run noise comp portions of the pipeline DEFAULT: FALSE")

parser.add_option("-c", "--check-datafind-jobs",action="store_true",default=False,\
    help="check to see if previous datafind caches are okay, and don't bother with those jobs DEFAULT: FALSE")

parser.add_option("-x", "--write-dax",action="store_true",default=False,\
    help="write out a dax workflow DEFAULT: FALSE")

parser.add_option("-C", "--cat-noise-jobs",action="store_true",default=False,\
    help="Cat together the output of noise jobs (run this after they complete)")

parser.add_option("-p", "--plot",action="store_true",default=False,\
    help="Plot the output of noise jobs (run this after they complete)")

parser.add_option("-s", "--plot-systematics",action="store_true",default=False,\
    help="plot the systematic errors")

parser.add_option("-V", "--veto-list",action="store_true",default=False,\
    help="Write a veto list of bad times")

# parse the command line
(opts,args) = parser.parse_args()


# Should this be an option to set in the INI FILE ?
df_pad=128

if not opts.config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not opts.segment_filename:
  print >> sys.stderr, "No segment filename specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not opts.basename:
  print >> sys.stderr, "No dag file base name specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

# create the config parser object and read in the ini file
cp = ConfigParser.ConfigParser()
cp.read(opts.config_file)

# create a log file that the Condor jobs will write to
tempfile.tempdir = opts.log_path
tempfile.template = opts.basename + '.log'
logfile = tempfile.mktemp()
fh = open( logfile, "w" )
fh.close()

# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile,opts.write_dax)
dag.set_dag_file(opts.basename)

subsuffix = '.sub'

# create the Condor jobs that will be used in the DAG
mkdir_job = strain.MkdirJob('logs',cp)
mkdir_node = strain.MkdirNode(mkdir_job,'cache')
if opts.write_dax: dag.add_node(mkdir_node)

# try and make a directory to store the cache files and job logs
try: os.mkdir('logs')
except: pass
#
#try: os.mkdir('cache')
#except: pass

df_job = pipeline.LSCDataFindJob('cache','logs',cp,opts.write_dax)
noise_job = strain.NoiseJob(cp,opts.write_dax)
qscan_job = strain.qscanJob(cp)

# submit files
df_job.set_sub_file( opts.basename + '.datafind'+ subsuffix )
noise_job.set_sub_file( opts.basename + '.noisecomp' + subsuffix )

# get the pad and chunk lengths from the values in the ini file
length = int(cp.get('pipeline', 'segment-length'))

# get the ifo to filter
ifo = cp.get('pipeline','ifo')
datatype_hoft = cp.get('input','type-hoft')
datatype_derr = cp.get('input','type-derr')

# read in data that is specific to calibration epochs
epochs = strain.EpochData(cp,opts)

epoch_cnt = 0;

# loop over the segments defined by the calibration epochs
print "\n"
for epoch in epochs.epoch_segs():
  noise_output_files = []
  noise_output_files2 = []
  print "setting up jobs for calibration epoch " + str(epoch[1])+" - "+str(epoch[2]) + "..."
  #output the epochs in their own directories
  epoch_dir = 'EPOCH'+'-'+str(epoch[1])+'-'+str(epoch[2])
  mkdir_node2 = strain.MkdirNode(mkdir_job,epoch_dir)
  if opts.write_dax: dag.add_node(mkdir_node2)
  if opts.cat_noise_jobs: catfile = strain.open_noise_cat_file(epoch_dir)
  # Make a ScienceData class for the calibration epochs
  epoch_data = pipeline.ScienceData()
  epoch_data.append_from_tuple(epoch)
  # read science segs that are greater or equal to a chunk from the input file
  data = pipeline.ScienceData()
  data.read(opts.segment_filename,0)
  # intersect the science segments with the calibration epoch
  data.intersection(epoch_data)
  # create the chunks from the science segments
  data.make_chunks(length,0,0,0,0)
  data.make_short_chunks_from_unused(0,0,0,0,0)

  # create all the LSCdataFind jobs to run in sequence
  prev_df1 = None
  prev_df2 = None
  # only do data find jobs if requested
  # find all the h(t) data
  df1 = pipeline.LSCDataFindNode(df_job)
  df1.set_start(int(epoch[1])-df_pad)
  df1.set_end(int(epoch[2]) +df_pad)
  df1.set_observatory(ifo[0])
  df1.set_type(datatype_hoft)
  df1.set_name("df1_"+ifo+"_"+str(epoch_cnt))
  # see if the cache files laying around are still okay
  if opts.check_datafind_jobs:
    try: df1cache = lal.Cache.fromfile(open(df1.get_output(),'r'))
    except: df1cache = None
    if df1cache: found,missed = df1cache.checkfilesexist("ignore")
    else: missed = True
  else: missed = True

  if opts.data_find and missed and opts.write_dax: df1.add_parent(mkdir_node)
  if prev_df1 and opts.data_find and missed:
    df1.add_parent(prev_df1)
  if opts.data_find and missed: dag.add_node(df1)
  prev_df1 = df1

  # find all the DARM_ERR data
  df2 = pipeline.LSCDataFindNode(df_job)
  df2.set_start(int(epoch[1])-df_pad)
  df2.set_end(int(epoch[2])+df_pad)
  df2.set_observatory(ifo[0])
  df2.set_type(datatype_derr)
  df2.set_name("df2_"+ifo+"_"+str(epoch_cnt))
  #df2.add_var_opt("no-proxy"," ")
  #df2.add_output_file(df2.get_output_cache())
  # see if the cache files laying around are still okay
  if opts.check_datafind_jobs:
    try: df2cache = lal.Cache.fromfile(open(df2.get_output(),'r'))
    except: df2cache = None
    if df2cache: found,missed = df2cache.checkfilesexist("ignore")
    else: missed = True
  else: missed = True
  if opts.data_find and missed and opts.write_dax: df2.add_parent(mkdir_node)
  if prev_df2 and opts.data_find and missed:
    df2.add_parent(prev_df2)
  if opts.data_find and missed: dag.add_node(df2)
  prev_df2 = df2

  segment_no = -1
  for seg in data:
    segment_no = segment_no + 1

    #noise jobs
    chunk_number=-1
    # only add noise jobs if requested
    for chunk in seg:
      chunk_number=chunk_number+1
      # throw out small jobs
      if chunk.end()-chunk.start() < int(cp.get('noisecomp','time')):
        continue
      #make the directory where the data's going to go
      gps_str=str(chunk.start())

      #Noise job for first ifo
      noise1 = strain.NoiseNode(noise_job)
      noise1.set_start(chunk.start())
      noise1.set_end(chunk.end())
      noise1.add_var_opt('olg-re',epochs.epoch_data[epoch_cnt][4])
      noise1.add_var_opt('olg-im',epochs.epoch_data[epoch_cnt][5])
      noise1.add_var_opt('servo-re',epochs.epoch_data[epoch_cnt][6])
      noise1.add_var_opt('servo-im',epochs.epoch_data[epoch_cnt][7])
      noise1.add_var_opt('whitener-re',epochs.epoch_data[epoch_cnt][8])
      noise1.add_var_opt('whitener-im',epochs.epoch_data[epoch_cnt][9])
      noise1.add_var_opt('olg-file',epochs.epoch_data[epoch_cnt][11])
      noise1.add_input_file(epochs.epoch_data[epoch_cnt][11])
      noise1.add_var_opt('sensing-file',epochs.epoch_data[epoch_cnt][12])
      noise1.add_input_file(epochs.epoch_data[epoch_cnt][12])
      noise1.add_cache(df1.get_output(),'hoft-cache')
      noise1.add_cache(df2.get_output(),'derr-cache')
      noise1.add_input_file(cp.get('noisecomp','freq-file'))
      noise1.add_var_opt('output-file',epoch_dir+'/out-'+gps_str+'-'+str(len(chunk))+'.txt')
      noise1.add_output_file(epoch_dir+'/out-'+gps_str+'-'+str(len(chunk))+'.txt')
      noise1.add_var_opt('output-bin',epoch_dir+'/out-bin-'+gps_str+'-'+str(len(chunk))+'.txt')
      noise1.add_output_file(epoch_dir+'/out-bin-'+gps_str+'-'+str(len(chunk))+'.txt')
      noise1.set_name("noise_"+ifo+"_"+str(epoch_cnt)+str(segment_no)+"_"+str(chunk_number))
      # only add data find jobs as parents if requested on command line
      if opts.write_dax: noise1.add_parent(mkdir_node2)
      if opts.data_find:
        noise1.add_parent(df1)
        noise1.add_parent(df2)
      if opts.noise_comp: dag.add_node(noise1)
      noise_output_files.extend(noise1.get_output_files())
      if opts.cat_noise_jobs: strain.cat_noise_jobs(catfile,noise1)
  #increment the epoch counter
  if opts.cat_noise_jobs: catfile.close()
  if opts.plot_systematics: strain.plot_systematics(noise_output_files,cp,epoch_dir,epochs.epoch_data[epoch_cnt],dag,opts)
  if opts.plot or opts.veto_list: strain.plot_noise_jobs(noise_output_files,cp,epoch_dir,epochs.epoch_data[epoch_cnt],dag,qscan_job,opts)
  noise_output_files = []
  epoch_cnt+=1

# write out the DAG
if not opts.cat_noise_jobs:
  dag.write_sub_files()
  dag.write_dag()
  if opts.write_dax: dag.write_pegasus_rls_cache(cp.get("ldgsubmitdax","gsiftp"),cp.get("ldgsubmitdax","pool"))
  if opts.write_script: dag.write_script()

  print "\nDAG contains " + str(len(dag.get_nodes())) + " nodes.\n"

  # write out a log file for this script
  log_fh = open(opts.basename + '.pipeline.log', 'w')

  # FIXME: the following code uses obsolete CVS ID tags.
  # It should be modified to use git version information.
  log_fh.write( "$Id$" + "\n\n" )
  log_fh.write( "Invoked with arguments:\n" )

  classattrs = dir(optparse.Values)
  for name in [_ for _ in dir(opts) if _ not in classattrs]:
      log_fh.write("--"+name.replace('_','-')+"="+ str(getattr(opts, name))+"\n")


  log_fh.write( "\n" )
  log_fh.write( "Parsed " + str(len(data)) + " science segments\n" )
  total_data = 0
  for seg in data:
    for chunk in seg:
      total_data += len(chunk)
  print >> log_fh, "total data =", total_data

  print >> log_fh, "\n===========================================\n"
  print >> log_fh, data
  for seg in data:
    print >> log_fh, seg
    for chunk in seg:
      print >> log_fh, chunk, 'length', int(chunk.end())-int(chunk.start())
      endgps=chunk.end()

if not opts.cat_noise_jobs:
  # write a message telling the user that the DAG has been written
  print "\nCreated a DAG file which can be submitted by executing"
  print "\n   condor_submit_dag", dag.get_dag_file()
  print """\nfrom a condor submit machine (e.g. hydra.phys.uwm.edu)\n
  If you are running LSCdataFind jobs, do not forget to initialize your grid
  proxy certificate on the condor submit machine by running the commands

  unset X509_USER_PROXY
  grid-proxy-init -hours 72

  Enter your pass phrase when promted. The proxy will be valid for 72 hours.
  If you expect the LSCdataFind jobs to take longer to complete, increase the
  time specified in the -hours option to grid-proxy-init. You can check that
  the grid proxy has been sucessfully created by executing the command:

  grid-cert-info -all -file /tmp/x509up_u`id -u`

  This will also give the expiry time of the proxy. You should also make sure
  that the environment variable LSC_DATAFIND_SERVER is set the hostname and
  optional port of server to query. For example on the UWM medusa cluster this
  you should use

  export LSC_DATAFIND_SERVER=dataserver.phys.uwm.edu

  Contact the administrator of your cluster to find the hostname and port of the
  LSCdataFind server.
  """

sys.exit(0)

