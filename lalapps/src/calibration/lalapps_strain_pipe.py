"""
strain_pipe - standalone strain pipeline driver script

This script produces the condor submit and dag files to run
the standalone strain
"""

__author__ = 'Xavier Siemens<siemens@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'


# import standard modules
import sys, os, shutil
import getopt, re, string
import tempfile
import ConfigParser

# import the modules we need to build the pipeline
from glue import pipeline
import strain

def usage():
  msg = """\
Usage: lalapps_ring_pipe [options]

  -h, --help               display this message
  -v, --version            print version information and exit
  -s, --gps-start-time     GPS start time
  -e, --gps-end-time       GPS end time
  -S, --segment-file       segments file
  -a,  IGNORED
  -b,  IGNORED
  -f, --dag-file           basename for .dag file (excluding the .dag)
  -t, --aux-path           path to auxiliary files
"""
  print >> sys.stderr, msg

# pasrse the command line options to figure out what we should do
shortop = "hv:s:e:S:f:t:a:b:"
longop = [
  "help",
  "version",
  "gps-start-time=",
  "gps-end-time=",
  "segment-file=",
  "dag-file=",
  "aux-path=",
  "trig-start",
  "trig-end"
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

config_file = None
log_path = None
GPSStart = None
GPSEnd = None
basename = None
aux_path = None
segment_filename = None

for o, a in opts:
  if o in ("-v", "--version"):
    print "$Id$"
    sys.exit(0)
  elif o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-s", "--gps-start-time"):
    GPSStart = a
  elif o in ("-e", "--gps-end-time"):
    GPSEnd = a
  elif o in ("-S", "--segment-file"):
    segment_filename = a
  elif o in ("-f", "--dag-file"):
    basename = os.path.splitext(a)[0]
  elif o in ("-t", "--aux-path"):
    aux_path = a
  elif o in ("-a", "--trig-start"):
    pass
  elif o in ("-b", "--trig-end"):
    pass
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

config_file = aux_path + '/strainDAGH1.ini'
log_path = '/usr1/xsiemens/'
df_pad=128

if not config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not log_path:
  print >> sys.stderr, "No log file path specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (not GPSStart or not GPSEnd) and not segment_filename:
  print >> sys.stderr, "No GPS start time and end times or segment filename specified."
  print >> sys.stderr, "Either GPS start time and end times, or a segment filename must be specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not basename:
  print >> sys.stderr, "No dag file base name specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not aux_path:
  print >> sys.stderr, "No auxiliary file path specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

# try and make a directory to store the cache files and job logs
try: os.mkdir('logs')
except: pass

try: os.mkdir('cache')
except: pass

# here determone whether I am running using a segment file or not
running_online = 0
if not segment_filename:
  running_online = 1

# create the config parser object and read in the ini file
cp = ConfigParser.ConfigParser()
cp.read(config_file)

# create a log file that the Condor jobs will write to
tempfile.tempdir = log_path
tempfile.template = basename + '.log'
logfile = tempfile.mktemp()
fh = open( logfile, "w" )
fh.close()

# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile)
dag.set_dag_file(basename)

# create the Condor jobs that will be used in the DAG
df_job = pipeline.LSCDataFindJob('cache','logs',cp)
strain_job = strain.StrainJob(cp)

# submit files
subsuffix = '.sub'
df_job.set_sub_file( basename + '.datafind'+ subsuffix )
strain_job.set_sub_file( basename + '.strain' + subsuffix )

# if runnign on-line make segments filename
if running_online:
  segment_file=open('strain_segment.txt',mode='w')
  print >> segment_file, '1',GPSStart,' ',GPSEnd,' ',int(GPSEnd)-int(GPSStart)
  segment_file.close()
  segment_filename='strain_segment.txt'

# get the pad and chunk lengths from the values in the ini file
length = int(cp.get('pipeline', 'segment-length'))
overlap = int(cp.get('strain', 'wings'))

# read science segs that are greater or equal to a chunk from the input file
data = pipeline.ScienceData()
data.read(segment_filename,0)

# create the chunks from the science segments
data.make_chunks(length+2*overlap,2*overlap,0,0,0)
data.make_short_chunks_from_unused(0,2*overlap,0,0,0)
#make_short_chunks_from_unused(self, min_length, overlap=0, play=0, sl=0, excl_play=0)
#I don't really know why the min length needs to be set to 0 but it does.

# get the ifo to filter
ifo = cp.get('pipeline','ifo')
datatype = cp.get('input','type')
base_data_dirL1 = cp.get('pipeline','data-dirL1')
base_data_dirL2 = cp.get('pipeline','data-dirL2')
frametype = cp.get('strain','frame-type')
frametypeL1=frametype+"_L1"
frametypeL2=frametype+"_L2"

#open a file that will contain all the names of frame files to be generated
framelist_file=open('frame_files.txt',mode='w')
#open a file that when run under bash will check all the cache files
cachecheck_file=open('cache_check.txt',mode='w')

# create all the LSCdataFind jobs to run in sequence
prev_df = None
segment_no = -1
for seg in data:
  segment_no = segment_no + 1
  # find all the data for first ifo
  df = pipeline.LSCDataFindNode(df_job)
  df.set_start(str( int(seg.start())-df_pad ) )
  df.set_end(str( int(seg.end())+df_pad )  )
  df.set_observatory(ifo[0])
  df.set_type(datatype)
  df.set_name("df_"+ifo+"_"+str(segment_no))

  command = "/archive/home/xsiemens/lscsoft/glue/bin/LSCdataFindcheck --gps-start-time "+str(seg.start())+\
  " --gps-end-time "+str(seg.end())+" "+df.get_output()
  print >> cachecheck_file, command

  if prev_df:
    df.add_parent(prev_df)
  dag.add_node(df)
  prev_df = df

  #strain jobs
  chunk_number=-1
  for chunk in seg:
    chunk_number=chunk_number+1

    #make the directory where the data's going to go
    gps_str=str(chunk.start())
    gps_time_first_four=gps_str[0]+gps_str[1]+gps_str[2]+gps_str[3]
    try: os.mkdir(base_data_dirL1+'/'+ifo[0]+'-'+frametypeL1+'-'+gps_time_first_four)
    except OSError, err:
      import errno
      #print "Warning:", err
    try: os.mkdir(base_data_dirL2+'/'+ifo[0]+'-'+frametypeL2+'-'+gps_time_first_four)
    except OSError, err:
      import errno
      #print "Warning:", err

    #Strain job for first ifo
    strain1 = strain.StrainNode(strain_job)
    strain1.set_start(chunk.start())
    strain1.set_end(chunk.end())
    strain1.set_cache(df.get_output())
    strain1.set_name("strain_"+ifo+"_"+str(segment_no)+"_"+str(chunk_number))
    strain1.add_var_opt('data-dirL1',base_data_dirL1+'/'+ifo[0]+'-'+frametypeL1+'-'+gps_time_first_four+'/')
    strain1.add_var_opt('data-dirL2',base_data_dirL2+'/'+ifo[0]+'-'+frametypeL2+'-'+gps_time_first_four+'/')
    strain1.add_parent(df)
    dag.add_node(strain1)

    print >> framelist_file, 'ls '+directory+'/'+ifo[0] \
    +'-'+frametype+'-'+str(int(chunk.start())+int(overlap))+'-' \
    +str(int(chunk.end())-int(chunk.start())-2*int(overlap))+'.gwf > /dev/null'

cachecheck_file.close()
framelist_file.close()

# write out the DAG
dag.write_sub_files()
dag.write_dag()

# write out a log file for this script
log_fh = open(basename + '.pipeline.log', 'w')

# FIXME: the following code uses obsolete CVS ID tags.
# It should be modified to use git version information.
log_fh.write( "$Id$" + "\n\n" )
log_fh.write( "Invoked with arguments:\n" )
for o, a in opts:
  log_fh.write( o + ' ' + a + '\n' )
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

if running_online:
  print >> sys.stdout, seg.start()+overlap,int(chunk.end())-overlap

if not running_online:
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

