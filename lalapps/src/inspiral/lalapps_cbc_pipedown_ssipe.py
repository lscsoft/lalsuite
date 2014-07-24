##############################################################################
# import standard modules
import os, sys, re, copy
from optparse import OptionParser
import tempfile

##############################################################################
# import the modules we need to build the pipeline
from glue import pipeline
from glue import lal
from glue import iterutils
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import process
from glue import git_version
import inspiral
from lalapps import inspiralutils

__prog__ = 'cbc_pipedown'
__author__ = 'Collin Capano <cdcapano@physics.syr.edu>'

description = \
"Program to construct post-processing dag."

##############################################################################
# Function Definitions

def get_veto_cat_from_tag( tag ):
  """
  Returns the veto category number in a tag.
  Assumes tag has the form _CAT_N_VETO for N>1 and
  nothing for N=1.
  """
  if 'VETO' in tag:
    cat_num = int(tag.split('_')[-2])
  else:
    # is category 1 veto
    cat_num = 1

  return cat_num

def get_veto_segments_name( veto_cat_num, cumulative = True ):
  """
  Given a category number, returns a veto segments name 
  as set by segs_from_cats.

  @veto_cat_num: integer representing the category veto
  @cumulative: If set to True, will add CUMULATIVE to the name.
  """
  if cumulative:
    return ''.join([ 'VETO_CAT', str(veto_cat_num), '_CUMULATIVE' ])
  else:
    return ''.join([ 'VETO_CAT', str(veto_cat_num) ])
    

##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################

##############################################################################
# parse command-line arguments
parser = OptionParser( 
  version = git_version.verbose_msg,
  usage   = "%prog [options]",
  description = description
  )

parser.add_option( "", "--ihope-cache", action = "store", type = "string",
  default = None,
  help =
    "The ihope cache to read."
  )
parser.add_option( "", "--config-file", action = "store", type = "string",
  default = None,
  help =
    "The .ini file to use."
  )
parser.add_option( "", "--log-path", action = "store", type = "string",
  default = None,
  help =
    "Directory to write condor log file and perform SQLite operation in. " +
    "Should be a local directory."
  )
parser.add_option( "", "--gps-start-time", action = "store", type = "string",
  default = None,
  help =
    "GPS start time of the ihope run."
  )
parser.add_option( "", "--gps-end-time", action = "store", type = "string",
  default = None,
  help =
    "GPS end time of the ihope run."
  )
parser.add_option( "", "--generate-all-data-plots", action = "store_true",
  default = False,
  help =
    "Turn on if want to open the box. Otherwise, only plots of playground " +
    "data will be made. WARNING: Even if this option is off, all_data and " +
    "and exclude_play coincs will still exist in the resulting databases " +
    "(they just won't be plotted)."
  )

(options, args) = parser.parse_args()

##############################################################################
# Sanity check of input arguments
if not options.ihope_cache:
  raise ValueError, "An ihope-cache file is required."
if not options.config_file:
  raise ValueError, "A config-file is required."
if not options.log_path:
  raise ValueError, "A log-path is required."

##############################################################################
# Create log file
logdoc = ligolw.Document()
logdoc.appendChild(ligolw.LIGO_LW())
proc_id = process.register_to_xmldoc(logdoc, __prog__, options.__dict__, version = git_version.id)

##############################################################################
# parse the ini file and initialize
cp = pipeline.DeepCopyableConfigParser()
cp.read(options.config_file)

tmp_space = cp.get('pipeline', 'node-tmp-dir')

experiment_start = options.gps_start_time
experiment_duration = str(int(options.gps_end_time) - int(options.gps_start_time))

try:
  do_cats = cp.get('pipeline', 'pipedown-cats')
  do_cats = do_cats.split(',')
except:
  do_cats = False

# if logs directory not present, create it
try:
  os.mkdir('logs/')
except:
  pass

##############################################################################
# create a log file that the Condor jobs will write to
basename = re.sub(r'\.ini', r'', options.config_file)

tempfile.tempdir = options.log_path
tempfile.template = '.'.join([ basename, 'dag.log.' ])
logfile = tempfile.mktemp()
fh = open( logfile, "w" )
fh.close()

##############################################################################
# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile)
dag.set_dag_file(basename)

##############################################################################
# Open the ihope cache and create THINCA cache

print "Parsing the ihope cache..."

coinc_tag = cp.get('pipeline', 'coinc-file-tag')
ihope_cache = [line for line in file(options.ihope_cache) \
  if coinc_tag in line or " INJECTIONS" in line]

thinca_cache = lal.Cache([lal.CacheEntry(entry) for entry in ihope_cache \
  if coinc_tag in entry])
inj_cache = lal.Cache([lal.CacheEntry(entry) for entry in ihope_cache if \
  " INJECTIONS" in entry])

del ihope_cache

# get the USERTAGS from the thinca_cache
# for single stage runs with ssipe, the thinca's output is of the form
# IFOs_THINCA_UserTag_StartTime_Duration.xml.gz
# where UserTag = TiSiNum_RunName_CAT_X_VETO

skip_tags = 2
user_tags = set([ '_'.join([ entry.description.split('_')[ii] 
  for ii in range( skip_tags, len(entry.description.split('_')) ) ])
  for entry in thinca_cache ])

# only do the cats requested
if do_cats:
  user_tags = set([tag for tag in user_tags if str(get_veto_cat_from_tag(tag)) in do_cats])

# create a list to store non-injection databases for later caching
non_sim_dbs = []

##############################################################################
# Get the time column name (for inspiral, it is end_time; for ringdown, it is 
# start_time)

print "Getting the name of the time column..."
time_column = cp.get('pipeline', 'time-column')

##############################################################################
# Create the needed jobs to be run

# ligolw_sqlite jobs; there are 4 different types:
# 1. sql_replace_job: reads from a cache and writes to a database
#   with the --replace option set (only for adding thinca_to_coinc files )
# 2. sql_fromcache_job: reads from a cache and writes to a database without
#   the replace option
# 3. sql_frominput_job: adds a single file to a databse (for veto xml files)
# 4. sql_extract_job: extracts xmls from a database (for turning simulation
#   databases into xml files)
sql_replace_job = pipeline.LigolwSqliteJob(cp)
sql_replace_job.set_replace()
sql_fromcache_job = pipeline.LigolwSqliteJob(cp)
sql_frominput_job = pipeline.LigolwSqliteJob(cp)
sql_extract_job = pipeline.LigolwSqliteJob(cp)

# for plotting jobs with plot_playground_only option, only create jobs without
# if generate-all-data-plots set
plotslides_play_job = inspiral.PlotSlidesJob(cp)
plotslides_play_job.set_plot_playground_only()
plotcumhist_play_job = inspiral.PlotCumhistJob(cp)
plotcumhist_play_job.set_plot_playground_only()
if options.generate_all_data_plots:
  plotslides_job = inspiral.PlotSlidesJob(cp)
  plotcumhist_job = inspiral.PlotCumhistJob(cp)

# since some files won't be passed through minifollowups, need two sets of
# printsims and printmissed jobs, one with columns, one without
printsims_job = inspiral.LigolwCBCPrintJob(cp, 'printsims', ['cbc_print', 'printsims'])
printsims_nominifup_job = inspiral.LigolwCBCPrintJob(cp, 'printsims', ['cbc_print', 'printsims'])
printmissed_job = inspiral.LigolwCBCPrintJob(cp, 'printmissed', ['cbc_print', 'printmissed'])
printmissed_nominifup_job = inspiral.LigolwCBCPrintJob(cp, 'printmissed', ['cbc_print', 'printmissed'])
# time-slide minifollowup jobs require a separate time-slides argument
# and missed injection minifollowups don't use full xml files
minifup_job = inspiral.MiniFollowupsJob(cp)
minifup_ts_job = inspiral.MiniFollowupsJob(cp)
minifup_ts_job.set_time_slides()
minifup_cm_job = inspiral.MiniFollowupsJob(cp)

# the number of and options given to plotfm job is optional and is set by the ini file
# to handle this, we do a similar thing as is done in ihope with inspinj: the 'plot_found/missed'
# section in the ini file gives the names of sections to read for plotfm jobs; these sections are cycled over
# and temporarily re-named to 'plotfm' then passed to the PlotFMJob class to create the job
plotfm_jobs = []
if cp.has_section('plotfm'):
  for plotfm_tag in cp.options('plotfm'):
    # create a copy of the config parser, blank out the plotfm section, and replace
    # with a new plotfm section with the options taken from plotfm_tag section
    fmcp = copy.deepcopy(cp)
    fmcp.remove_section('plotfm')
    fmcp.add_section('plotfm')
    for fmopt, fmval in cp.items(plotfm_tag):
      fmcp.set('plotfm', fmopt, fmval)
    plotfm_jobs.append( (plotfm_tag, inspiral.PlotFMJob(fmcp)) )

# ditto repop coinc
repop_coinc_jobs = []
if cp.has_section('repop_coinc'):
  for repop_coinc_tag in cp.options('repop_coinc'):
    repopcp = copy.deepcopy(cp)
    repopcp.remove_section('repop_coinc')
    repopcp.add_section('repop_coinc')
    for repopopt, repopval in cp.items(repop_coinc_tag):
      repopcp.set('repop_coinc', repopopt, repopval)
    repop_coinc_jobs.append( (repop_coinc_tag, inspiral.RepopCoincJob(repopcp)) )

# ditto dbinjfind
dbinjfind_jobs = []
if cp.has_section('dbinjfind'):
  for dbinjfind_tag in cp.options('dbinjfind'):
    dbinjcp = copy.deepcopy(cp)
    dbinjcp.remove_section('dbinjfind')
    dbinjcp.add_section('dbinjfind')
    for dbinjopt, dbinjval in cp.items(dbinjfind_tag):
      dbinjcp.set('dbinjfind', dbinjopt, dbinjval)
    dbinjfind_jobs.append( (dbinjfind_tag, inspiral.DBInjFindJob(dbinjcp)) )

# ditto cluster_coincs
cluster_coincs_jobs = []
if cp.has_section('cluster_coincs'):
  for cluster_coincs_tag in cp.options('cluster_coincs'):
    clustercp = copy.deepcopy(cp)
    clustercp.remove_section('cluster_coincs')
    clustercp.add_section('cluster_coincs')
    for clusteropt, clusterval in cp.items(cluster_coincs_tag):
      clustercp.set('cluster_coincs', clusteropt, clusterval)
    cluster_coincs_jobs.append( (cluster_coincs_tag, inspiral.ClusterCoincsJob(clustercp)) )

# all other jobs: all the nodes share the same static options for these
dbsimplify_job = inspiral.DBSimplifyJob(cp)
dbaddinj_job = inspiral.DBAddInjJob(cp)
injfind_job = inspiral.InjFindJob(cp)
comp_durs_job = inspiral.ComputeDurationsJob(cp)
ucfar_job = inspiral.CFarJob(cp, ['cfar-uncombined'])
ccfar_job = inspiral.CFarJob(cp, ['cfar-combined'])
printlc_job = inspiral.LigolwCBCPrintJob(cp, 'printlc', ['cbc_print', 'printlc'])
plotifar_job = inspiral.PlotIfarJob(cp)
search_volume_job = inspiral.SearchVolumeJob(cp)
search_upper_limit_job = inspiral.SearchUpperLimitJob(cp)


# set submit file names
sql_replace_job.set_sub_file( '.'.join([ basename, 'replace_from_cache', 'ligolw_sqlite', 'sub' ]) )
sql_fromcache_job.set_sub_file( '.'.join([ basename, 'add_from_cache', 'ligolw_sqlite', 'sub' ]) )
sql_frominput_job.set_sub_file( '.'.join([ basename, 'add_from_input', 'ligolw_sqlite', 'sub' ]) )
sql_extract_job.set_sub_file( '.'.join([ basename, 'extract', 'ligolw_sqlite', 'sub' ]) )
dbsimplify_job.set_sub_file( '.'.join([ basename, 'dbsimplify', 'sub' ]) )
dbaddinj_job.set_sub_file( '.'.join([ basename, 'dbaddinj', 'sub' ]) )
injfind_job.set_sub_file( '.'.join([ basename, 'injfind', 'sub' ]) )
comp_durs_job.set_sub_file( '.'.join([ basename, 'compute_durations', 'sub' ]) )
for repoptag, repopjob in repop_coinc_jobs:
  repopjob.set_sub_file( '.'.join([ basename, repoptag, 'repop_coinc', 'sub' ]) )
for dbinjtag, dbinjjob in dbinjfind_jobs:
  dbinjjob.set_sub_file( '.'.join([ basename, dbinjtag, 'dbinjfind', 'sub' ]) )
for clustertag, clusterjob in cluster_coincs_jobs:
  clusterjob.set_sub_file( '.'.join([ basename, clustertag, 'cluster_coincs', 'sub' ]) )
ucfar_job.set_sub_file( '.'.join([ basename, 'uncombined', 'cfar', 'sub' ]) )
ccfar_job.set_sub_file( '.'.join([ basename, 'combined', 'cfar', 'sub' ]) )
printlc_job.set_sub_file( '.'.join([ basename, 'printlc', 'sub' ]) )
printsims_job.set_sub_file( '.'.join([ basename, 'printsims', 'sub' ]) )
printsims_nominifup_job.set_sub_file( '.'.join([ basename, 'no_minifollowups', 'printsims', 'sub' ]) )
printmissed_job.set_sub_file( '.'.join([ basename, 'printmissed', 'sub' ]) )
printmissed_nominifup_job.set_sub_file( '.'.join([ basename, 'no_minifollowups', 'printmissed', 'sub' ]) )
minifup_job.set_sub_file( '.'.join([ basename, 'minifollowups', 'sub' ]) )
minifup_ts_job.set_sub_file( '.'.join([ basename, 'time_slides', 'minifollowups', 'sub' ]) )
minifup_cm_job.set_sub_file( '.'.join([ basename, 'closest_missed', 'minifollowups', 'sub' ]) )
plotifar_job.set_sub_file( '.'.join([ basename, 'plotifar', 'sub' ]) )
plotslides_play_job.set_sub_file( '.'.join([ basename, 'playground_only', 'plotslides', 'sub' ]) )
plotcumhist_play_job.set_sub_file( '.'.join([ basename, 'playground_only', 'plotcumhist', 'sub' ]) )
search_volume_job.set_sub_file( '.'.join([ basename, 'search_volume', 'sub' ]) )
search_upper_limit_job.set_sub_file( '.'.join([ basename, 'search_upper_limit', 'sub' ]) )
for fmtag, fmjob in plotfm_jobs:
  fmjob.set_sub_file( '.'.join([ basename, fmtag, 'plotfm', 'sub' ]) )
if options.generate_all_data_plots:
  plotslides_job.set_sub_file( '.'.join([ basename, 'plotslides', 'sub' ]) )
  plotcumhist_job.set_sub_file( '.'.join([ basename, 'plotcumhist', 'sub' ]) )

# set memory requirments for memory intensive jobs
sql_replace_job.add_condor_cmd("request_memory", "1000")
sql_fromcache_job.add_condor_cmd("request_memory", "1000")
sql_frominput_job.add_condor_cmd("request_memory", "1000")
sql_extract_job.add_condor_cmd("request_memory", "1000")
injfind_job.add_condor_cmd("request_memory", "1000")
printlc_job.add_condor_cmd("request_memory", "1000")
printsims_job.add_condor_cmd("request_memory", "1000")
printsims_nominifup_job.add_condor_cmd("request_memory", "1000")
minifup_job.add_condor_cmd("request_memory", "2000")
minifup_ts_job.add_condor_cmd("request_memory", "2000")
minifup_cm_job.add_condor_cmd("request_memory", "2000")

##############################################################################
# Cycle over the user tags, creating databases for each 

veto_files = {}
cluster_nodes = {}
injfind_parents = {}
sim_tags = []

for tag in user_tags:

  print "Creating jobs for %s..." % tag
  
  # determine whether or not this was an injection run by checking if there
  # is an injection file for this tag

  file_sieve = '_'.join([ r'INJECTIONS_[0-9]{1,}', tag.split('_CAT_')[0] ]) + r'$'
  inj_file = [ entry for entry in inj_cache if re.match(file_sieve, entry.description) is not None ]
  if len(inj_file) == 0:
    simulation = False
  elif len(inj_file) == 1:
    simulation = True
    inj_file = inj_file[0].url
    sim_tags.append(tag.split('_CAT_')[0])
  else:
    raise ValueError, "More than one injection file found for %s" % tag

  ############################################################################
  # Creating the thinca_user_tag cache file and writing it to disk

  print "\tcreating the thinca_usertag cache file..."

  # sieve thinca_cache for THINCA files with this tag
  file_sieve = '*' + tag
  thinca_usertag_cache = thinca_cache.sieve( description = file_sieve, exact_match = True )

  # get distinct on_instruments in the thinca_usertag cache
  distinct_instrument_sets = set([ entry.observatory for entry in thinca_usertag_cache ])
  # from the distinct_instrument_sets figure out what distinct ifos there are
  distinct_ifo_set = set()
  for on_instruments in distinct_instrument_sets:
    distinct_ifo_set |= lsctables.instrument_set_from_ifos(on_instruments)
  distinct_ifos = ''.join(sorted(distinct_ifo_set))

  # determine location of veto-file relevant to tag
  cat_num = str(get_veto_cat_from_tag( tag ))
  veto_file_path = cp.get('input', 'ihope-segments-directory')
  veto_file_name = ''.join([
    distinct_ifos, '-VETOTIME_CAT_', cat_num, '-', experiment_start, '-', experiment_duration,
    '.xml' ])
  veto_file = '/'.join([ veto_file_path, veto_file_name ])
  if not os.path.exists( veto_file ):
    raise ValueError, "Veto file %s could not be found." % veto_file
  # store the veto file for additional later use
  veto_cat = '_'.join(['CAT', cat_num, 'VETO'])
  veto_files[veto_cat] = veto_file

  # create cache name from what's in the output cache
  cache_type = thinca_usertag_cache[0].description
  thinca_cache_name = '.'.join([ 
    '-'.join([ distinct_ifos, cache_type, experiment_start, experiment_duration ]), 
    'cache' ])
  # write output cache
  fp = open( thinca_cache_name, 'w' )
  thinca_usertag_cache.tofile( fp )
  fp.close()

  ############################################################################
  # Setup a LigolwSqliteNode for putting thinca into a sql db

  print "\tsetting up node to put thinca files into a SQLite database..."
  
  # set node options
  t2sql_node = pipeline.LigolwSqliteNode( sql_replace_job )
  t2sql_node.set_category('ligolw_sqlite')
  t2sql_node.set_input_cache( thinca_cache_name )
  t2sql_node.set_tmp_space( tmp_space )
  # database name has form: 
  # distinct_ifos-USER_TAG_RAW_CBC_RESULTS-cache_start-cache_duration.sqlite
  db_type = tag + '_RAW_CBC_RESULTS'
  raw_result_db = '.'.join([
    '-'.join([ distinct_ifos, db_type, experiment_start, experiment_duration ]),
    'sqlite' ])
  # check to make sure the database doesn't already exist
  if os.path.exists( raw_result_db ):
    print "WARNING: Raw result database %s already exists; " % raw_result_db + \
    "if it isn't moved, it will be overwritten when DAG is submitted."
    
  t2sql_node.set_database( raw_result_db )
  
  dag.add_node( t2sql_node )
  
  ############################################################################
  # Setup a DBSimplifyNode to clean up the output of the t2sql_node 
  
  print "\tsetting up dbsimplify node to clean the database..."
  
  # set node options
  dbsimplify_node = inspiral.DBSimplifyNode( dbsimplify_job )
  dbsimplify_node.set_category('dbsimplify')
  dbsimplify_node.set_tmp_space( tmp_space )
  dbsimplify_node.set_database( raw_result_db )
  
  # set parent node
  dbsimplify_node.add_parent( t2sql_node )
  
  dag.add_node( dbsimplify_node )

  ############################################################################
  # Make the detection statistic nicer and Kipp and Duncan happier
 
  print "\tsetting up repop_coinc node to recalculate detection statistic..."

  # set node options
  last_node = dbsimplify_node
  for repop_coinc_tag, repop_coinc_job in repop_coinc_jobs:
      repop_coinc_node = inspiral.RepopCoincNode( repop_coinc_job )
      repop_coinc_node.set_category('repop_coinc')
      repop_coinc_node.set_tmp_space( tmp_space )
      repop_coinc_node.set_input( raw_result_db )
      repop_coinc_node.set_output( raw_result_db )

      # set parent node
      repop_coinc_node.add_parent( last_node )

      dag.add_node( repop_coinc_node )
      # update last node to be the last repop_coinc node
      last_node = repop_coinc_node

  #############################################################################
  # Setup a ClusterCoincsNode to cluster the output of dbsimplify_node
  
  print "\tsetting up cluster node to cluster coincs in the database..."
  
  # set node options
  # for the first clustering job, we want the input to be the raw_result_db,
  # and in the odd case that someone doesn't want any clustering (bad idea!)
  # the result_db is just the raw_result db
  result_db = raw_result_db
  for cluster_coincs_tag, cluster_coincs_job in cluster_coincs_jobs:
      cluster_coincs_node = inspiral.ClusterCoincsNode( cluster_coincs_job )
      cluster_coincs_node.set_category('cluster_coincs')
      cluster_coincs_node.set_tmp_space( tmp_space )
      cluster_coincs_node.set_input( result_db )
      cluster_coincs_node.add_var_opt('time-column',time_column)
      # output database name has form:
      # distinct_ifos-CBC_TRIGDB_CLUSTERED-USER_TAG-gps_start_time-durations.sqlite
      db_type = tag + '_CLUSTERED_CBC_RESULTS'
      result_db = '.'.join([
        '-'.join([ distinct_ifos, db_type, experiment_start, experiment_duration ]),
        'sqlite' ])
      cluster_coincs_node.set_output(result_db)

      # set parent node; this is the last repop coinc node (if there were any; 
      # if not, the dbsimplify node)
      cluster_coincs_node.add_parent( last_node )

      dag.add_node( cluster_coincs_node )
      # update last node to be the last cluster_coincs node
      last_node = cluster_coincs_node

  # add to list of parents for veto2sql nodes
  cluster_nodes[tag] = last_node
  
  ############################################################################
  # Do additional jobs for injection tags

  if simulation:
    # add dbaddinj node
    print "\tsetting up dbaddinj node to add the injection file..."

    # set node options
    dbaddinj_node = inspiral.DBAddInjNode( dbaddinj_job )
    dbaddinj_node.set_category('dbaddinj')
    dbaddinj_node.set_tmp_space( tmp_space )
    dbaddinj_node.set_database( result_db )
    dbaddinj_node.set_injection_file( inj_file )
    dbaddinj_node.set_inj_tag(sim_tags[-1])
    dbaddinj_node.add_parent( last_node )

    dag.add_node( dbaddinj_node )

    # add sqlite extract node
    print "\tsetting up ligolw_sqlite node to extract the injection database to an xml..."
    
    # set node options
    simxml_node = pipeline.LigolwSqliteNode( sql_extract_job )
    simxml_node.set_category('ligolw_sqlite')
    simxml_node.set_tmp_space( tmp_space )
    simxml_node.set_database( result_db )
    sim_xml = result_db.replace('.sqlite', '.xml')
    simxml_node.set_xml_output( sim_xml )

    simxml_node.add_parent( dbaddinj_node )
    dag.add_node( simxml_node )

    # add to list of injfind parents
    if veto_cat not in injfind_parents:
      injfind_parents[veto_cat] = []
    injfind_parents[veto_cat].append( simxml_node )

  else:
    # just cache the result_db
    non_sim_dbs.append( result_db )
    

##############################################################################
# done cycling over tags: Create injfind job and node

# cache the sim xmls by veto category
print "Creating injfind nodes..."

injfind_nodes = {}
sim_caches = {}
sim_tags = set(sim_tags+['ALLINJ'])
for veto_cat, node_list in injfind_parents.items():
  
  # create a injfind node for each veto_category
  injfind_node = inspiral.InspInjFindNode( injfind_job )
  injfind_node.set_category('injfind')

  # add input files and parents
  for simxml_node in node_list:
    injfind_node.add_file_arg( simxml_node.get_output() )
    injfind_node.add_parent( simxml_node )
  dag.add_node( injfind_node )

  injfind_nodes[veto_cat] = injfind_node

  # Cache the files
  sim_cache = lal.Cache().from_urls( injfind_node.get_input_files() )
  if veto_cat == 'CAT_1_VETO':
    cache_type = 'ALLINJ_CLUSTERED_CBC_RESULTS'
  else:
    cache_type = '_'.join(['ALLINJ', veto_cat, 'CLUSTERED_CBC_RESULTS' ])
  instruments = set()
  for entry in sim_cache:
    instruments |= lsctables.instrument_set_from_ifos(entry.observatory)
  instruments = ''.join(sorted(instruments))
  sim_cache_name = '.'.join([ 
    '-'.join([ instruments, cache_type, experiment_start, experiment_duration ]), 
    'cache' ])
  sim_cache.tofile( open(sim_cache_name, 'w') )
  sim_caches[veto_cat] = sim_cache_name


##############################################################################
# now cycle over non-injection databases, 
# carrying out the rest of the pipeline

# cache the non_sim_dbs
result_dbs_cache = lal.Cache().from_urls( non_sim_dbs )

for result_db in result_dbs_cache:
  
  # get tag and veto_cat
  print "Creating jobs for %s database..." % result_db.description
  tag = result_db.description.replace('_CLUSTERED_CBC_RESULTS', '')
  cat_num = get_veto_cat_from_tag( tag )
  veto_cat = '_'.join([ 'CAT', str(cat_num), 'VETO' ])

  # get all possible instruments_on in this database
  instruments = lsctables.instrument_set_from_ifos(result_db.observatory)
  distinct_instrument_sets = [instruments]
  distinct_instrument_sets.extend( set(sub_combo)
    for nn in range(2, len(instruments))
    for sub_combo in iterutils.choices( list(instruments), nn ) )

  # add the injection xmls to the FULL_DATA databases
  if 'FULL_DATA' in tag and veto_cat in sim_caches:

    # create a sqlite node to add the injetion results
    sim2fulldb_node = pipeline.LigolwSqliteNode( sql_fromcache_job )
    sim2fulldb_node.set_category('ligolw_sqlite')
    sim2fulldb_node.set_input_cache( sim_caches[veto_cat] )
    sim2fulldb_node.set_database( result_db.path )
    sim2fulldb_node.set_tmp_space( tmp_space )
    
    sim2fulldb_node.add_parent( injfind_nodes[veto_cat] )
    sim2fulldb_node.add_parent( cluster_nodes[tag] )
    dag.add_node( sim2fulldb_node )

    # create a dbsimplify node to clean the database
    dbsimplify2_node = inspiral.DBSimplifyNode( dbsimplify_job )
    dbsimplify2_node.set_category('dbsimplify')
    dbsimplify2_node.set_tmp_space( tmp_space )
    dbsimplify2_node.set_database( result_db.path )

    dbsimplify2_node.add_parent( sim2fulldb_node )
    dag.add_node( dbsimplify2_node )

    print "\tsetting up dbinjfind node to add exact/nearby definitions..."

    # set dbinjfind node options
    last_node = dbsimplify2_node
    for dbinjfind_tag, dbinjfind_job in dbinjfind_jobs:
        dbinjfind_node = inspiral.DBInjFindNode( dbinjfind_job )
        dbinjfind_node.set_category('dbinjfind')
        dbinjfind_node.set_tmp_space( tmp_space )
        dbinjfind_node.set_input( result_db.path )
        dbinjfind_node.set_output( result_db.path )

        # set parent node
        dbinjfind_node.add_parent( last_node )

        dag.add_node( dbinjfind_node )
        # update last node to be the last dbinjfind node
        last_node = dbinjfind_node

  else:
    last_node = cluster_nodes[tag]

  ############################################################################
  # Setup a LigolwSqliteNode for putting the veto-segments file into the
  # database 
  
  # write node to add the veto file
  veto2sql_node = pipeline.LigolwSqliteNode( sql_frominput_job )
  veto2sql_node.set_category('ligolw_sqlite')
  veto2sql_node.add_file_arg( veto_files[veto_cat] )
  veto2sql_node.set_tmp_space( tmp_space )
  veto2sql_node.set_database( result_db.path )

  # set parent node
  veto2sql_node.add_parent( last_node )

  dag.add_node( veto2sql_node )
                     
  ############################################################################
  # Compute durations in the database

  print "\tsetting up compute_durations node..."

  # set node options
  comp_durs_node = inspiral.ComputeDurationsNode( comp_durs_job)
  comp_durs_node.set_category('compute_durations')
  comp_durs_node.set_tmp_space( tmp_space )
  comp_durs_node.set_database( result_db.path )

  # set parent node
  comp_durs_node.add_parent( veto2sql_node )

  dag.add_node(comp_durs_node)
  
  ############################################################################
  # MVSC Calculation

  if 'FULL_DATA' in tag and veto_cat in sim_caches:
    print "\tsetting up MVSC dag..."
    mvsc_dag_name = options.config_file.replace('.ini','')+'_mvsc_'+tag+'_n'+cp.get("mvsc_dag","number-of-trees")+'_l'+cp.get("mvsc_dag","leaf-size")+'_s'+cp.get("mvsc_dag","sampled-parameters")+'_c'+cp.get("mvsc_dag","criterion-for-optimization")+'.dag'
    mvsc_dag_generator_job = inspiral.MVSCDagGenerationJob(cp)
    for key,val in cp.items("mvsc_dag"):
      mvsc_dag_generator_job.add_opt(key,val)
    mvsc_dag_generator_job.add_opt("ini-file", options.config_file)
    mvsc_dag_generator_node = inspiral.MVSCDagGenerationNode(mvsc_dag_generator_job)
    mvsc_dag_generator_node.set_user_tag(tag)
    mvsc_dag_generator_node.set_database(result_db.path)
    mvsc_dag_generator_node.add_parent(comp_durs_node)
    dag.add_node(mvsc_dag_generator_node)
    mvsc_dag_job = pipeline.CondorDAGManJob(mvsc_dag_name, os.getcwd(), None)
    mvsc_dag_node = pipeline.CondorDAGManNode(mvsc_dag_job)
    mvsc_dag_node.add_parent(mvsc_dag_generator_node)
    dag.add_node(mvsc_dag_node)

  ############################################################################
  # Compute the uncombined false alarm rates
  
  print "\tsetting up cfar nodes:"
  print "\t\tfor uncombined false alarm rates..."
  
  # set node options: output database is same as input
  ucfar_node = inspiral.CFarNode( ucfar_job )
  ucfar_node.set_category('cfar')
  ucfar_node.set_tmp_space( tmp_space )
  ucfar_node.set_input( result_db.path )
  ucfar_node.set_output( result_db.path )
  
  # set parent node
  ucfar_node.add_parent( comp_durs_node )
  if 'FULL_DATA' in tag and veto_cat in sim_caches:
    ucfar_node.add_parent( mvsc_dag_node )
  
  dag.add_node( ucfar_node )
  
  ############################################################################
  # Compute the combined false alarm rates
  
  print "\t\tfor combined false alarm rates..."
  
  # set node options: output database is same as input
  ccfar_node = inspiral.CFarNode( ccfar_job )
  ccfar_node.set_category('cfar')
  ccfar_node.set_tmp_space( tmp_space )
  ccfar_node.set_input( result_db.path )
  ccfar_node.set_output( result_db.path )
  
  # set parent node
  ccfar_node.add_parent( ucfar_node )
  
  dag.add_node( ccfar_node )

  ############################################################################
  # UPPERLIMIT CALCULATIONS
  ############################################################################

  # Compute the upper limits as a function of mass for these types of mass binnings
  if 'FULL_DATA' in tag:
    #for full data, use FAR of loudest event
    output_cache_name = ''.join(sorted(instruments)) + "-SEARCH_VOLUME_" + tag + ".cache"
    search_volume_node = inspiral.SearchVolumeNode( search_volume_job )
    search_volume_node.set_open_box()
    search_volume_node.add_database( result_db.path )
    search_volume_node.set_user_tag( tag )
    search_volume_node.set_output_cache( output_cache_name )
    search_volume_node.set_veto_segments_name( "VETO_CAT"+str(get_veto_cat_from_tag(tag))+"_CUMULATIVE" )
    search_volume_node.add_parent( ccfar_node )

    search_upper_limit_node = inspiral.SearchUpperLimitNode( search_upper_limit_job )
    search_upper_limit_node.set_open_box()
    search_upper_limit_node.add_input_cache( output_cache_name  )
    search_upper_limit_node.set_user_tag( tag )
    search_upper_limit_node.add_parent( search_volume_node )
    dag.add_node( search_volume_node )
    dag.add_node( search_upper_limit_node )

    # make corresponding jobs for playground (using FAR of expected loudest event)
    playtag = tag.replace("FULL_DATA","PLAYGROUND")
    search_volume_node = inspiral.SearchVolumeNode( search_volume_job )
    output_cache_name = output_cache_name.replace("FULL_DATA","PLAYGROUND")
    search_volume_node.add_database( result_db.path )
    search_volume_node.set_user_tag( playtag )
    search_volume_node.set_output_cache( output_cache_name )
    search_volume_node.set_veto_segments_name( "VETO_CAT"+str(get_veto_cat_from_tag(tag))+"_CUMULATIVE" )
    search_volume_node.add_parent( ccfar_node )

    search_upper_limit_node = inspiral.SearchUpperLimitNode( search_upper_limit_job )
    search_upper_limit_node.add_input_cache( output_cache_name )
    search_upper_limit_node.set_user_tag( playtag )
    search_upper_limit_node.add_parent( search_volume_node )
    dag.add_node( search_volume_node )
    dag.add_node( search_upper_limit_node )


  ############################################################################
  # Summary: Setup PrintLC and MiniFollowup Nodes to generate a summary of 
  # loudest non-simulation events

  print "\tsetting up printlc and minifollowup nodes..."

  # set datatypes to generate files for
  if 'PLAYGROUND' in tag:
    datatypes = ['playground', 'slide']
  else:
    datatypes = ['all_data', 'playground', 'slide']

  for datatype in datatypes:
    print "\t\tfor %s..." % datatype

    # set file naming type
    type_prefix = tag
    type = '_'.join([ type_prefix, 'LOUDEST', datatype.upper(), 'EVENTS_BY', cp.get('printlc', 'ranking-stat').upper()])
    # cycle over all ifos times, creating different tables for each
    for on_instruments in distinct_instrument_sets:
      on_instruments = lsctables.ifos_from_instrument_set(on_instruments)
      # set output and extracted xml file names
      if cp.has_option('printlc', 'output-format'):
        summ_file_type = cp.get('printlc', 'output-format')
      else:
        summ_file_type = 'xml'
      summary_filename = '.'.join([
        '-'.join([ ''.join(on_instruments.split(',')), type + '_SUMMARY', experiment_start, experiment_duration ]),
        summ_file_type ])
      full_filename = '-'.join([ ''.join(on_instruments.split(',')), type, experiment_start, experiment_duration ])
      # set node options
      printlc_node = inspiral.PrintLCNode( printlc_job )
      printlc_node.set_category('printlc')
      printlc_node.set_tmp_space( tmp_space )
      printlc_node.set_input( result_db.path )
      printlc_node.set_output( summary_filename )
      printlc_node.set_extract_to_xml( '.'.join([ full_filename, 'xml' ]) )
      printlc_node.set_extract_to_database( '.'.join([ full_filename, 'sqlite' ]) )
      printlc_node.set_include_only_coincs( '[ALLin' + on_instruments + ']' )
      printlc_node.set_datatype( datatype )
      printlc_node.add_var_opt('time-column',time_column)

      # set parent node
      printlc_node.add_parent( ccfar_node )

      dag.add_node( printlc_node )

      # create the minifollowups nodes (not for CAT_1)
      if veto_cat != 'CAT_1_VETO':
        prefix = '-'.join([ ''.join(on_instruments.split(',')), tag ])
        suffix = re.sub(prefix + '_', '', summary_filename).rstrip('.xml')
  
        if datatype == 'slide':
          minifup_node = inspiral.MiniFollowupsNode( minifup_ts_job )
        else:
          minifup_node = inspiral.MiniFollowupsNode( minifup_job )
        minifup_node.set_category('minifollowups')
        minifup_node.set_cache_file( options.ihope_cache )
        minifup_node.set_cache_string( tag.split('_CAT_')[0] )
        minifup_node.set_prefix( prefix )
        minifup_node.set_suffix( suffix )
        minifup_node.set_input_xml( '.'.join([ full_filename, 'xml' ]) )
        minifup_node.set_input_xml_summary( summary_filename )
        minifup_node.set_output_html_table( re.sub('.xml', '.html', summary_filename ) )
        minifup_node.set_table_name( 'loudest_events' )
  
        minifup_node.add_parent( printlc_node )
        dag.add_node( minifup_node )
      
  ############################################################################
  # Injection summary: Setup printsims, minifollowup, and printmissed nodes

  if 'PLAYGROUND' not in tag:
    print "\tsetting up injection summary nodes..."

    # cycle over all the different types of injections
    for sim_tag in sim_tags:
  
      # 
      # Printsims
      #
      datatypes = ['playground']
      # TODO: if you're wondering why the list only has one entry... this originally cycled over
      # all_data and playground. As the rank column was dropped from minifollowups, and since
      # this would require wip to pick up two different files (playground for closed box, all_data for open)
      # I'm just dropping the all_data for now. If we decide to put it back later, I'll do so; otherwise,
      # if it isn't a useful feature, I'll drop the list and hardcode a datatype below at a future date.
  
      for datatype in datatypes:
        # set file naming type
        type = '_'.join([ sim_tag, veto_cat, 'QUIETEST_FOUND_COMPARED_TO', datatype.upper(), 'BY', cp.get('printsims', 'ranking-stat').upper()])
        # cycle over all ifos times, creating different tables for each
        for on_instruments in distinct_instrument_sets:
          on_instruments = lsctables.ifos_from_instrument_set(on_instruments)
          # set output-format to html for ALLINJ since won't run minifup on these types
          if sim_tag == 'ALLINJ':
            summ_file_type = 'html'
            columns = ','.join([
              'injected_decisive_distance',
              'injected_gps_time',
              'injected_event_time_utc__Px_click_for_daily_ihope_xP_',
              'elogs',
              'injected_mass1',
              'injected_mass2',
              'sim_tag',
              'recovered_combined_far',
              'recovered_snr'])
          elif cp.has_option('printsims', 'output-format'):
            summ_file_type = cp.get('printsims', 'output-format')
          else:
            summ_file_type = 'xml'
          summary_filename = '.'.join([
            '-'.join([ ''.join(on_instruments.split(',')), type + '_SUMMARY', experiment_start, experiment_duration ]),
            summ_file_type ])
          full_filename = '-'.join([ ''.join(on_instruments.split(',')), type, experiment_start, experiment_duration ])
          # set node options
          if sim_tag == 'ALLINJ':
            printsims_node = inspiral.PrintSimsNode( printsims_nominifup_job )
            printsims_node.set_columns( columns )
          else:
            printsims_node = inspiral.PrintSimsNode( printsims_job )
          printsims_node.set_category('printsims')
          printsims_node.set_tmp_space( tmp_space )
          printsims_node.set_input( result_db.path )
          printsims_node.set_output( summary_filename )
          printsims_node.set_extract_to_xml( '.'.join([ full_filename, 'xml' ]) )
          printsims_node.set_include_only_coincs( '[ALLin' + on_instruments + ']' )
          printsims_node.set_comparison_datatype( datatype )
          printsims_node.set_sim_tag( sim_tag )
          printsims_node.set_output_format( summ_file_type )
          printsims_node.add_var_opt('time-column',time_column)
    
          # set parent node
          printsims_node.add_parent( ccfar_node )
    
          dag.add_node( printsims_node )
    
          # create the minifollowups nodes
          # this is not done for ALLINJ because there is no ALLINJ tag in the ihope_cache
          if not (veto_cat == 'CAT_1_VETO' or sim_tag == 'ALLINJ'):
            prefix = '-'.join([ ''.join(on_instruments.split(',')), '_'.join([sim_tag, veto_cat]) ])
            suffix = re.sub(prefix + '_', '', summary_filename).rstrip('.xml')
      
            minifup_node = inspiral.MiniFollowupsNode( minifup_job )
            minifup_node.set_category('minifollowups')
            minifup_node.set_cache_file( options.ihope_cache )
            minifup_node.set_cache_string( sim_tag )
            minifup_node.set_prefix( prefix )
            minifup_node.set_suffix( suffix )
            minifup_node.set_input_xml( '.'.join([ full_filename, 'xml' ]) )
            minifup_node.set_input_xml_summary( summary_filename )
            minifup_node.set_output_html_table( re.sub('.xml', '.html', summary_filename ) )
            minifup_node.set_table_name( 'selected_found_injections' )
      
            minifup_node.add_parent( printsims_node )
            dag.add_node( minifup_node )
      
      #
      #  Printmissed
      #
      type = '_'.join([ sim_tag, veto_cat, 'CLOSEST_MISSED_INJECTIONS' ])
      # cycle over all ifos times, creating different tables for each
      for on_instruments in distinct_instrument_sets:
        on_instruments = lsctables.ifos_from_instrument_set(on_instruments)
        # set output and extracted xml file names
        if sim_tag == 'ALLINJ':
          summ_file_type = 'html'
          columns = ','.join([
            'rank',
            'decisive_distance',
            'gps_time',
            'injection_time_utc__Px_click_for_daily_ihope_xP_',
            'elogs',
            'mchirp',
            'mass1',
            'mass2',
            'eff_dist_h',
            'eff_dist_l',
            'eff_dist_v',
            'sim_tag',
            'mini_followup'
            'mass1',
            'mass2',
            'eff_dist_h',
            'eff_dist_l',
            'eff_dist_v',
            'sim_tag'])
        elif cp.has_option('printmissed', 'output-format'):
          summ_file_type = cp.get('printmissed', 'output-format')
        else:
          summ_file_type = 'xml'
        summary_filename = '.'.join([
          '-'.join([ ''.join(on_instruments.split(',')), type + '_SUMMARY', experiment_start, experiment_duration ]),
          summ_file_type ])

        # set node options
        if sim_tag == 'ALLINJ':
          printmissed_node = inspiral.PrintMissedNode( printmissed_nominifup_job )
          printmissed_node.set_columns( columns )
        else:
          printmissed_node = inspiral.PrintMissedNode( printmissed_job )
        printmissed_node.set_category('printmissed')
        printmissed_node.set_tmp_space( tmp_space )
        printmissed_node.set_input( result_db.path )
        printmissed_node.set_output( summary_filename )
        printmissed_node.set_include_only_coincs( '[ALLin' + on_instruments + ']' )
        printmissed_node.set_sim_tag( sim_tag )
        printmissed_node.set_output_format( summ_file_type )
  
        # set parent node
        printmissed_node.add_parent( ccfar_node )
  
        dag.add_node( printmissed_node )

        # create the minifollowups nodes
        # this is not done for ALLINJ because there is no ALLINJ tag in the ihope_cache
        if not (veto_cat == 'CAT_1_VETO' or sim_tag == 'ALLINJ'):
          prefix = '-'.join([ ''.join(on_instruments.split(',')), '_'.join([sim_tag, veto_cat]) ])
          suffix = re.sub(prefix + '_', '', summary_filename).rstrip('.xml')
    
          minifup_node = inspiral.MiniFollowupsNode( minifup_cm_job )
          minifup_node.set_category('minifollowups')
          minifup_node.set_cache_file( options.ihope_cache )
          minifup_node.set_cache_string( sim_tag )
          minifup_node.set_prefix( prefix )
          minifup_node.set_suffix( suffix )
          minifup_node.set_input_xml( summary_filename )
          minifup_node.set_output_html_table( re.sub('.xml', '.html', summary_filename ) )
          minifup_node.set_table_name( 'close_missed_injections' )
    
          minifup_node.add_parent( printmissed_node )
          dag.add_node( minifup_node )

      # create the plotfm node for each plotfm job
      # only do this once per sim_tag as a plotfm job covers all instrument sets
      for plotfm_tag, plotfm_job in plotfm_jobs:
        plotfm_node = inspiral.PlotFMNode( plotfm_job )
        plotfm_node.set_category('plotfm')
        plotfm_node.set_tmp_space( tmp_space )
        plotfm_node.add_file_arg( result_db.path )
        plotfm_node.set_sim_tag( sim_tag )
        plotfm_node.set_user_tag( '_'.join([ plotfm_tag, sim_tag, veto_cat ]) )
           
        plotfm_node.add_parent( ccfar_node )
        dag.add_node( plotfm_node )



  ############################################################################
  # Plotting: Generate all result plots
  
  print "\tsetting up plotting jobs..."

  # Write plotslides node
  print "\t\tcreating plotslides node..."
  if not options.generate_all_data_plots or 'PLAYGROUND' in tag:
    plotslides_node = inspiral.PlotSlidesNode( plotslides_play_job )
  else:
    plotslides_node = inspiral.PlotSlidesNode( plotslides_job )

  plotslides_node.set_category('plotslides')
  plotslides_node.set_tmp_space( tmp_space )
  plotslides_node.set_input( result_db.path )
  plotslides_node.set_user_tag( tag )

  plotslides_node.add_parent( ccfar_node )
  dag.add_node( plotslides_node )

  # create plotcumhist node
  print "\t\tcreating plotcumhist node..."
  if not options.generate_all_data_plots or 'PLAYGROUND' in tag:
    plotcumhist_node = inspiral.PlotCumhistNode( plotcumhist_play_job )
  else:
    plotcumhist_node = inspiral.PlotCumhistNode( plotcumhist_job )

  plotcumhist_node.set_tmp_space( tmp_space )
  plotcumhist_node.set_input( result_db.path )
  plotcumhist_node.set_user_tag( tag )

  plotcumhist_node.add_parent( ccfar_node )

  dag.add_node( plotcumhist_node )

  # Write plotifar nodes for different datatypes
  print "\t\tcreating plotifar node for datatypes:"
  for datatype in ['all_data', 'playground', 'exclude_play']: 
    # only create nodes for non-playground if options.plot-playground-only not set
    if (not options.generate_all_data_plots or 'PLAYGROUND' in tag)  and datatype != 'playground':
      continue
    print "\t\t\t%s..." % datatype
    plotifar_node = inspiral.PlotIfarNode( plotifar_job )
    plotifar_node.set_category('plotifar')
    plotifar_node.set_tmp_space( tmp_space )
    plotifar_node.set_input( result_db.path )
    plotifar_node.set_datatype( datatype )
    plotifar_node.set_user_tag( tag )
    
    # set parent node
    plotifar_node.add_parent( ccfar_node )

    dag.add_node( plotifar_node )
 
  
##############################################################################
##############################################################################
# Final Step: Write the DAG

print "Writing DAG and sub files..."

# set max-jobs: currently, only minifollowups is set
dag.add_maxjobs_category('minifollowups', 15)

dag.write_sub_files()
dag.write_script()
dag.write_dag()

# write process end time to log file and write log file
process.set_process_end_time(proc_id)
utils.write_filename(logdoc, basename+'.log.xml', xsl_file = "ligolw.xsl")

print "Finished!"
print "Now run:\n\tcondor_submit_dag %s" % os.path.basename(dag.get_dag_file())

sys.exit(0)

