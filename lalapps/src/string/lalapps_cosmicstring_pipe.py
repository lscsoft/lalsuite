#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Standalone ring pipeline driver script



This script produces the condor submit and dag files to run
the standalone ring code on LIGO data
"""


import sys, os
import tempfile
import ConfigParser
from optparse import OptionParser


from glue import iterutils
from glue import pipeline
from glue import segments
from glue import segmentsUtils
from glue.lal import CacheEntry
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import segments as ligolwsegments
from glue import offsetvector
from pylal import ligolw_tisi
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from lalapps import cosmicstring
from lalapps import power


__author__ = 'Xavier Siemens<siemens@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'


#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#


def parse_command_line():
	parser = OptionParser(
		usage = "%prog [options] ...",
		description = "FIXME"
	)
	parser.add_option("-f", "--config-file", metavar = "filename", help = "Use this configuration file (required).")
	parser.add_option("-l", "--log-path", metavar = "path", help = "Make condor put log files in this directory (required).")
	parser.add_option("--background-time-slides", metavar = "filename", action = "append", help = "Set the name of the file from which to obtain the time slide table for use in the background branch of the pipeline (required).  This option can be given multiple times to parallelize the background analysis across time slides.  You will want to make sure the time slide files have distinct vectors to not repeat the same analysis multiple times, and in particular you'll want to make sure only one of them has a zero-lag vector in it.") 
	parser.add_option("--injection-time-slides", metavar = "filename", help = "Set the name of the file from which to obtain the time slide table for use in the injection branch of the pipeline (required).")
	parser.add_option("--segments-file", metavar = "filename", help = "Set the name of the LIGO Light-Weight XML file from which to obtain segment lists (required).  See ligolw_segments and ligolw_segment_query for more information on constructing an XML-format segments file.  See also --segments-name.")
	parser.add_option("--segments-name", metavar = "name", default = "segments", help = "Set the name of the segment lists to retrieve from the segments file (default = \"segments\").  See also --segments-file.")
	parser.add_option("--vetoes-file", metavar = "filename", help = "Set the name of the LIGO Light-Weight XML file from which to obtain veto segment lists (optional).  See ligolw_segments and ligolw_segment_query for more information on constructing an XML-format segments file.  See also --vetos-name.")
	parser.add_option("--vetoes-name", metavar = "name", default = "vetoes", help = "Set the name of the segment lists to retrieve from the veto segments file (default = \"vetoes\").  See also --vetoes-file.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")

	options, filenames = parser.parse_args()

	required_options = ["log_path", "config_file", "background_time_slides", "injection_time_slides", "segments_file"]
	missing_options = [option for option in required_options if getattr(options, option) is None]
	if missing_options:
		raise ValueError, "missing required options %s" % ", ".join(sorted("--%s" % option.replace("_", "-") for option in missing_options))

	if options.vetoes_file is not None:
		options.vetoes_cache = set([CacheEntry(None, None, None, "file://localhost" + os.path.abspath(options.vetoes_file))])
	else:
		options.vetoes_cache = set()

	options.injection_time_slides = [options.injection_time_slides]

	return options, filenames

options, filenames = parse_command_line()


#
# =============================================================================
#
#                                  Initialize
#
# =============================================================================
#


#
# start a log file for this script
#

basename = os.path.splitext(os.path.basename(options.config_file))[0]
log_fh = open(basename + '.pipeline.log', 'w')
# FIXME: the following code uses obsolete CVS ID tags.
# It should be modified to use git version information.
print >>log_fh, "$Id$\nInvoked with arguments:"
for name_value in options.__dict__.items():
	print >>log_fh, "%s %s" % name_value
print >>log_fh

#
# create the config parser object and read in the ini file
#

config_parser = ConfigParser.ConfigParser()
config_parser.read(options.config_file)

#
# initialize lalapps.power and lalapps.cosmicstring modules
#

power.init_job_types(config_parser, job_types = ("datafind", "binj", "lladd", "binjfind", "burca", "sqlite"))
cosmicstring.init_job_types(config_parser, job_types = ("string", "meas_likelihood", "calc_likelihood", "runsqlite"))

#
# make directories to store the cache files, job logs, and trigger output
#

def make_dag_directories(top_level_directory, config_parser):
	cwd = os.getcwd()
	power.make_dir_if_not_exists(top_level_directory)
	os.chdir(top_level_directory)
	power.make_dag_directories(config_parser)
	# FIXME:  move this into make_dag_directories().  requires update
	# of excess power and gstlal dags
	power.make_dir_if_not_exists(power.get_triggers_dir(config_parser))
	os.chdir(cwd)

power.make_dag_directories(config_parser)
injection_folders = []
for i in range(config_parser.getint('pipeline', 'injection-runs')):
	injection_folders.append(os.path.abspath("injections%d" % i))
	make_dag_directories(injection_folders[-1], config_parser)
noninjection_folders = []
noninjection_folders.append(os.path.abspath("noninjections"))
make_dag_directories(noninjection_folders[-1], config_parser)

#
# create a log file that the Condor jobs will write to
#

logfile = tempfile.mkstemp(prefix = basename, suffix = '.log', dir = options.log_path)[1]

#
# create the DAG writing the log to the specified directory
#

dag = pipeline.CondorDAG(logfile)
dag.set_dag_file(basename)
clipsegments_sql_filename = os.path.abspath("clipsegments.sql")

#
# get chunk lengths from the values in the ini file
#

short_segment_duration = config_parser.getint('lalapps_StringSearch', 'short-segment-duration')
pad = config_parser.getint('lalapps_StringSearch', 'pad')
min_segment_length = config_parser.getint('pipeline', 'segment-length') # not including pad at each end
trig_overlap = config_parser.getint('pipeline', 'trig_overlap')
overlap = short_segment_duration / 2 + 2 * pad	# FIXME:  correct?

#
# get the instruments and raw segments
#

instruments = lsctables.instrument_set_from_ifos(config_parser.get('pipeline','ifos'))
seglists = ligolwsegments.segmenttable_get_by_name(utils.load_filename(options.segments_file, gz = (options.segments_file or "stdin").endswith(".gz"), verbose = options.verbose), options.segments_name).coalesce()
# remove extra instruments
for instrument in set(seglists) - instruments:
	if options.verbose:
		print >>sys.stderr, "warning: ignoring segments for '%s' found in '%s'" % (instrument, options.segments_file)
	del seglists[instrument]
# check for missing instruments
if not instruments.issubset(set(seglists)):
	raise ValueError, "segment lists retrieved from '%s' missing segments for instruments %s" % (options.segments_file, ", ".join(instruments - set(seglists)))
# now rely on seglists' keys to provide the instruments
del instruments

#
# Using time slide information, construct segment lists describing times
# requiring trigger construction.
#

if options.verbose:
	print >>sys.stderr, "Computing segments for which lalapps_StringSearch jobs are required ..."

background_time_slides = {}
background_seglists = segments.segmentlistdict()
for filename in options.background_time_slides:
	cache_entry = CacheEntry(None, None, None, "file://localhost" + os.path.abspath(filename))
	background_time_slides[cache_entry] = ligolw_tisi.load_time_slides(filename, verbose = options.verbose, gz = filename.endswith(".gz")).values()
	for i in range(len(background_time_slides[cache_entry])):
		background_time_slides[cache_entry][i] = offsetvector.offsetvector((instrument, LIGOTimeGPS(offset)) for instrument, offset in background_time_slides[cache_entry][i].items())
	background_seglists |= cosmicstring.compute_segment_lists(seglists, ligolw_tisi.time_slide_component_vectors(background_time_slides[cache_entry], 2), min_segment_length, pad)

injection_time_slides = {}
injection_seglists = segments.segmentlistdict()
for filename in options.injection_time_slides:
	cache_entry = CacheEntry(None, None, None, "file://localhost" + os.path.abspath(filename))
	injection_time_slides[cache_entry] = ligolw_tisi.load_time_slides(filename, verbose = options.verbose, gz = filename.endswith(".gz")).values()
	for i in range(len(injection_time_slides[cache_entry])):
		injection_time_slides[cache_entry][i] = offsetvector.offsetvector((instrument, LIGOTimeGPS(offset)) for instrument, offset in injection_time_slides[cache_entry][i].items())
	injection_seglists |= cosmicstring.compute_segment_lists(seglists, ligolw_tisi.time_slide_component_vectors(injection_time_slides[cache_entry], 2), min_segment_length, pad)


#
# check for duplicate offset vectors
#

def check_for_reused_offsetvectors(background_time_slides, injection_time_slides):
	# make sure none of the injection run offset vectors is also used
	# for a non-injection run

	for background_cache_entry, background_offsetvector in [(cache_entry, offsetvector) for cache_entry, offsetvectors in background_time_slides.items() for offsetvector in offsetvectors]:
		for injection_cache_entry, injection_offsetvector in [(cache_entry, offsetvector) for cache_entry, offsetvectors in injection_time_slides.items() for offsetvector in offsetvectors]:
			if not cmp(background_offsetvector, injection_offsetvector):
				raise ValueError, "injections offset vector %s from %s is the same as non-injections offset vector %s from %s.  to avoid a self-selection bias, injections must not be performed at the same relative time shifts as a non-injection run" % (str(injection_offsetvector), injection_cache_entry.url, str(background_offsetvector), background_cache_entry.url)

check_for_reused_offsetvectors(background_time_slides, injection_time_slides)


#
# =============================================================================
#
#                                  Build DAG
#
# =============================================================================
#


#
# datafinds
#


datafinds = power.make_datafind_stage(dag, injection_seglists | background_seglists, verbose = options.verbose)


#
# trigger production, coincidence, injection identification, likelihood
# ratio probability density data collection
#


def make_coinc_branch(dag, datafinds, seglists, time_slides, min_segment_length, pad, overlap, short_segment_duration, tag, vetoes_cache = set(), do_injections = False, injections_offset = 0.0, verbose = False):
	#
	# injection job
	#

	binjnodes = set()
	if do_injections:
		# don't know what to do with more than one list of offset
		# vectors
		assert len(time_slides) == 1

		# get the largest injection offset's magnitude
		maxoffset = max(abs(offset) for offsetvectorlist in time_slides.values() for offsetvector in offsetvectorlist for offset in offsetvector.values())

		# to save disk space and speed the dag along we don't
		# generate a single injection list for the entire analysis
		# run, instead a separate list is constructed for each
		# block of data to be analyzed.  we need to be careful that
		# two nearby injection lists don't contain injections for
		# the same time, so we protract the segments by the time
		# step and coalesce so that only gaps between segments
		# larger than twice the time step result in separate files
		# being generated.  we could allow smaller gaps to survive,
		# but this way we don't have to worry about it.

		# injections_offset is a number between 0 and 1 in units of
		# the period between injections

		for seg in seglists.union(seglists).protract(power.binjjob.time_step + maxoffset).coalesce().contract(power.binjjob.time_step + maxoffset):
			binjnodes |= power.make_binj_fragment(dag, seg.protract(maxoffset), time_slides.keys()[0], tag, offset = injections_offset)

		# artificial parent-child relationship to induce dagman to
		# submit binj jobs as the corresponding datafinds complete
		# instead of submiting all of one kind before any of the next.
		# makes dag run faster because it allows string search jobs to
		# start moving onto the cluster without waiting for all the
		# datafinds and/or all the binjs to complete

		for datafindnode in datafinds:
			seg = segments.segment(datafindnode.get_start(), datafindnode.get_end())
			for binjnode in binjnodes:
				if seg.intersects(power.cache_span(binjnode.get_output_cache())):
					binjnode.add_parent(datafindnode)

	#
	# trigger generator jobs
	#

	# set max job length to ~3600 s (will be clipped to an allowed
	# size)
	trigger_nodes = cosmicstring.make_single_instrument_stage(dag, datafinds, seglists, tag, min_segment_length, pad, overlap, short_segment_duration, max_job_length = 3600, binjnodes = binjnodes, verbose = verbose)

	#
	# coincidence analysis
	#

	coinc_nodes = []
	for n, (time_slides_cache_entry, these_time_slides) in enumerate(time_slides.items()):
		if verbose:
			print >>sys.stderr, "%s %d/%d (%s):" % (tag, n + 1, len(time_slides), time_slides_cache_entry.path)
		coinc_nodes.append(set())

		#
		# ligolw_cafe & ligolw_add
		#

		tisi_cache = set([time_slides_cache_entry])
		lladd_nodes = set()
		for seg, parents, cache, clipseg in power.group_coinc_parents(trigger_nodes, these_time_slides, extentlimit = 50000000.0 / (len(these_time_slides) or 1), verbose = verbose):
			binj_cache = set(cache_entry for node in binjnodes for cache_entry in node.get_output_cache() if cache_entry.segment.intersects(seg))
			# otherwise too many copies of the offset vector
			# will be fed into burca
			assert len(binj_cache) < 2
			if do_injections:
				# lalapps_binj has already copied the time
				# slide document into its own output
				extra_input_cache = vetoes_cache
			else:
				# ligolw_add needs to copy the time slide
				# document into its output
				extra_input_cache = tisi_cache | vetoes_cache
			these_lladd_nodes = power.make_lladd_fragment(dag, parents | binjnodes, "%s_%d" % (tag, n), segment = seg, input_cache = cache | binj_cache, extra_input_cache = extra_input_cache, remove_input = do_injections and clipseg is not None, preserve_cache = binj_cache | tisi_cache | vetoes_cache)
			if clipseg is not None:
				#
				# this is a fragment of a too-large burca
				# job, construct it specially and add the
				# command-line option needed to clip the
				# output
				#

				assert len(these_lladd_nodes) == 1
				coinc_nodes[-1] |= power.make_burca_fragment(dag, these_lladd_nodes, "%s_%d" % (tag, n), coincidence_segments = segments.segmentlist([clipseg]), verbose = verbose)

			else:
				#
				# this is not a fragment of a too-large
				# burca job, add it to the pool of files to
				# be processed by the burcas that don't
				# require special clipping command line
				# options
				#

				lladd_nodes |= these_lladd_nodes

		#
		# ligolw_burca pool.  these are the burca jobs that don't
		# require special clipping command line options, and so can
		# bulk-process many files with each job
		#

		if verbose:
			print >>sys.stderr, "building burca jobs ..."
		coinc_nodes[-1] |= power.make_burca_fragment(dag, lladd_nodes, "%s_%d" % (tag, n), verbose = verbose)
		if verbose:
			print >>sys.stderr, "done %s %d/%d" % (tag, n + 1, len(time_slides))

	#
	# ligolw_binjfind
	#

	if do_injections:
		if verbose:
			print >>sys.stderr, "building binjfind jobs ..."
		coinc_nodes = [power.make_binjfind_fragment(dag, these_coinc_nodes, "%s_%d" % (tag, n), verbose = verbose) for n, these_coinc_nodes in enumerate(coinc_nodes)]

	#
	# ligolw_sqlite and lalapps_run_sqlite
	#

	if verbose:
		print >>sys.stderr, "building sqlite jobs ..."
	coinc_nodes = [power.make_sqlite_fragment(dag, these_coinc_nodes, "%s_%d" % (tag, n), verbose = verbose) for n, these_coinc_nodes in enumerate(coinc_nodes)]
	coinc_nodes = [cosmicstring.make_run_sqlite_fragment(dag, these_coinc_nodes, "%s_%d" % (tag, n), clipsegments_sql_filename) for n, these_coinc_nodes in enumerate(coinc_nodes)]

	#
	# lalapps_string_meas_likelihood
	#

	if verbose:
		print >>sys.stderr, "building lalapps_string_meas_likelihood jobs ..."
	likelihood_nodes = [cosmicstring.make_meas_likelihood_fragment(dag, these_coinc_nodes, "%s_%d" % (tag, n)) for n, these_coinc_nodes in enumerate(coinc_nodes)]

	#
	# write output cache
	#

	if verbose:
		print >>sys.stderr, "writing output cache ..."
	for n, (these_coinc_nodes, these_likelihood_nodes) in enumerate(zip(coinc_nodes, likelihood_nodes)):
		power.write_output_cache(these_coinc_nodes | these_likelihood_nodes, "%s_%s_output.cache" % (os.path.splitext(dag.get_dag_file())[0], "%s_%d" % (tag, n)))

	#
	# done
	#

	return coinc_nodes, likelihood_nodes


user_tag = config_parser.get('pipeline', 'user_tag')


injection_coinc_nodes = []
injection_likelihood_nodes = []
for i, folder in enumerate(injection_folders):
	cwd = os.getcwd()
	os.chdir(folder)
	# note:  unpacking enforces rule that each injection run uses a
	# single time-slides file
	[these_coinc_nodes], [these_likelihood_nodes] = make_coinc_branch(dag, datafinds, injection_seglists, injection_time_slides, min_segment_length, pad, overlap, short_segment_duration, "%s_INJ_%d" % (user_tag, i), vetoes_cache = options.vetoes_cache, do_injections = True, injections_offset = float(i) / len(injection_folders), verbose = options.verbose)
	os.chdir(cwd)
	injection_coinc_nodes.append(these_coinc_nodes)
	injection_likelihood_nodes.append(these_likelihood_nodes)


for i, folder in enumerate(noninjection_folders):
	assert i == 0	# loop only works for one iteration
	cwd = os.getcwd()
	os.chdir(folder)
	background_coinc_nodes, background_likelihood_nodes = make_coinc_branch(dag, datafinds, background_seglists, background_time_slides, min_segment_length, pad, overlap, short_segment_duration, "%s_%d" % (user_tag, i), vetoes_cache = options.vetoes_cache, do_injections = False, verbose = options.verbose)
	os.chdir(cwd)


def flatten_node_groups(node_groups):
	return reduce(lambda a, b: a | b, node_groups, set())


all_background_likelihood_nodes = flatten_node_groups(background_likelihood_nodes)
all_injection_likelihood_nodes = flatten_node_groups(injection_likelihood_nodes)


#
# round-robin parameter distribution density data for likelihood ratio
# computation.
#


if options.verbose:
	print >>sys.stderr, "building lalapps_string_calc_likelihood jobs ..."

def round_robin_and_flatten(injection_coinc_node_groups, injection_likelihood_node_groups):
	# see the documentation for glue.iterutils.choices() for an
	# explanation of the procedure used here to round-robin the node
	# lists
	A = list(iterutils.choices(injection_coinc_node_groups, 1))
	B = list(iterutils.choices(injection_likelihood_node_groups, len(injection_likelihood_node_groups) - 1))
	B.reverse()
	A = [flatten_node_groups(node_groups) for node_groups in A]
	B = [flatten_node_groups(node_groups) for node_groups in B]
	return zip(A, B)

coinc_nodes = set()
for n, (these_inj_coinc_nodes, these_inj_likelihood_nodes) in enumerate(round_robin_and_flatten(injection_coinc_nodes, injection_likelihood_nodes)):
	these_inj_coinc_nodes = cosmicstring.make_calc_likelihood_fragment(dag, these_inj_coinc_nodes, these_inj_likelihood_nodes | all_background_likelihood_nodes, "%s_INJ_%d" % (user_tag, n), verbose = options.verbose)
	# to prevent race conditions, none of the calc_likelihood jobs
	# should start before all the meas_likelihood jobs have completed,
	# even those that produce data not required by the calc_likelihood
	# jobs because the calc_likelihood jobs might overwrite the
	# database files just as the meas_likelihood jobs are opening them
	# for reading.  this is only a problem for injection
	# calc_likelihood jobs because of the round-robining, the
	# non-injection ones already have all meas_likelihood jobs as
	# parents.
	for extra_parent in all_injection_likelihood_nodes - these_inj_likelihood_nodes:
		for coinc_node in these_inj_coinc_nodes:
			coinc_node.add_parent(extra_parent)
	coinc_nodes |= these_inj_coinc_nodes
coinc_nodes |= cosmicstring.make_calc_likelihood_fragment(dag, flatten_node_groups(background_coinc_nodes), all_injection_likelihood_nodes | all_background_likelihood_nodes, user_tag, files_per_calc_likelihood = 1, verbose = options.verbose)



#
# =============================================================================
#
#                                     Done
#
# =============================================================================
#

#
# write DAG
#

if options.verbose:
	print >>sys.stderr, "writing dag ..."
dag.write_sub_files()
dag.write_dag()
cosmicstring.write_clip_segment_sql_file(clipsegments_sql_filename)

#
# write a message telling the user that the DAG has been written
#

print """Created a DAG file which can be submitted by executing

$ condor_submit_dag %s

from a condor submit machine (e.g. hydra.phys.uwm.edu)

Do not forget to initialize your grid proxy certificate on the condor
submit machine by running the commands

$ unset X509_USER_PROXY
$ grid-proxy-init -hours $((24*7)):00

""" % dag.get_dag_file()
