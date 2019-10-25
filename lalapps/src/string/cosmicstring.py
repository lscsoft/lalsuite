#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Classes needed for the cosmic string analysis pipeline.
"""


from __future__ import print_function


import math
import os
import sys


from ligo import segments
from glue import iterutils
from glue import pipeline
from lal import LIGOTimeGPS
from lal.utils import CacheEntry
from lalburst import cafe
from lalapps import power


__author__ = 'Xavier Siemens<siemens@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'


#
# =============================================================================
#
#                                Configuration
#
# =============================================================================
#


def get_files_per_meas_likelihood(config_parser):
	return config_parser.getint("pipeline", "files_per_meas_likelihood")


def get_files_per_calc_likelihood(config_parser):
	return config_parser.getint("pipeline", "files_per_calc_likelihood")


def get_files_per_run_sqlite(config_parser):
	return config_parser.getint("pipeline", "files_per_run_sqlite")


#
# =============================================================================
#
#                            DAG Node and Job Class
#
# =============================================================================
#


class MeasLikelihoodJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", power.get_executable(config_parser, "lalapps_string_meas_likelihood"))
		self.set_sub_file("lalapps_string_meas_likelihood.sub")
		self.set_stdout_file(os.path.join(power.get_out_dir(config_parser), "lalapps_string_meas_likelihood-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(power.get_out_dir(config_parser), "lalapps_string_meas_likelihood-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", power.get_accounting_group(config_parser))
		self.add_ini_opts(config_parser, "lalapps_string_meas_likelihood")

		self.cache_dir = power.get_cache_dir(config_parser)
		self.output_dir = "."
		self.files_per_meas_likelihood = get_files_per_meas_likelihood(config_parser)
		if self.files_per_meas_likelihood < 1:
			raise ValueError("files_per_meas_likelihood < 1")


class MeasLikelihoodNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = []

		self._CondorDAGNode__macros["initialdir"] = os.getcwd()
		self.cache_dir = os.path.join(os.getcwd(), self.job().cache_dir)
		self.output_dir = os.path.join(os.getcwd(), self.job().output_dir)

	def set_name(self, *args):
		pipeline.CondorDAGNode.set_name(self, *args)
		self.cache_name = os.path.join(self.cache_dir, "%s.cache" % self.get_name())
		self.add_var_opt("input-cache", self.cache_name)

	def add_input_cache(self, cache):
		if self.output_cache:
			raise AttributeError("cannot change attributes after computing output cache")
		self.input_cache.extend(cache)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def set_output(self, description):
		if self.output_cache:
			raise AttributeError("cannot change attributes after computing output cache")
		cache_entry = power.make_cache_entry(self.input_cache, description, "")
		filename = os.path.join(self.output_dir, "%s-STRING_LIKELIHOOD_%s-%d-%d.xml.gz" % (cache_entry.observatory, cache_entry.description, int(cache_entry.segment[0]), int(abs(cache_entry.segment))))
		cache_entry.url = "file://localhost" + os.path.abspath(filename)
		self.add_var_opt("output", filename)
		del self.output_cache[:]
		self.output_cache.append(cache_entry)
		return filename

	def get_input_cache(self):
		return  self.input_cache

	def get_output_cache(self):
		if not self.output_cache:
			raise AttributeError("must call set_output(description) first")
		return self.output_cache

	def write_input_files(self, *args):
		f = file(self.cache_name, "w")
		for c in self.input_cache:
			print(str(c), file=f)
		pipeline.CondorDAGNode.write_input_files(self, *args)

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class CalcLikelihoodJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", power.get_executable(config_parser, "lalapps_string_calc_likelihood"))
		self.set_sub_file("lalapps_string_calc_likelihood.sub")
		self.set_stdout_file(os.path.join(power.get_out_dir(config_parser), "lalapps_string_calc_likelihood-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(power.get_out_dir(config_parser), "lalapps_string_calc_likelihood-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", power.get_accounting_group(config_parser))
		self.add_ini_opts(config_parser, "lalapps_string_calc_likelihood")
		self.cache_dir = power.get_cache_dir(config_parser)
		self.files_per_calc_likelihood = get_files_per_calc_likelihood(config_parser)
		if self.files_per_calc_likelihood < 1:
			raise ValueError("files_per_calc_likelihood < 1")


class CalcLikelihoodNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.likelihood_cache = []
		self.output_cache = self.input_cache
		self._CondorDAGNode__macros["initialdir"] = os.getcwd()
		self.cache_dir = os.path.join(os.getcwd(), self.job().cache_dir)

	def set_name(self, *args):
		pipeline.CondorDAGNode.set_name(self, *args)
		self.cache_name = os.path.join(self.cache_dir, "%s.cache" % self.get_name())
		self.add_var_opt("input-cache", self.cache_name)
		self.likelihood_cache_name = os.path.join(self.cache_dir, "%s_likelihood.cache" % self.get_name())
		self.add_var_opt("likelihood-cache", self.likelihood_cache_name)

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)
		for c in cache:
			self.add_output_file(c.path)

	def add_likelihood_cache(self, cache):
		self.likelihood_cache.extend(cache)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def get_input_cache(self):
		return  self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def get_likelihood_cache(self):
		return self.likelihood_cache

	def write_input_files(self, *args):
		f = file(self.cache_name, "w")
		for c in self.input_cache:
			print(str(c), file=f)
		f = file(self.likelihood_cache_name, "w")
		for c in self.likelihood_cache:
			print(str(c), file=f)
		pipeline.CondorDAGNode.write_input_files(self, *args)

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class StringJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_StringSearch job used by the string pipeline. The static options
  are read from the section in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specified in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,config_parser):
    """
    config_parser = ConfigParser object from which options are read.
    """
    pipeline.CondorDAGJob.__init__(self, power.get_universe(config_parser), power.get_executable(config_parser, "lalapps_StringSearch"))
    pipeline.AnalysisJob.__init__(self, config_parser)
    self.add_ini_opts(config_parser, "lalapps_StringSearch")
    self.set_stdout_file(os.path.join(power.get_out_dir(config_parser), "lalapps_StringSearch-$(cluster)-$(process).out"))
    self.set_stderr_file(os.path.join(power.get_out_dir(config_parser), "lalapps_StringSearch-$(cluster)-$(process).err"))
    self.add_condor_cmd("getenv", "True")
    self.add_condor_cmd("accounting_group", power.get_accounting_group(config_parser))
    self.set_sub_file("lalapps_StringSearch.sub")
    #self.add_condor_cmd("Requirements", "Memory > 1100")

    self.output_dir = power.get_triggers_dir(config_parser)


class StringNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  A RingNode runs an instance of the ring code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_StringSearch.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('pipeline','user_tag')
    self.output_cache = []
    self._CondorDAGNode__macros["initialdir"] = os.getcwd()
    self.output_dir = os.path.join(os.getcwd(), self.job().output_dir)

  def set_ifo(self, instrument):
    """
    Load additional options from the per-instrument section in
    the config file.
    """
    if self.output_cache:
      raise AttributeError("cannot change attributes after computing output cache")
    pipeline.AnalysisNode.set_ifo(self, instrument)
    for optvalue in self.job()._AnalysisJob__cp.items("lalapps_StringSearch_%s" % instrument):
      self.add_var_arg("--%s %s" % optvalue)

  def set_user_tag(self, tag):
    if self.output_cache:
      raise AttributeError("cannot change attributes after computing output cache")
    self.__usertag = tag
    self.add_var_opt("user-tag", self.__usertag)

  def get_user_tag(self):
    return self.__usertag

  def get_output_cache(self):
    """
    Returns a LAL cache of the output file name.  Calling this
    method also induces the output name to get set, so it must
    be at least once.
    """
    if not self.output_cache:
      self.output_cache = [CacheEntry(self.get_ifo(), self.__usertag, segments.segment(LIGOTimeGPS(self.get_start()), LIGOTimeGPS(self.get_end())), "file://localhost" + os.path.abspath(self.get_output()))]
    return self.output_cache

  def get_output_files(self):
    raise NotImplementedError

  def get_output(self):
    """
    Returns the file name of output from the ring code. This must be kept
    synchronized with the name of the output file in ring.c.
    """
    if self._AnalysisNode__output is None:
      if None in (self.get_start(), self.get_end(), self.get_ifo(), self.__usertag):
        raise ValueError("start time, end time, ifo, or user tag has not been set")
      seg = segments.segment(LIGOTimeGPS(self.get_start()), LIGOTimeGPS(self.get_end()))
      self.set_output(os.path.join(self.output_dir, "%s-STRINGSEARCH_%s-%d-%d.xml.gz" % (self.get_ifo(), self.__usertag, int(self.get_start()), int(self.get_end()) - int(self.get_start()))))

    return self._AnalysisNode__output

  def set_injection_file(self, file):
    """
    Set the name of the XML file from which to read a list of
    software injections.
    """
    self.add_var_opt("injection-file", file)
    self.add_input_file(file)


class RunSqliteJob(pipeline.CondorDAGJob):
	"""
	A lalapps_run_sqlite job used by the gstlal pipeline. The static
	options are read from the [lalapps_run_sqlite] section in the ini
	file.  The stdout and stderr from the job are directed to the logs
	directory.  The job runs in the universe specified in the ini file.
	The path to the executable is determined from the ini file.
	"""
	def __init__(self, config_parser):
		"""
		config_parser = ConfigParser object
		"""
		pipeline.CondorDAGJob.__init__(self, "vanilla", power.get_executable(config_parser, "lalapps_run_sqlite"))
		self.add_ini_opts(config_parser, "lalapps_run_sqlite")
		self.set_stdout_file(os.path.join(power.get_out_dir(config_parser), "lalapps_run_sqlite-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(power.get_out_dir(config_parser), "lalapps_run_sqlite-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", power.get_accounting_group(config_parser))
		self.set_sub_file("lalapps_run_sqlite.sub")
		self.files_per_run_sqlite = get_files_per_run_sqlite(config_parser)
		if self.files_per_run_sqlite < 1:
			raise ValueError("files_per_run_sqlite < 1")


class RunSqliteNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = self.input_cache
		self._CondorDAGNode__macros["initialdir"] = os.getcwd()

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)
		for c in cache:
			filename = c.path
			pipeline.CondorDAGNode.add_file_arg(self, filename)
			self.add_output_file(filename)

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def set_sql_file(self, filename):
		self.add_var_opt("sql-file", filename)


#
# =============================================================================
#
#                                 Segmentation
#
# =============================================================================
#


def clip_segment_length(segment_length, pad, short_segment_duration):
  # clip segment to the length required by lalapps_StringSearch.  if
  #
  #   duration = segment length - padding
  #
  # then
  #
  #   duration / short_segment_duration - 0.5
  #
  # must be an odd integer, therefore
  #
  #   2 * duration + short_segment_duration
  #
  # must be divisble by (4 * short_segment_duration)
  assert segment_length >= 2 * pad
  duration = segment_length - 2 * pad
  extra = (2 * duration + short_segment_duration) % (4 * short_segment_duration)
  extra /= 2

  # clip
  segment_length -= extra

  # done. negative return value not allowed
  assert segment_length >= 0
  return segment_length


def segment_ok(seg, min_segment_length, pad):
  """
  Return True if the segment seg is long enough to be analyzed by
  lalapps_StringSearch.
  """
  return float(abs(seg)) - 2 * pad >= min_segment_length


def remove_too_short_segments(seglists, min_segment_length, pad):
  """
  Remove segments from the segmentlistdict seglists that are too short to
  analyze.

  CAUTION:  this function modifies seglists in place.
  """
  for seglist in seglists.values():
    iterutils.inplace_filter(lambda seg: segment_ok(seg, min_segment_length, pad), seglist)


def compute_segment_lists(seglists, offset_vectors, min_segment_length, pad):
  # don't modify original
  seglists = seglists.copy()

  # ignore offset vectors referencing instruments we don't have
  offset_vectors = [offset_vector for offset_vector in offset_vectors if set(offset_vector.keys()).issubset(set(seglists.keys()))]

  # cull too-short single-instrument segments from the input
  # segmentlist dictionary;  this can significantly increase
  # the speed of the get_coincident_segmentlistdict()
  # function when the input segmentlists have had many data
  # quality holes poked out of them
  remove_too_short_segments(seglists, min_segment_length, pad)

  # extract the segments that are coincident under the time
  # slides
  new = cafe.get_coincident_segmentlistdict(seglists, offset_vectors)

  # round to integer boundaries because lalapps_StringSearch can't accept
  # non-integer start/stop times
  # FIXME:  fix that in lalapps_StringSearch
  for seglist in new.values():
    for i in range(len(seglist)):
      seglist[i] = segments.segment(int(math.floor(seglist[i][0])), int(math.ceil(seglist[i][1])))
  # intersect with original segments to ensure we haven't expanded beyond
  # original bounds
  new &= seglists

  # again remove too-short segments
  remove_too_short_segments(new, min_segment_length, pad)

  # done
  return new


#
# =============================================================================
#
#                                DAG Job Types
#
# =============================================================================
#


stringjob = None
meas_likelihoodjob = None
calc_likelihoodjob = None
runsqlitejob = None


def init_job_types(config_parser, job_types = ("string", "meas_likelihoodjob", "calc_likelihood", "runsqlite")):
  """
  Construct definitions of the submit files.
  """
  global stringjob, meas_likelihoodjob, calc_likelihoodjob, runsqlitejob

  # lalapps_StringSearch
  if "string" in job_types:
    stringjob = StringJob(config_parser)

  # lalapps_string_meas_likelihood
  if "meas_likelihood" in job_types:
    meas_likelihoodjob = MeasLikelihoodJob(config_parser)

  # lalapps_string_calc_likelihood
  if "calc_likelihood" in job_types:
    calc_likelihoodjob = CalcLikelihoodJob(config_parser)

  # lalapps_run_sqlite
  if "runsqlite" in job_types:
    runsqlitejob = RunSqliteJob(config_parser)


#
# =============================================================================
#
#                          lalapps_StringSearch Jobs
#
# =============================================================================
#


#
# one job
#


def make_string_fragment(dag, parents, instrument, seg, tag, framecache, injargs = {}):
	node = StringNode(stringjob)
	node.set_name("lalapps_StringSearch_%s_%s_%d_%d" % (tag, instrument, int(seg[0]), int(abs(seg))))
	map(node.add_parent, parents)
	# FIXME:  StringNode should not be subclassed from AnalysisNode,
	# because that class is too hard-coded.  For example, there is no
	# way to switch to analysing gaussian noise except to comment out
	# this line in the code.
	node.set_cache(framecache)
	node.set_ifo(instrument)
	node.set_start(seg[0])
	node.set_end(seg[1])
	node.set_user_tag(tag)
	for arg, value in injargs.items():
		# this is a hack, but I can't be bothered
		node.add_var_arg("--%s %s" % (arg, value))
	dag.add_node(node)
	return set([node])


#
# one segment
#


def split_segment(seg, min_segment_length, pad, overlap, short_segment_duration, max_job_length):
	# avoid infinite loop
	if min_segment_length + 2 * pad <= overlap:
		raise ValueError("infinite loop: min_segment_length + 2 * pad must be > overlap")

	# clip max_job_length down to an allowed size
	max_job_length = clip_segment_length(max_job_length, pad, short_segment_duration)

	seglist = segments.segmentlist()
	while abs(seg) >= min_segment_length + 2 * pad:
		# try to use max_job_length each time
		if abs(seg) >= max_job_length:
			seglist.append(segments.segment(seg[0], seg[0] + max_job_length))
		else:
			seglist.append(segments.segment(seg[0], seg[0] + clip_segment_length(abs(seg), pad, short_segment_duration)))
		assert abs(seglist[-1]) != 0	# safety-check for no-op
		# bounds must be integers
		if abs((int(seglist[-1][0]) - seglist[-1][0]) / seglist[-1][0]) > 1e-14 or abs((int(seglist[-1][1]) - seglist[-1][1]) / seglist[-1][1]) > 1e-14:
			raise ValueError("segment %s does not have integer boundaries" % str(seglist[-1]))
		# advance segment
		seg = segments.segment(seglist[-1][1] - overlap, seg[1])
	if not seglist:
		raise ValueError("unable to use segment %s" % str(seg))
	return seglist


def make_string_segment_fragment(dag, datafindnodes, instrument, seg, tag, min_segment_length, pad, overlap, short_segment_duration, max_job_length, binjnodes = set(), verbose = False):
	"""
	Construct a DAG fragment for an entire segment, splitting the
	segment into multiple trigger generator jobs.
	"""
	# figure out which binj nodes, if any, produce output for this job
	binjnodes = set(node for node in binjnodes if power.cache_span(node.get_output_cache()).intersects(seg))

	# only one frame cache file can be provided as input, and only one
	# injection description file can be provided as input.
	# the unpacking indirectly tests that the file count is correct
	[framecache] = [node.get_output() for node in datafindnodes]
	if binjnodes:
		[simfile] = [cache_entry.path for node in binjnodes for cache_entry in node.get_output_cache()]
		injargs = {"injection-file": simfile}
	else:
		injargs = {}
	seglist = split_segment(seg, min_segment_length, pad, overlap, short_segment_duration, max_job_length)
	if verbose:
		print("Segment split: " + str(seglist), file=sys.stderr)
	nodes = set()
	for seg in seglist:
		nodes |= make_string_fragment(dag, datafindnodes | binjnodes, instrument, seg, tag, framecache, injargs = injargs)
	return nodes


#
# all segments
#


def make_single_instrument_stage(dag, datafinds, seglistdict, tag, min_segment_length, pad, overlap, short_segment_duration, max_job_length, binjnodes = set(), verbose = False):
	nodes = set()
	for instrument, seglist in seglistdict.items():
		for seg in seglist:
			if verbose:
				print("generating %s fragment %s" % (instrument, str(seg)), file=sys.stderr)

			# find the datafind job this job is going to need
			dfnodes = set([node for node in datafinds if (node.get_ifo() == instrument) and (seg in segments.segment(node.get_start(), node.get_end()))])
			if len(dfnodes) != 1:
				raise ValueError("error, not exactly 1 datafind is suitable for trigger generator job at %s in %s" % (str(seg), instrument))

			# trigger generator jobs
			nodes |= make_string_segment_fragment(dag, dfnodes, instrument, seg, tag, min_segment_length, pad, overlap, short_segment_duration, max_job_length, binjnodes = binjnodes, verbose = verbose)

	# done
	return nodes


#
# =============================================================================
#
#                           lalapps_run_sqlite Jobs
#
# =============================================================================
#


def write_clip_segment_sql_file(filename):
	code = """DELETE FROM
	segment
WHERE
	(end_time + 1e-9 * end_time_ns < (SELECT MIN(in_start_time + 1e-9 * in_start_time_ns) FROM search_summary NATURAL JOIN process WHERE program == 'StringSearch'))
	OR
	(start_time + 1e-9 * start_time_ns > (SELECT MAX(in_end_time + 1e-9 * in_end_time_ns) FROM search_summary NATURAL JOIN process WHERE program == 'StringSearch'));

VACUUM;"""

	print(code, file=file(filename, "w"))

	return filename


def make_run_sqlite_fragment(dag, parents, tag, sql_file, files_per_run_sqlite = None):
	if files_per_run_sqlite is None:
		files_per_run_sqlite = runsqlitejob.files_per_run_sqlite
	nodes = set()
	input_cache = power.collect_output_caches(parents)
	while input_cache:
		node = RunSqliteNode(runsqlitejob)
		node.set_sql_file(sql_file)
		node.add_input_cache([cache_entry for cache_entry, parent in input_cache[:files_per_run_sqlite]])
		for parent in set(parent for cache_entry, parent in input_cache[:files_per_run_sqlite]):
			node.add_parent(parent)
		del input_cache[:files_per_run_sqlite]
		seg = power.cache_span(node.get_output_cache())
		node.set_name("lalapps_run_sqlite_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
		dag.add_node(node)
		nodes.add(node)
	return nodes


#
# =============================================================================
#
#                          lalapps_string_meas_likelihood Jobs
#
# =============================================================================
#


def make_meas_likelihood_fragment(dag, parents, tag, files_per_meas_likelihood = None):
	if files_per_meas_likelihood is None:
		files_per_meas_likelihood = meas_likelihoodjob.files_per_meas_likelihood
	nodes = set()
	input_cache = power.collect_output_caches(parents)
	while input_cache:
		node = MeasLikelihoodNode(meas_likelihoodjob)
		node.add_input_cache([cache_entry for cache_entry, parent in input_cache[:files_per_meas_likelihood]])
		for parent in set(parent for cache_entry, parent in input_cache[:files_per_meas_likelihood]):
			node.add_parent(parent)
		del input_cache[:files_per_meas_likelihood]
		seg = power.cache_span(node.get_input_cache())
		node.set_name("lalapps_string_meas_likelihood_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
		node.set_output(tag)
		dag.add_node(node)
		nodes.add(node)
	return nodes


#
# =============================================================================
#
#                          lalapps_string_calc_likelihood Jobs
#
# =============================================================================
#


def make_calc_likelihood_fragment(dag, parents, likelihood_parents, tag, files_per_calc_likelihood = None, verbose = False):
  if files_per_calc_likelihood is None:
    files_per_calc_likelihood = calc_likelihoodjob.files_per_calc_likelihood
  input_cache = power.collect_output_caches(parents)
  likelihood_cache = power.collect_output_caches(likelihood_parents)
  nodes = set()
  while input_cache:
    node = CalcLikelihoodNode(calc_likelihoodjob)
    node.add_input_cache([cache_entry for cache_entry, parent in input_cache[:files_per_calc_likelihood]])
    for parent in set(parent for cache_entry, parent in input_cache[:files_per_calc_likelihood]):
      node.add_parent(parent)
    del input_cache[:files_per_calc_likelihood]
    seg = power.cache_span(node.get_input_cache())
    node.set_name("lalapps_string_calc_likelihood_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
    for cache_entry, parent in likelihood_cache:
      node.add_parent(parent)
      node.add_likelihood_cache([cache_entry])
    dag.add_node(node)
    nodes.add(node)
  return nodes
