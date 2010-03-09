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


import math
import os
import sys


from glue import iterutils
from glue import pipeline
from glue import segments
from glue.lal import CacheEntry, LIGOTimeGPS
from pylal import llwapp
from lalapps import power


__author__ = 'Xavier Siemens<siemens@gravity.phys.uwm.edu>'
__date__ = '$Date$'[7:-2]
__version__ = '$Revision$'[11:-2]


#
# =============================================================================
#
#                            DAG Node and Job Class
#
# =============================================================================
#


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
    self.set_sub_file("lalapps_StringSearch.sub")
    self.add_condor_cmd("Requirements", "Memory > 800")


class StringNode(pipeline.AnalysisNode):
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

  def set_ifo(self, instrument):
    """
    Load additional options from the per-instrument section in
    the config file.
    """
    if self.output_cache:
      raise AttributeError, "cannot change attributes after computing output cache"
    pipeline.AnalysisNode.set_ifo(self, instrument)
    for optvalue in self.job()._AnalysisJob__cp.items("lalapps_StringSearch_%s" % instrument):
      self.add_var_arg("--%s %s" % optvalue)

  def set_user_tag(self, tag):
    if self.output_cache:
      raise AttributeError, "cannot change attributes after computing output cache"
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
        raise ValueError, "start time, end time, ifo, or user tag has not been set"
      seg = segments.segment(LIGOTimeGPS(self.get_start()), LIGOTimeGPS(self.get_end()))
      self.set_output("triggers/%s-STRINGSEARCH_%s-%d-%d.xml.gz" % (self.get_ifo(), self.__usertag, int(self.get_start()), int(self.get_end()) - int(self.get_start())))

    return self._AnalysisNode__output

  def set_injection_file(self, file):
    """
    Set the name of the XML file from which to read a list of
    software injections.
    """
    self.add_var_opt("injection-file", file)
    self.add_input_file(file)


#
# =============================================================================
#
#                                 Segmentation
#
# =============================================================================
#


def clip_segment(seg, pad, short_segment_duration):
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
  duration = float(abs(seg)) - 2 * pad
  extra = (2 * duration + short_segment_duration) % (4 * short_segment_duration)
  extra /= 2

  # clip segment
  seg = segments.segment(seg[0], seg[1] - extra)

  # bounds must be integers
  if abs((int(seg[0]) - seg[0]) / seg[0]) > 1e-14 or abs((int(seg[1]) - seg[1]) / seg[1]) > 1e-14:
    raise ValueError, "segment %s does not have integer boundaries" % str(seg)
  seg = segments.segment(int(seg[0]), int(seg[1]))

  # done
  return seg


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
  # the speed of the llwapp.get_coincident_segmentlistdict()
  # function when the input segmentlists have had many data
  # quality holes poked out of them
  remove_too_short_segments(seglists, min_segment_length, pad)

  # extract the segments that are coincident under the time
  # slides
  new = llwapp.get_coincident_segmentlistdict(seglists, offset_vectors)

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


def init_job_types(config_parser, job_types = ("string",)):
  """
  Construct definitions of the submit files.
  """
  global stringjob

  # lalapps_StringSearch
  if "string" in job_types:
    stringjob = StringJob(config_parser)


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


def split_segment(seg, min_segment_length, pad, overlap, short_segment_duration):
	seglist = segments.segmentlist()
	while abs(seg) >= min_segment_length + 2 * pad:
		# try to use 5* min_segment_length each time (must be an odd
		# multiple!).
		if abs(seg) >= 5 * min_segment_length + 2 * pad:
			seglist.append(segments.segment(seg[0], seg[0] + 5 * min_segment_length + 2 * pad))
		else:
			seglist.append(clip_segment(seg, pad, short_segment_duration))
		# advance segment
		seg = segments.segment(seglist[-1][1] - overlap, seg[1])
	return seglist


def make_string_segment_fragment(dag, datafindnodes, instrument, segment, tag, min_segment_length, pad, overlap, short_segment_duration, binjnodes = set(), verbose = False):
	"""
	Construct a DAG fragment for an entire segment, splitting the
	segment into multiple trigger generator jobs.
	"""
	# only one frame cache file can be provided as input, and only one
	# injection description file can be provided as input.
	# the unpacking indirectly tests that the file count is correct
	[framecache] = [node.get_output() for node in datafindnodes]
	if binjnodes:
		[simfile] = [cache_entry.path() for node in binjnodes for cache_entry in node.get_output_cache()]
		injargs = {"injection-file": simfile}
	else:
		injargs = {}
	seglist = split_segment(segment, min_segment_length, pad, overlap, short_segment_duration)
	if verbose:
		print >>sys.stderr, "Segment split: " + str(seglist)
	nodes = set()
	for seg in seglist:
		nodes |= make_string_fragment(dag, datafindnodes | binjnodes, instrument, seg, tag, framecache, injargs = injargs)
	return nodes


#
# all segments
#


def make_single_instrument_stage(dag, datafinds, seglistdict, tag, min_segment_length, pad, overlap, short_segment_duration, binjnodes = set(), verbose = False):
	nodes = set()
	for instrument, seglist in seglistdict.items():
		for seg in seglist:
			if verbose:
				print >>sys.stderr, "generating %s fragment %s" % (instrument, str(seg))

			# find the datafind job this job is going to need
			dfnodes = set([node for node in datafinds if (node.get_ifo() == instrument) and (seg in segments.segment(node.get_start(), node.get_end()))])
			if len(dfnodes) != 1:
				raise ValueError, "error, not exactly 1 datafind is suitable for trigger generator job at %s in %s" % (str(seg), instrument)

			# trigger generator jobs
			nodes |= make_string_segment_fragment(dag, dfnodes, instrument, seg, tag, min_segment_length, pad, overlap, short_segment_duration, binjnodes = binjnodes, verbose = verbose)

	# done
	return nodes
