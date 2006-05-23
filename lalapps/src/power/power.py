# $Id$

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
Classes needed for the excess power analysis pipeline.
"""

__author__ = "Duncan Brown <duncan@gravity.phys.uwm.edu>, Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]

import math
import os
import sys
import time

from glue import segments
from glue.lal import CacheEntry
from glue import pipeline


#
# =============================================================================
#
#                                   Helpers
#
# =============================================================================
#

def get_universe(config_parser):
	return config_parser.get("condor", "universe")

def get_executable(config_parser, name):
	return config_parser.get("condor", name)

def get_out_dir(config_parser):
	return config_parser.get("pipeline", "out_dir")

def get_cache_dir(config_parser):
	return config_parser.get("pipeline", "cache_dir")

def make_dag_directories(config_parser):
	os.mkdir(get_cache_dir(config_parser))
	os.mkdir(get_out_dir(config_parser))


#
# =============================================================================
#
#                            DAG Node and Job Class
#
# =============================================================================
#

class BurstInjJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
	"""
	A lalapps_binj job used by the power pipeline. The static options
	are read from the [lalapps_binj] section in the ini file. The
	stdout and stderr from the job are directed to the logs directory.
	The job runs in the universe specified in the ini file.  The path
	to the executable is determined from the ini file.
	"""
	def __init__(self, config_parser):
		"""
		config_parser = ConfigParser object from which options are
		read.
		"""
		pipeline.CondorDAGJob.__init__(self, get_universe(config_parser), get_executable(config_parser, "lalapps_binj"))
		pipeline.AnalysisJob.__init__(self, config_parser)

		self.add_ini_opts(config_parser, "lalapps_binj")

		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "binj-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "binj-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err"))
		self.set_sub_file("lalapps_binj.sub")


class BurstInjNode(pipeline.AnalysisNode):
	"""
	A BurstInjNode runs an instance of lalapps_binj.
	"""
	def __init__(self, job):
		"""
		job = A CondorDAGJob that can run an instance of
		lalapps_power.
		"""
		pipeline.CondorDAGNode.__init__(self, job)
		pipeline.AnalysisNode.__init__(self)
		self.__usertag = None

	def set_user_tag(self, tag):
		self.__usertag = tag
		self.add_var_opt("user-tag", self.__usertag)

	def get_user_tag(self):
		return self.__usertag

	def get_output_files(self):
		"""
		Returns the file name of output from the power code. This
		must be kept synchronized with the name of the output file
		in binj.c.  Note in particular the calculation of the
		"start" and "duration" parts of the name.
		"""
		if not self.get_start() or not self.get_end():
			raise ValueError, "start time or end time has not been set"

		if self.__usertag:
			filename = "HL-INJECTIONS_%s-%d-%d.xml" % (self.__usertag, int(self.get_start()), int(self.get_end() - self.get_start()))
		else:
			filename = "HL-INJECTIONS-%d-%d.xml" % (int(self.get_start()), int(self.get_end() - self.get_start()))
		return [filename]

	def get_output(self):
		# FIXME: use get_output_files() instead
		return self.get_output_files()[0]


class PowerJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
	"""
	A lalapps_power job used by the power pipeline. The static options
	are read from the [lalapps_power] and [lalapps_power_<inst>]
	sections in the ini file. The stdout and stderr from the job are
	directed to the logs directory.  The job runs in the universe
	specified in the ini file. The path to the executable is determined
	from the ini file.
	"""
	def __init__(self, config_parser):
		"""
		config_parser = ConfigParser object from which options are
		read.
		"""
		pipeline.CondorDAGJob.__init__(self, get_universe(config_parser), get_executable(config_parser, "lalapps_power"))
		pipeline.AnalysisJob.__init__(self, config_parser)

		self.add_ini_opts(config_parser, "lalapps_power")

		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "power-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "power-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err"))
		self.set_sub_file("lalapps_power.sub")


class PowerNode(pipeline.AnalysisNode):
	"""
	A PowerNode runs an instance of the power code in a Condor DAG.
	"""
	def __init__(self, job):
		"""
		job = A CondorDAGJob that can run an instance of
		lalapps_power.
		"""
		pipeline.CondorDAGNode.__init__(self, job)
		pipeline.AnalysisNode.__init__(self)
		self.__usertag = None

	def set_ifo(self, instrument):
		"""
		Load additional options from the per-instrument section in
		the config file.
		"""
		pipeline.AnalysisNode.set_ifo(self, instrument)
		for optvalue in self.job()._AnalysisJob__cp.items("lalapps_power_%s" % instrument):
			self.add_var_arg("--%s %s" % optvalue)

	def set_user_tag(self, tag):
		self.__usertag = tag
		self.add_var_opt("user-tag", self.__usertag)

	def get_user_tag(self):
		return self.__usertag

	def get_output_files(self):
		"""
		Returns the file name of output from the power code. This
		must be kept synchronized with the name of the output file
		in power.c.  Note in particular the calculation of the
		"start" and "duration" parts of the name.
		"""
		if None in [self.get_start(), self.get_end(), self.get_ifo(), self.__usertag]:
			raise ValueError, "start time, end time, ifo, or user tag has not been set"
		filename = "%s-POWER_%s-%d-%d.xml" % (self.get_ifo(), self.__usertag, int(self.get_start()), int(self.get_end()) - int(self.get_start()))
		return [filename]

	def get_output(self):
		# FIXME: use get_output_files() instead
		return self.get_output_files()[0]

	def set_mdccache(self, file):
		"""
		Set the LAL frame cache to to use. The frame cache is
		passed to the job with the --frame-cache argument.  @param
		file: calibration file to use.
		"""
		self.add_var_opt("mdc-cache", file)
		self.add_input_file(file)

	def set_burstinj(self, file):
		"""
		Set the LAL frame cache to to use. The frame cache is
		passed to the job with the --frame-cache argument.  @param
		file: calibration file to use.
		"""
		self.add_var_opt("burstinjection-file", file)
		self.add_input_file(file)

	def set_inspinj(self, file):
		"""
		Set the LAL frame cache to to use. The frame cache is
		passed to the job with the --frame-cache argument.  @param
		file: calibration file to use.
		"""
		self.add_var_opt("inspiralinjection-file", file)
		self.add_input_file(file)

	def set_siminj(self, file):
		"""
		Set the LAL frame cache to to use. The frame cache is
		passed to the job with the --frame-cache argument.  @param
		file: calibration file to use.
		"""
		self.add_var_opt("siminjection-file", file)
		self.add_input_file(file)


class BucutJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_bucut"))
		self.set_sub_file("ligolw_bucut.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_bucut-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_bucut-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_bucut")


class BucutNode(pipeline.CondorDAGNode):
	def add_file_arg(self, filename):
		pipeline.CondorDAGNode.add_file_arg(self, filename)
		self.add_output_file(filename)

	def get_output(self):
		# FIXME: use get_output_files() instead
		return self.get_output_files()[0]


class BuclusterJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_bucluster"))
		self.set_sub_file("ligolw_bucluster.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_bucluster-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_bucluster-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_bucluster")


class BuclusterNode(pipeline.CondorDAGNode):
	def add_file_arg(self, filename):
		pipeline.CondorDAGNode.add_file_arg(self, filename)
		self.add_output_file(filename)

	def get_output(self):
		# FIXME: use get_output_files() instead
		return self.get_output_files()[0]


class BinjfindJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_binjfind"))
		self.set_sub_file("ligolw_binjfind.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_binjfind-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_binjfind-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_binjfind")


class BinjfindNode(pipeline.CondorDAGNode):
	def add_file_arg(self, filename):
		pipeline.CondorDAGNode.add_file_arg(self, filename)
		self.add_output_file(filename)

	def get_output(self):
		# FIXME: use get_output_files() instead
		return self.get_output_files()[0]


class TisiJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_tisi"))
		self.set_sub_file("ligolw_tisi.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_tisi-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_tisi-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_tisi")


class TisiNode(pipeline.CondorDAGNode):
	def add_file_arg(self, filename):
		pipeline.CondorDAGNode.add_file_arg(self, filename)
		self.add_output_file(filename)

	def get_output(self):
		# FIXME: use get_output_files() instead
		return self.get_output_files()[0]


class BurcaJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_burca"))
		self.set_sub_file("ligolw_burca.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_burca-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_burca-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_burca")


class BurcaNode(pipeline.CondorDAGNode):
	def add_file_arg(self, filename):
		pipeline.CondorDAGNode.add_file_arg(self, filename)
		self.add_output_file(filename)

	def get_output(self):
		# FIXME: use get_output_files() instead
		return self.get_output_files()[0]


#
# =============================================================================
#
#                                DAG Job Types
#
# =============================================================================
#

# This is *SUCH* a hack I don't know where to begin.  Please, shoot me.

datafindjob = None
binjjob = None
powerjob = None
lladdjob = None
tisijob = None
binjfindjob = None
bucutjob = None
buclusterjob = None
burcajob = None

def init_job_types(config_parser, types = ["datafind", "binj", "power", "lladd", "tisi", "binjfind", "bucluster", "bucut"]):
	"""
	Construct definitions of the submit files.
	"""
	global datafindjob, binjjob, powerjob, lladdjob, tisijob, binjfindjob, buclusterjob, llb2mjob, bucutjob

	# LSCdataFind
	if "datafind" in types:
		datafindjob = pipeline.LSCDataFindJob(get_cache_dir(config_parser), get_out_dir(config_parser), config_parser)

	# lalapps_binj
	if "binj" in types:
		binjjob = BurstInjJob(config_parser)

	# lalapps_power
	if "power" in types:
		powerjob = PowerJob(config_parser)

	# ligolw_add
	if "lladd" in types:
		lladdjob = pipeline.LigolwAddJob(get_out_dir(config_parser), config_parser)
		lladdjob.cache_dir = get_cache_dir(config_parser)

	# ligolw_tisi
	if "tisi" in types:
		tisijob = TisiJob(config_parser)

	# ligolw_binjfind
	if "binjfind" in types:
		binjfindjob = BinjfindJob(config_parser)

	# ligolw_bucut
	if "bucut" in types:
		bucutjob = BucutJob(config_parser)

	# ligolw_bucluster
	if "bucluster" in types:
		buclusterjob = BuclusterJob(config_parser)

	# ligolw_burca
	if "burca" in types:
		burcajob = BurcaJob(config_parser)


#
# =============================================================================
#
#                                 Segmentation
#
# =============================================================================
#

def split_segment(powerjob, segment, psds_per_job):
	"""
	Split the data segment into correctly-overlaping segments.  We try
	to have the numbers of PSDs in each segment be equal to
	psds_per_job, but with a short segment at the end if needed.
	"""
	psd_length = float(powerjob.get_opts()["psd-average-points"]) / float(powerjob.get_opts()["resample-rate"])
	window_length = float(powerjob.get_opts()["window-length"]) / float(powerjob.get_opts()["resample-rate"])
	window_shift = float(powerjob.get_opts()["window-shift"]) / float(powerjob.get_opts()["resample-rate"])
	filter_corruption = float(powerjob.get_opts()["filter-corruption"]) / float(powerjob.get_opts()["resample-rate"])

	psd_overlap = window_length - window_shift

	joblength = psds_per_job * (psd_length - psd_overlap) + psd_overlap + 2 * filter_corruption
	joboverlap = 2 * filter_corruption + psd_overlap

	# can't use range() with non-integers
	segs = segments.segmentlist()
	t = segment[0]
	while t < segment[1] - joboverlap:
		segs.append(segments.segment(t, t + joblength) & segment)
		t += joblength - joboverlap
	return segs


def test_segment(segment, psds_per_job):
	"""
	Return True if the segment can be analyzed using lalapps_power.
	"""
	# This is slow, but guaranteed to be using the same algorithm as
	# split_segment()
	return bool(split_segment(powerjob, segment, psds_per_job))


#
# =============================================================================
#
#                           LSCdataFind DAG Fragment
#
# =============================================================================
#

def make_datafind_fragment(dag, observatory, seg):
	datafind_pad = 512

	node = pipeline.LSCDataFindNode(datafindjob)
	node.set_name("LSCdataFind-%s-%s" % (int(seg[0]), int(seg.duration())))
	node.set_start(seg[0] - datafind_pad)
	node.set_end(seg[1] + 1)
	node.set_observatory(observatory)
	dag.add_node(node)

	return node


#
# =============================================================================
#
#                           ligolw_add DAG Fragment
#
# =============================================================================
#

def make_lladd_fragment(dag, parents, seg, tag):
	cache_name = os.path.join(lladdjob.cache_dir, "lladd-%s.cache" % tag)
	cachefile = file(cache_name, "w")

	node = pipeline.LigolwAddNode(lladdjob)
	node.set_name("lladd-%s" % tag)
	for parent in parents:
		node.add_parent(parent)
		cache = CacheEntry()
		cache.observatory = "ANY"
		cache.description = "EMPTY"
		cache.segment = seg
		cache.url = "file://localhost" + os.path.join(os.getcwd(), parent.get_output())
		print >>cachefile, str(cache)
	node.add_var_opt("input-cache", cache_name)
	dag.add_node(node)

	return node


#
# =============================================================================
#
#                          lalapps_power DAG Fragment
#
# =============================================================================
#

def make_power_fragment(dag, parents, framecache, seg, instrument, tag, injargs = {}):
	suffix = "%s-%s-%s" % (tag, int(seg[0]), int(seg.duration()))

	node = PowerNode(powerjob)
	node.set_name("lalapps_power-%s" % suffix)
	map(node.add_parent, parents)
	node.set_cache(framecache)
	node.set_ifo(instrument)
	node.set_start(seg[0])
	node.set_end(seg[1])
	node.set_user_tag(tag)
	for arg, value in injargs.iteritems():
		# this is a hack, but I can't be bothered
		node.add_var_arg("--%s %s" % (arg, value))
	dag.add_node(node)

	return node


#
# =============================================================================
#
#        DAG Fragment Combining Multiple lalapps_power With ligolw_add
#
# =============================================================================
#

def make_multipower_fragment(dag, powerparents, lladdparents, framecache, seglist, instrument, tag, injargs = {}):
	node = make_lladd_fragment(dag, [make_power_fragment(dag, powerparents, framecache, seg, instrument, tag, injargs) for seg in seglist] + lladdparents, seglist.extent(), "POWER_%s" % tag)
	node.set_output("%s-POWER_%s-%s-%s.xml" % (instrument, tag, int(seglist.extent()[0]), int(seglist.extent().duration())))
	return node


#
# =============================================================================
#
#             Analyze A Segment Using Multiple lalapps_power Jobs
#
# =============================================================================
#

def make_power_segment_fragment(dag, datafindnode, segment, instrument, psds_per_job, tag):
	"""
	Construct a DAG fragment for an entire segment, splitting the
	segment into multiple power jobs.
	"""
	seglist = split_segment(powerjob, segment, psds_per_job)
	print >>sys.stderr, "Segment split: " + str(seglist)

	lladdnode = make_multipower_fragment(dag, [datafindnode], [], datafindnode.get_output(), seglist, instrument, tag)

	return lladdnode


#
# =============================================================================
#
#                          lalapps_binj DAG Fragment
#
# =============================================================================
#

def make_binj_fragment(dag, tag, seg, offset, flow, fhigh, fratio, injection_bands):
	# one injection every time-step / pi seconds
	period = injection_bands * float(binjjob.get_opts()["time-step"]) / math.pi
	start = seg[0] - seg[0] % period + period * offset

	node = BurstInjNode(binjjob)
	node.set_start(start)
	node.set_end(seg[1])
	node.set_name("lalapps_binj-%d-%d" % (int(node.get_start()), int(flow)))
	node.set_user_tag(tag)
	node.add_macro("macroflow", flow)
	node.add_macro("macrofhigh", fhigh)
	node.add_macro("macrofratio", fratio)
	node.add_macro("macroseed", int(time.time() + node.get_start()))
	dag.add_node(node)

	return node


#
# =============================================================================
#
#         DAG Fragment Combining Multiple lalapps_binj With ligolw_add
#
# =============================================================================
#

def make_multibinj_fragment(dag, tag, seg):
	flow = float(powerjob.get_opts()["low-freq-cutoff"])
	fhigh = flow + float(powerjob.get_opts()["bandwidth"])

	# do this many injections between flow and fhigh inclusively
	injection_bands = binjjob._AnalysisJob__cp.getint("pipeline", "injection_bands")

	# determine frequency ratio from number of injections across band
	# (take a bit off to allow fhigh to be picked up despite round-off
	# errors)
	fratio = 0.9999999 * (fhigh / flow) ** (1.0/injection_bands)

	binjnodes = [make_binj_fragment(dag, tag, seg, 0.0, flow, fhigh, fratio, injection_bands)]

	node = make_lladd_fragment(dag, binjnodes, seg, tag)
	node.set_output("HL-%s-%d-%d.xml" % (tag, int(seg[0]), int(seg.duration())))

	return node


#
# =============================================================================
#
#                           ligolw_tisi DAG Fragment
#
# =============================================================================
#

def make_tisi_fragment(dag, tag):
	node = TisiNode(tisijob)
	node.set_name("ligolw_tisi-%s" % tag)
	node.add_file_arg("TISI_%s.xml" % tag)
	node.add_macro("macrocomment", tag)
	dag.add_node(node)

	return node


#
# =============================================================================
#
#                          ligolw_bucut DAG Fragment
#
# =============================================================================
#

def make_bucut_fragment(dag, tag, parent):
	node = BucutNode(bucutjob)
	node.set_name("ligolw_bucut-%s" % tag)
	node.add_parent(parent)
	node.add_file_arg(parent.get_output())
	node.add_macro("macrocomment", tag)
	dag.add_node(node)

	return node


#
# =============================================================================
#
#     Analyze A Segment Using Multiple lalapps_power Jobs With Injections
#
# =============================================================================
#

def make_injection_segment_fragment(dag, datafindnode, segment, instrument, psds_per_job, tag, binjfrag = None, tisifrag = None):
	seglist = split_segment(powerjob, segment, psds_per_job)
	print >>sys.stderr, "Injections split: " + str(seglist)

	if not binjfrag:
		binjfrag = make_multibinj_fragment(dag, "INJECTIONS_%s" % tag, seglist.extent())

	if not tisifrag:
		tisifrag = make_tisi_fragment(dag, "INJECTIONS_%s" % tag)

	lladdnode = make_multipower_fragment(dag, [datafindnode, binjfrag], [binjfrag, tisifrag], datafindnode.get_output(), seglist, instrument, "INJECTIONS_%s" % tag, injargs = {"burstinjection-file": binjfrag.get_output()})

	bucutnode = make_bucut_fragment(dag, "INJECTIONS_%s" % tag, lladdnode)

	return bucutnode


#
# =============================================================================
#
#                         ligolw_binjfind DAG Fragment
#
# =============================================================================
#

def make_binjfind_fragment(dag, parent, tag):
	cluster = BuclusterNode(buclusterjob)
	cluster.set_name("ligolw_bucluster-%s" % tag)
	cluster.add_parent(parent)
	cluster.add_file_arg(parent.get_output())
	cluster.add_macro("macrocomment", tag)
	dag.add_node(cluster)

	binjfind = BinjfindNode(binjfindjob)
	binjfind.set_name("ligolw_binjfind-%s" % tag)
	binjfind.add_parent(cluster)
	binjfind.add_file_arg(cluster.get_output())
	binjfind.add_macro("macrocomment", tag)
	dag.add_node(binjfind)

	return binjfind
