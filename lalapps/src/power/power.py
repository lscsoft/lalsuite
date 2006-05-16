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
import time

from glue.lal import CacheEntry
from glue import pipeline


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
	def __init__(self, cp):
		"""
		cp = ConfigParser object from which options are read.
		"""
		self.__executable = cp.get("condor", "lalapps_binj")
		self.__universe = cp.get("condor", "universe")
		pipeline.CondorDAGJob.__init__(self, self.__universe, self.__executable)
		pipeline.AnalysisJob.__init__(self, cp)

		self.add_ini_opts(cp, "lalapps_binj")

		self.set_stdout_file("logs/binj-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out")
		self.set_stderr_file("logs/binj-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err")
		self.set_sub_file("lalapps_binj.sub")


class BurstInjNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
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
		try:
			self.__usertag = job.get_config("lalapps_binj", "user-tag")
			self.add_var_opt("user-tag", self.__usertag)
		except:
			self.__usertag = None

	def set_user_tag(self, tag):
		self.__usertag = tag
		self.add_var_opt("user-tag", self.__usertag)

	def get_user_tag(self):
		return self.__usertag

	def get_output(self):
		"""
		Returns the file name of output from the power code. This
		must be kept synchronized with the name of the output file
		in binj.c.  Note in particular the calculation of the
		"start" and "duration" parts of the name.
		"""
		if not self.get_start() or not self.get_end():
			raise ValueError, "Start time or end time has not been set"

		if self.__usertag:
			return "HL-INJECTIONS_%s-%d-%d.xml" % (self.__usertag, int(self.get_start()), int(self.get_end() - self.get_start()))
		else:
			return "HL-INJECTIONS-%d-%d.xml" % (int(self.get_start()), int(self.get_end() - self.get_start()))


class PowerJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
	"""
	A lalapps_power job used by the power pipeline. The static options
	are read from the [lalapps_power] section in the ini file. The
	stdout and stderr from the job are directed to the logs directory.
	The job runs in the universe specified in the ini file. The path to
	the executable is determined from the ini file.
	"""
	def __init__(self, out_dir, cp):
		"""
		cp = ConfigParser object from which options are read.
		"""
		self.__executable = cp.get("condor", "lalapps_power")
		self.__universe = cp.get("condor", "universe")
		pipeline.CondorDAGJob.__init__(self, self.__universe, self.__executable)
		pipeline.AnalysisJob.__init__(self, cp)

		self.add_ini_opts(cp, "lalapps_power")

		self.set_stdout_file(os.path.join(out_dir, "power-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(out_dir, "power-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err"))
		self.set_sub_file("lalapps_power.sub")


class PowerNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
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
		try:
			self.__usertag = job.get_config("pipeline", "user-tag")
			self.add_var_opt("user-tag", self.__usertag)
		except:
			self.__usertag = None

	def set_user_tag(self, tag):
		self.__usertag = tag
		self.add_var_opt("user-tag", self.__usertag)

	def get_user_tag(self):
		return self.__usertag

	def get_output(self):
		"""
		Returns the file name of output from the power code. This
		must be kept synchronized with the name of the output file
		in power.c.  Note in particular the calculation of the
		"start" and "duration" parts of the name.
		"""
		if None in [self.get_start(), self.get_end(), self.get_ifo(), self.__usertag]:
			raise ValueError, "Start time, end time, ifo, or user tag has not been set"
		return "%s-POWER_%s-%d-%d.xml" % (self.get_ifo(), self.__usertag, int(self.get_start()), int(self.get_end()) - int(self.get_start()))

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


class BuclusterJob(pipeline.CondorDAGJob):
	def __init__(self, out_dir, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", config_parser.get("condor", "ligolw_bucluster"))
		self.set_sub_file("ligolw_bucluster.sub")
		self.set_stdout_file(os.path.join(out_dir, "ligolw_bucluster-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(out_dir, "ligolw_bucluster-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_bucluster")


class BuclusterNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
	pass


class BinjfindJob(pipeline.CondorDAGJob):
	def __init__(self, out_dir, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", config_parser.get("condor", "ligolw_binjfind"))
		self.set_sub_file("ligolw_binjfind.sub")
		self.set_stdout_file(os.path.join(out_dir, "ligolw_binjfind-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(out_dir, "ligolw_binjfind-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_binjfind")


class BinjfindNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
	pass


class TisiJob(pipeline.CondorDAGJob):
	def __init__(self, out_dir, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", config_parser.get("condor", "ligolw_tisi"))
		self.set_sub_file("ligolw_tisi.sub")
		self.set_stdout_file(os.path.join(out_dir, "ligolw_tisi-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(out_dir, "ligolw_tisi-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_tisi")


class TisiNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
	pass


class BurcaJob(pipeline.CondorDAGJob):
	def __init__(self, executable, out_dir, arguments):
		pipeline.CondorDAGJob.__init__(self, "vanilla", executable)
		self.set_sub_file("ligolw_burca.sub")
		self.set_stdout_file(os.path.join(out_dir, "ligolw_burca-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(out_dir, "ligolw_burca-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		if arguments:
			self.add_arg(arguments)


class BurcaNode(pipeline.CondorDAGNode):
	pass


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
buclusterjob = None

def init_job_types(cache_dir, out_dir, config_parser, types = ["datafind", "binj", "power", "lladd", "tisi", "binjfind", "bucluster"]):
	"""
	Construct definitions of the submit files.
	"""
	global datafindjob, binjjob, powerjob, lladdjob, tisijob, binjfindjob, buclusterjob, llb2mjob

	# LSCdataFind
	if "datafind" in types:
		datafindjob = pipeline.LSCDataFindJob(cache_dir, out_dir, config_parser)

	# lalapps_binj
	if "binj" in types:
		binjjob = BurstInjJob(config_parser)

	# lalapps_power
	if "power" in types:
		powerjob = PowerJob(out_dir, config_parser)

	# ligolw_add
	if "lladd" in types:
		lladdjob = pipeline.LigolwAddJob(out_dir, config_parser)

	# ligolw_tisi
	if "tisi" in types:
		tisijob = TisiJob(out_dir, config_parser)

	# ligolw_binjfind
	if "binjfind" in types:
		binjfindjob = BinjfindJob(out_dir, config_parser)

	# ligolw_bucluster
	if "bucluster" in types:
		buclusterjob = BuclusterJob(out_dir, config_parser)


#
# =============================================================================
#
#                                 Segmentation
#
# =============================================================================
#

def split_segment(powerjob, segment, psds):
	"""
	Split the data segment into correctly-overlaping segments.  We try
	to have the numbers of PSDs in each segment be equal to psds, but
	with a short segment at the end if needed.
	"""
	psd_length = float(powerjob.get_opts()["psd-average-points"]) / float(powerjob.get_opts()["resample-rate"])
	window_length = float(powerjob.get_opts()["window-length"]) / float(powerjob.get_opts()["resample-rate"])
	window_shift = float(powerjob.get_opts()["window-shift"]) / float(powerjob.get_opts()["resample-rate"])
	filter_corruption = float(powerjob.get_opts()["filter-corruption"]) / float(powerjob.get_opts()["resample-rate"])

	psd_overlap = window_length - window_shift

	joblength = psds * (psd_length - psd_overlap) + psd_overlap + 2 * filter_corruption
	joboverlap = 2 * filter_corruption + psd_overlap

	# can't use range() with non-integers
	segs = segments.segmentlist()
	t = segment[0]
	while t < segment[1] - joboverlap:
		segs.append(segments.segment(t, t + joblength) & segment)
		t += joblength - joboverlap
	return segs


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

def make_lladd_fragment(dag, cache_dir, parents, seg, tag):
	cache_name = os.path.join(cache_dir, "lladd-%s.cache" % tag)
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

def make_multipower_fragment(dag, cache_dir, powerparents, lladdparents, framecache, seglist, instrument, tag, injargs = {}):
	node = make_lladd_fragment(dag, cache_dir, [make_power_fragment(dag, powerparents, framecache, seg, instrument, tag, injargs) for seg in seglist] + lladdparents, seglist.extent(), "POWER_%s" % tag)
	node.set_output("%s-POWER_%s-%s-%s.xml" % (instrument, tag, int(seglist.extent()[0]), int(seglist.extent().duration())))
	return node


#
# =============================================================================
#
#             Analyze A Segment Using Multiple lalapps_power Jobs
#
# =============================================================================
#

def make_power_segment_fragment(dag, cache_dir, datafindnode, segment, instrument, psds_per_job, tag):
	"""
	Construct a DAG fragment for an entire segment, splitting the
	segment into multiple power jobs.
	"""
	seglist = split_segment(powerjob, segment, psds_per_job)
	print >>sys.stderr, "Segment split: " + str(seglist)

	lladdnode = make_multipower_fragment(dag, cache_dir, [datafindnode], [], datafindnode.get_output(), seglist, instrument, tag)

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

def make_multibinj_fragment(dag, cache_dir, tag, seg):
	flow = float(powerjob.get_opts()["low-freq-cutoff"])
	fhigh = flow + float(powerjob.get_opts()["bandwidth"])

	# do this many injections between flow and fhigh inclusively
	injection_bands = binjjob._AnalysisJob__cp.getint("pipeline", "injection_bands")

	# determine frequency ratio from number of injections across band
	# (take a bit off to allow fhigh to be picked up despite round-off
	# errors)
	fratio = 0.9999999 * (fhigh / flow) ** (1.0/injection_bands)

	binjnodes = [make_binj_fragment(dag, tag, seg, 0.0, flow, fhigh, fratio, injection_bands)]

	node = make_lladd_fragment(dag, cache_dir, binjnodes, seg, tag)
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
	node.set_output("TISI_%s.xml" % tag)
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

def make_injection_segment_fragment(dag, cache_dir, datafindnode, segment, instrument, calibration_cache, psds_per_job, tag):
	seglist = split_segment(powerjob, segment, psds_per_job)
	print >>sys.stderr, "Injections split: " + str(seglist)

	binjfrag = make_multibinj_fragment(dag, cache_dir, "INJECTIONS_%s" % tag, seglist.extent())

	tisifrag = make_tisi_fragment(dag, "INJECTIONS_%s" % tag)

	lladdnode = make_multipower_fragment(dag, cache_dir, [datafindnode, binjfrag], [binjfrag, tisifrag], datafindnode.get_output(), seglist, instrument, "INJECTIONS_%s" % tag, injargs = {"burstinjection-file": binjfrag.get_output(), "calibration-cache": calibration_cache})

	return lladdnode


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
	cluster.set_input(parent.get_output())
	cluster.set_output(parent.get_output())
	cluster.add_macro("macrocomment", tag)
	dag.add_node(cluster)

	binjfind = BinjfindNode(binjfindjob)
	binjfind.set_name("ligolw_binjfind-%s" % tag)
	binjfind.add_parent(cluster)
	binjfind.add_file_arg(cluster.get_output())
	binjfind.add_macro("macrocomment", tag)
	dag.add_node(binjfind)

	return binjfind
