#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Excess power pipeline construction tools.
"""


import errno
import os
import sys
import time


import igwn_segments as segments
from igwn_segments import utils as segmentsUtils
import lal
from lal import iterutils
from lal import pipeline
from lal.utils import CacheEntry
import lalburst
from . import cafe


__author__ = "Duncan Brown <duncan@gravity.phys.uwm.edu>, Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"


#
# =============================================================================
#
#                                   Helpers
#
# =============================================================================
#


def get_universe(config_parser):
	return config_parser.get("condor", "universe")


def get_accounting_group(config_parser):
	return config_parser.get("condor", "accounting_group")


def get_executable(config_parser, name):
	return config_parser.get("condor", name)


def get_out_dir(config_parser):
	return config_parser.get("pipeline", "out_dir")


def get_cache_dir(config_parser):
	return config_parser.get("pipeline", "cache_dir")


def get_triggers_dir(config_parser):
	return config_parser.get("pipeline", "triggers_dir")


def make_dir_if_not_exists(dir):
	try:
		os.mkdir(dir)
	except OSError as e:
		if e.errno != errno.EEXIST:
			# OK if directory exists, otherwise report error
			raise e


def make_dag_directories(config_parser):
	make_dir_if_not_exists(get_cache_dir(config_parser))
	make_dir_if_not_exists(get_out_dir(config_parser))


def get_files_per_bucluster(config_parser):
	return config_parser.getint("pipeline", "files_per_bucluster")


def get_files_per_bucut(config_parser):
	return config_parser.getint("pipeline", "files_per_bucut")


def get_files_per_burca(config_parser):
	return config_parser.getint("pipeline", "files_per_burca")


def get_files_per_binjfind(config_parser):
	return config_parser.getint("pipeline", "files_per_binjfind")


class TimingParameters(object):
	"""
	A class to hold timing parameter values.
	"""
	def __init__(self, config_parser):
		# initialize from config file
		self.resample_rate = config_parser.getfloat("lalapps_power", "resample-rate")
		self.window_length = config_parser.getint("lalapps_power", "window-length")
		self.max_tile_length = int(config_parser.getfloat("lalapps_power", "max-tile-duration") * self.resample_rate)
		self.tile_stride_fraction = config_parser.getfloat("lalapps_power", "tile-stride-fraction")
		self.filter_corruption = config_parser.getint("lalapps_power", "filter-corruption")
		self.max_tile_bandwidth = config_parser.getfloat("lalapps_power", "max-tile-bandwidth")

		# populate additional computed parameters from library code
		self.psd_length, self.psd_shift, self.window_shift, self.window_pad, self.tiling_length = lalburst.EPGetTimingParameters(
			self.window_length,
			self.max_tile_length,
			self.tile_stride_fraction,
			config_parser.getint("lalapps_power", "psd-average-points")
		)


def make_cache_entry(input_cache, description, path):
	# summarize segment information
	seglists = segments.segmentlistdict()
	for c in input_cache:
		seglists |= c.segmentlistdict

	# obtain instrument list
	instruments = seglists.keys()
	if None in instruments:
		instruments.remove(None)
	instruments.sort()

	# remove empty segment lists to allow extent_all() to work
	for instrument in seglists.keys():
		if not seglists[instrument]:
			del seglists[instrument]

	# make the URL
	if path:
		url = "file://localhost%s" % os.path.abspath(path)
	else:
		# FIXME:  old version of CacheEntry allowed None for URL,
		# new version doesn't.  correct fix is to modify calling
		# code to not try to initialize the output cache until
		# after the input is known, but for now we'll just do this
		# stupid hack.
		url = "file://localhost/dev/null"

	# construct a cache entry from the instruments and
	# segments that remain
	return CacheEntry("+".join(instruments) or None, description, seglists.extent_all(), url)


def collect_output_caches(parents):
	cache = [(cache_entry, parent) for parent in parents for cache_entry in parent.get_output_cache()]
	cache.sort(key = lambda x: x[0].segment)
	return cache


def match_nodes_to_caches(nodes, caches):
	"""
	For each cache, get the set of nodes whose output files it
	contains.  A node is allowed to provide more than one output file,
	and thus can be listed in more than one set.
	"""
	# cache_entry --> node loop-up table
	nodes = set(nodes)
	index = {}
	for node in nodes:
		for cache_entry in node.get_output_cache():
			index[cache_entry] = node

	# can't use [set()] * len(caches) for the normal reason
	node_groups = [set() for cache in caches]

	# form node groups matching input caches
	for node_group, cache in zip(node_groups, caches):
		for cache_entry in cache:
			node_group.add(index[cache_entry])

	# how many nodes didn't get used?
	unused = len(nodes) - len(set.union(*node_groups))

	# done
	return node_groups, unused


def cache_span(cache):
	a = min([cache_entry.segment[0] for cache_entry in cache])
	b = max([cache_entry.segment[1] for cache_entry in cache])
	return segments.segment(a, b)


#
# How to write an output cache
#


def write_output_cache(nodes, filename):
	f = file(filename, "w")
	for cache_entry, node in collect_output_caches(nodes):
		print(str(cache_entry), file=f)


#
# =============================================================================
#
#                            DAG Node and Job Class
#
# =============================================================================
#


class RMJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		"""
		config_parser = ConfigParser object
		"""
		pipeline.CondorDAGJob.__init__(self, "local", "/bin/rm")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "rm-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "rm-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", get_accounting_group(config_parser))
		self.add_opt("force", "")
		self.set_sub_file("rm.sub")


class RMNode(pipeline.CondorDAGNode):
	def __init__(self, job):
		pipeline.CondorDAGNode.__init__(self, job)
		self.input_cache = set()
		self.output_cache = set()
		self._CondorDAGNode__macros["initialdir"] = os.getcwd()

	def add_input_cache(self, cache):
		self.input_cache |= cache
		for cache_entry in cache:
			pipeline.CondorDAGNode.add_file_arg(self, cache_entry.path)

	def get_output_cache(self):
		return set()


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
		config_parser = ConfigParser object
		"""
		pipeline.CondorDAGJob.__init__(self, get_universe(config_parser), get_executable(config_parser, "lalapps_binj"))
		pipeline.AnalysisJob.__init__(self, config_parser)

		# do this many injections between flow and fhigh inclusively
		if config_parser.has_option("pipeline", "injection_bands"):
			self.injection_bands = config_parser.getint("pipeline", "injection_bands")
		else:
			self.injection_bands = None

		self.add_ini_opts(config_parser, "lalapps_binj")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "lalapps_binj-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "lalapps_binj-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", get_accounting_group(config_parser))
		self.set_sub_file("lalapps_binj.sub")

		self.output_dir = "."

		# one injection every time-step seconds
		self.time_step = config_parser.getfloat("lalapps_binj", "time-step")


class BurstInjNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
	def __init__(self, job):
		pipeline.CondorDAGNode.__init__(self, job)
		pipeline.AnalysisNode.__init__(self)
		self.__usertag = None
		self.output_cache = []
		self.output_dir = os.path.join(os.getcwd(), self.job().output_dir)
		self._CondorDAGNode__macros["initialdir"] = os.getcwd()

	def set_user_tag(self, tag):
		self.__usertag = tag
		self.add_var_opt("user-tag", self.__usertag)

	def get_user_tag(self):
		if self.output_cache:
			raise AttributeError("cannot change attributes after computing output cache")
		return self.__usertag

	def set_time_slide_file(self, filename):
		self.add_var_opt("time-slide-file", filename)

	def get_time_slide_file(self):
		return self.get_opts().get("macrotimeslidefile", None)

	def set_start(self, start):
		if self.output_cache:
			raise AttributeError("cannot change attributes after computing output cache")
		self.add_var_opt("gps-start-time", start)

	def get_start(self):
		return self.get_opts().get("macrogpsstarttime", None)

	def set_end(self, end):
		if self.output_cache:
			raise AttributeError("cannot change attributes after computing output cache")
		self.add_var_opt("gps-end-time", end)

	def get_end(self):
		return self.get_opts().get("macrogpsendtime", None)

	def get_output_cache(self):
		"""
		Returns a LAL cache of the output file name.  Calling this
		method also induces the output name to get set, so it must
		be at least once.
		"""
		if not self.output_cache:
			# FIXME:  instruments hardcoded to "everything"
			self.output_cache = [CacheEntry("G1+H1+H2+L1+T1+V1", self.__usertag, segments.segment(lal.LIGOTimeGPS(self.get_start()), lal.LIGOTimeGPS(self.get_end())), "file://localhost" + os.path.abspath(self.get_output()))]
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		if self._AnalysisNode__output is None:
			if None in (self.get_start(), self.get_end(), self.__usertag):
				raise ValueError("start time, end time, ifo, or user tag has not been set")
			seg = segments.segment(lal.LIGOTimeGPS(self.get_start()), lal.LIGOTimeGPS(self.get_end()))
			self.set_output(os.path.join(self.output_dir, "G1+H1+H2+L1+T1+V1-INJECTIONS_%s-%d-%d.xml.gz" % (self.__usertag, int(self.get_start()), int(self.get_end() - self.get_start()))))
		return self._AnalysisNode__output


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
		config_parser = ConfigParser object
		"""
		pipeline.CondorDAGJob.__init__(self, get_universe(config_parser), get_executable(config_parser, "lalapps_power"))
		pipeline.AnalysisJob.__init__(self, config_parser)
		self.add_ini_opts(config_parser, "lalapps_power")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "lalapps_power-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "lalapps_power-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", get_accounting_group(config_parser))
		self.set_sub_file("lalapps_power.sub")

		self.output_dir = "."


class PowerNode(pipeline.AnalysisNode):
	def __init__(self, job):
		pipeline.CondorDAGNode.__init__(self, job)
		pipeline.AnalysisNode.__init__(self)
		self.__usertag = None
		self.output_cache = []
		self.output_dir = os.path.join(os.getcwd(), self.job().output_dir)
		self._CondorDAGNode__macros["initialdir"] = os.getcwd()

	def set_ifo(self, instrument):
		"""
		Load additional options from the per-instrument section in
		the config file.
		"""
		if self.output_cache:
			raise AttributeError("cannot change attributes after computing output cache")
		pipeline.AnalysisNode.set_ifo(self, instrument)
		for optvalue in self.job()._AnalysisJob__cp.items("lalapps_power_%s" % instrument):
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
			self.output_cache = [CacheEntry(self.get_ifo(), self.__usertag, segments.segment(lal.LIGOTimeGPS(self.get_start()), lal.LIGOTimeGPS(self.get_end())), "file://localhost" + os.path.abspath(self.get_output()))]
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		if self._AnalysisNode__output is None:
			if None in (self.get_start(), self.get_end(), self.get_ifo(), self.__usertag):
				raise ValueError("start time, end time, ifo, or user tag has not been set")
			seg = segments.segment(lal.LIGOTimeGPS(self.get_start()), lal.LIGOTimeGPS(self.get_end()))
			self.set_output(os.path.join(self.output_dir, "%s-POWER_%s-%d-%d.xml.gz" % (self.get_ifo(), self.__usertag, int(self.get_start()), int(self.get_end()) - int(self.get_start()))))
		return self._AnalysisNode__output

	def set_mdccache(self, file):
		"""
		Set the LAL frame cache to to use. The frame cache is
		passed to the job with the --frame-cache argument.  @param
		file: calibration file to use.
		"""
		self.add_var_opt("mdc-cache", file)
		self.add_input_file(file)

	def set_injection_file(self, file):
		"""
		Set the name of the XML file from which to read a list of
		software injections.
		"""
		self.add_var_opt("injection-file", file)
		self.add_input_file(file)


class LigolwAddNode(pipeline.LigolwAddNode):
	def __init__(self, job, remove_input, *args):
		pipeline.LigolwAddNode.__init__(self, job, *args)
		self.input_cache = []
		self.output_cache = []
		self.cache_dir = os.path.join(os.getcwd(), self.job().cache_dir)
		self.output_dir = os.path.join(os.getcwd(), ".")	# "." == self.job().output_dir except the job class doesn't yet have this info
		self._CondorDAGNode__macros["initialdir"] = os.getcwd()
		self.remove_input = bool(remove_input)
		if self.remove_input:
			self.add_var_arg("--remove-input")

	def __update_output_cache(self, observatory = None, segment = None):
		del self.output_cache[:]
		cache_entry = make_cache_entry(self.input_cache, None, self._AnalysisNode__output)
		if observatory is not None:
			cache_entry.observatory = observatory
		if segment is not None:
			cache_entry.segment = segment
		cache_entry = self.output_cache.append(cache_entry)

	def set_name(self, *args):
		pipeline.LigolwAddNode.set_name(self, *args)
		self.cache_name = os.path.join(self.cache_dir, "%s.cache" % self.get_name())
		self.add_var_opt("input-cache", self.cache_name)

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)
		self.__update_output_cache()

	def set_output(self, path = None, observatory = None, segment = None):
		pipeline.LigolwAddNode.set_output(self, path)
		self.__update_output_cache(observatory = observatory, segment = segment)

	def add_preserve_cache(self, cache):
		if self.remove_input:
			for c in cache:
				self.add_var_arg("--remove-input-except %s" % c.path)

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def write_input_files(self, *args):
		f = file(self.cache_name, "w")
		for c in self.input_cache:
			print(str(c), file=f)
		pipeline.LigolwAddNode.write_input_files(self, *args)

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class BucutJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "lalburst_cut"))
		self.set_sub_file("lalburst_cut.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "lalburst_cut-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "lalburst_cut-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", get_accounting_group(config_parser))
		self.add_condor_cmd("Requirements", "Memory > 1100")
		self.add_ini_opts(config_parser, "lalburst_cut")

		self.files_per_bucut = get_files_per_bucut(config_parser)
		if self.files_per_bucut < 1:
			raise ValueError("files_per_bucut < 1")


class BucutNode(pipeline.CondorDAGNode):
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

	def add_file_arg(self, filename):
		raise NotImplementedError

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class BuclusterJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "lalburst_cluster"))
		self.set_sub_file("lalburst_cluster.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "lalburst_cluster-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "lalburst_cluster-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", get_accounting_group(config_parser))
		self.add_condor_cmd("Requirements", "Memory > 1100")
		self.add_ini_opts(config_parser, "lalburst_cluster")

		self.cache_dir = get_cache_dir(config_parser)

		self.files_per_bucluster = get_files_per_bucluster(config_parser)
		if self.files_per_bucluster < 1:
			raise ValueError("files_per_bucluster < 1")


class BuclusterNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = self.input_cache
		self.cache_dir = os.path.join(os.getcwd(), self.job().cache_dir)
		self._CondorDAGNode__macros["initialdir"] = os.getcwd()

	def set_name(self, *args):
		pipeline.CondorDAGNode.set_name(self, *args)
		self.cache_name = os.path.join(self.cache_dir, "%s.cache" % self.get_name())
		self.add_var_opt("input-cache", self.cache_name)

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def write_input_files(self, *args):
		f = file(self.cache_name, "w")
		for c in self.input_cache:
			print(str(c), file=f)
		pipeline.CondorDAGNode.write_input_files(self, *args)

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class BinjfindJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "lalburst_injfind"))
		self.set_sub_file("lalburst_injfind.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "lalburst_injfind-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "lalburst_injfind-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", get_accounting_group(config_parser))
		self.add_ini_opts(config_parser, "lalburst_injfind")

		self.files_per_binjfind = get_files_per_binjfind(config_parser)
		if self.files_per_binjfind < 1:
			raise ValueError("files_per_binjfind < 1")


class BinjfindNode(pipeline.CondorDAGNode):
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

	def add_file_arg(self, filename):
		raise NotImplementedError

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class BurcaJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "lalburst_coinc"))
		self.set_sub_file("lalburst_coinc.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "lalburst_coinc-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "lalburst_coinc-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", get_accounting_group(config_parser))
		self.add_condor_cmd("Requirements", "Memory >= $(macrominram)")
		self.add_ini_opts(config_parser, "lalburst_coinc")

		self.files_per_burca = get_files_per_burca(config_parser)
		if self.files_per_burca < 1:
			raise ValueError("files_per_burca < 1")


class Burca2Job(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "lalburst_coinc"))
		self.set_sub_file("lalburst_coinc2.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "lalburst_coinc2-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "lalburst_coinc2-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", get_accounting_group(config_parser))
		self.add_ini_opts(config_parser, "lalburst_coinc2")

		self.cache_dir = get_cache_dir(config_parser)


class BurcaNode(pipeline.CondorDAGNode):
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
		longest_duration = max(abs(cache_entry.segment) for cache_entry in self.input_cache)
		if longest_duration > 25000:
			# ask for >= 1300 MB
			self.add_macro("macrominram", 1300)
		elif longest_duration > 10000:
			# ask for >= 800 MB
			self.add_macro("macrominram", 800)
		else:
			# run on any node
			self.add_macro("macrominram", 0)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError

	def set_coincidence_segments(self, seglist):
		self.add_var_arg("--coincidence-segments %s" % ",".join(segmentsUtils.to_range_strings(seglist)))


class SQLiteJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_sqlite"))
		self.set_sub_file("ligolw_sqlite.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_sqlite-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_sqlite-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", get_accounting_group(config_parser))
		self.add_ini_opts(config_parser, "ligolw_sqlite")


class SQLiteNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = []
		self._CondorDAGNode__macros["initialdir"] = os.getcwd()

	def add_input_cache(self, cache):
		if self.output_cache:
			raise AttributeError("cannot change attributes after computing output cache")
		self.input_cache.extend(cache)
		for c in cache:
			filename = c.path
			pipeline.CondorDAGNode.add_file_arg(self, filename)
			self.add_output_file(filename)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def set_output(self, filename):
		if self.output_cache:
			raise AttributeError("cannot change attributes after computing output cache")
		self.add_macro("macrodatabase", filename)

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		if not self.output_cache:
			self.output_cache = [make_cache_entry(self.input_cache, None, self.get_opts()["macrodatabase"])]
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class BurcaTailorJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "lalburst_power_meas_likelihood"))
		self.set_sub_file("lalburst_power_meas_likelihood.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "lalburst_power_meas_likelihood-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "lalburst_power_meas_likelihood-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("accounting_group", get_accounting_group(config_parser))
		self.add_ini_opts(config_parser, "lalburst_power_meas_likelihood")

		self.cache_dir = get_cache_dir(config_parser)
		self.output_dir = "."


class BurcaTailorNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = []
		self.cache_dir = os.path.join(os.getcwd(), self.job().cache_dir)
		self.output_dir = os.path.join(os.getcwd(), self.job().output_dir)
		self._CondorDAGNode__macros["initialdir"] = os.getcwd()

	def set_name(self, *args):
		pipeline.CondorDAGNode.set_name(self, *args)
		self.cache_name = os.path.join(self.cache_dir, "%s.cache" % self.get_name())

	def add_input_cache(self, cache):
		if self.output_cache:
			raise AttributeError("cannot change attributes after computing output cache")
		self.input_cache.extend(cache)
		for c in cache:
			filename = c.path
			pipeline.CondorDAGNode.add_file_arg(self, filename)
		self.add_output_file(filename)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def set_output(self, description):
		if self.output_cache:
			raise AttributeError("cannot change attributes after computing output cache")
		cache_entry = make_cache_entry(self.input_cache, description, "")
		filename = os.path.join(self.output_dir, "%s-%s-%d-%d.xml.gz" % (cache_entry.observatory, cache_entry.description, int(cache_entry.segment[0]), int(abs(cache_entry.segment))))
		self.add_var_opt("output", filename)
		cache_entry.url = "file://localhost" + os.path.abspath(filename)
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
		# oh.  my.  god.  this is fscked.
		for arg in self.get_args():
			if "--add-from-cache" in arg:
				f = file(self.cache_name, "w")
				for c in self.input_cache:
					print(str(c), file=f)
				pipeline.CondorDAGNode.write_input_files(self, *args)
				break

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


#
# =============================================================================
#
#                                DAG Job Types
#
# =============================================================================
#


#
# This is *SUCH* a hack I don't know where to begin.  Please, shoot me.
#


datafindjob = None
binjjob = None
powerjob = None
lladdjob = None
binjfindjob = None
bucutjob = None
buclusterjob = None
burcajob = None
burca2job = None
sqlitejob = None
burcatailorjob = None


def init_job_types(config_parser, job_types = ("datafind", "rm", "binj", "power", "lladd", "binjfind", "bucluster", "bucut", "burca", "burca2", "sqlite", "burcatailor")):
	"""
	Construct definitions of the submit files.
	"""
	global datafindjob, rmjob, binjjob, powerjob, lladdjob, binjfindjob, buclusterjob, llb2mjob, bucutjob, burcajob, burca2job, sqlitejob, burcatailorjob

	# ligo_data_find
	if "datafind" in job_types:
		datafindjob = pipeline.LSCDataFindJob(os.path.join(os.getcwd(), get_cache_dir(config_parser)), os.path.join(os.getcwd(), get_out_dir(config_parser)), config_parser)

	# rm
	if "rm" in job_types:
		rmjob = RMJob(config_parser)

	# lalapps_binj
	if "binj" in job_types:
		binjjob = BurstInjJob(config_parser)

	# lalapps_power
	if "power" in job_types:
		powerjob = PowerJob(config_parser)

	# ligolw_add
	if "lladd" in job_types:
		lladdjob = pipeline.LigolwAddJob(os.path.join(get_out_dir(config_parser)), config_parser)
		lladdjob.cache_dir = get_cache_dir(config_parser)

	# lalburst_injfind
	if "binjfind" in job_types:
		binjfindjob = BinjfindJob(config_parser)

	# lalburst_cut
	if "bucut" in job_types:
		bucutjob = BucutJob(config_parser)

	# lalburst_cluster
	if "bucluster" in job_types:
		buclusterjob = BuclusterJob(config_parser)

	# lalburst_coinc
	if "burca" in job_types:
		burcajob = BurcaJob(config_parser)

	# lalburst_coinc2
	if "burca2" in job_types:
		burca2job = Burca2Job(config_parser)

	# ligolw_sqlite
	if "sqlite" in job_types:
		sqlitejob = SQLiteJob(config_parser)

	# lalburst_power_meas_likelihood
	if "burcatailor" in job_types:
		burcatailorjob = BurcaTailorJob(config_parser)


#
# =============================================================================
#
#                                 Segmentation
#
# =============================================================================
#


def psds_from_job_length(timing_params, t):
	"""
	Return the number of PSDs that can fit into a job of length t
	seconds.  In general, the return value is a non-integer.
	"""
	if t < 0:
		raise ValueError(t)
	# convert to samples, and remove filter corruption
	t = t * timing_params.resample_rate - 2 * timing_params.filter_corruption
	if t < timing_params.psd_length:
		return 0
	return (t - timing_params.psd_length) / timing_params.psd_shift + 1


def job_length_from_psds(timing_params, psds):
	"""
	From the analysis parameters and a count of PSDs, return the length
	of the job in seconds.
	"""
	if psds < 1:
		raise ValueError(psds)
	# number of samples
	result = (psds - 1) * timing_params.psd_shift + timing_params.psd_length
	# add filter corruption
	result += 2 * timing_params.filter_corruption
	# convert to seconds
	return result / timing_params.resample_rate


def split_segment(timing_params, segment, psds_per_job):
	"""
	Split the data segment into correctly-overlaping segments.  We try
	to have the numbers of PSDs in each segment be equal to
	psds_per_job, but with a short segment at the end if needed.
	"""
	# in seconds
	joblength = job_length_from_psds(timing_params, psds_per_job)
	# in samples
	joboverlap = 2 * timing_params.filter_corruption + (timing_params.psd_length - timing_params.psd_shift)
	# in seconds
	joboverlap /= timing_params.resample_rate

	segs = segments.segmentlist()
	t = segment[0]
	while t + joblength <= segment[1]:
		segs.append(segments.segment(t, t + joblength) & segment)
		t += joblength - joboverlap

	extra_psds = int(psds_from_job_length(timing_params, float(segment[1] - t)))
	if extra_psds:
		segs.append(segments.segment(t, t + job_length_from_psds(timing_params, extra_psds)))
	return segs


def segment_ok(timing_params, segment):
	"""
	Return True if the segment can be analyzed using lalapps_power.
	"""
	return psds_from_job_length(timing_params, float(abs(segment))) >= 1.0


def remove_too_short_segments(seglistdict, timing_params):
	"""
	Remove segments from seglistdict that are too short to analyze.

	CAUTION:  this function modifies seglistdict in place.
	"""
	for seglist in seglistdict.values():
		iterutils.inplace_filter(lambda seg: segment_ok(timing_params, seg), seglist)


#
# =============================================================================
#
#                            Single Node Fragments
#
# =============================================================================
#


datafind_pad = 512


def make_datafind_fragment(dag, instrument, seg):
	node = pipeline.LSCDataFindNode(datafindjob)
	node.set_name("ligo_data_find-%s-%d-%d" % (instrument, int(seg[0]), int(abs(seg))))
	node.set_start(seg[0] - datafind_pad)
	node.set_end(seg[1] + 1)
	# FIXME: argh, I need the node to know what instrument it's for,
	# but can't call set_ifo() because that adds a --channel-name
	# command line argument (!?)
	node._AnalysisNode__ifo = instrument
	node.set_observatory(instrument[0])
	if node.get_type() is None:
		node.set_type(datafindjob.get_config_file().get("datafind", "type_%s" % instrument))
	node.set_retry(3)
	dag.add_node(node)
	return set([node])


def make_lladd_fragment(dag, parents, tag, segment = None, input_cache = None, remove_input = False, preserve_cache = None, extra_input_cache = None):
	node = LigolwAddNode(lladdjob, remove_input = remove_input)

	# link to parents
	for parent in parents:
		node.add_parent(parent)

	# build input cache
	if input_cache is None:
		# default is to use all output files from parents
		for parent in parents:
			node.add_input_cache(parent.get_output_cache())
	else:
		# but calling code can provide its own collection
		node.add_input_cache(input_cache)
	if extra_input_cache is not None:
		# sometimes it helps to add some extra
		node.add_input_cache(extra_input_cache)
	if preserve_cache is not None:
		node.add_preserve_cache(preserve_cache)

	# construct names for the node and output file, and override the
	# segment if needed
	[cache_entry] = node.get_output_cache()
	if segment is None:
		segment = cache_entry.segment
	node.set_name("lladd_%s_%d_%d" % (tag, int(segment[0]), int(abs(segment))))
	node.set_output(os.path.join(node.output_dir, "%s-%s-%d-%d.xml.gz" % (cache_entry.observatory, tag, int(segment[0]), int(abs(segment)))), segment = segment)

	node.set_retry(3)
	dag.add_node(node)
	return set([node])


def make_power_fragment(dag, parents, instrument, seg, tag, framecache, injargs = {}):
	node = PowerNode(powerjob)
	node.set_name("lalapps_power_%s_%s_%d_%d" % (tag, instrument, int(seg[0]), int(abs(seg))))
	map(node.add_parent, parents)
	# FIXME:  PowerNode should not be subclassed from AnalysisNode,
	# because that class is too hard-coded.  For example, there is no
	# way to switch to analysing gaussian noise except to comment out
	# this line in the code.
	node.set_cache(framecache)
	node.set_ifo(instrument)
	node.set_start(seg[0])
	node.set_end(seg[1])
	node.set_user_tag(tag)
	for arg, value in injargs.iteritems():
		# this is a hack, but I can't be bothered
		node.add_var_arg("--%s %s" % (arg, value))
	dag.add_node(node)
	return set([node])


def make_binj_fragment(dag, seg, time_slides_cache_entry, tag, offset, flow = None, fhigh = None):
	# adjust start time to be commensurate with injection period
	start = seg[0] - seg[0] % binjjob.time_step + binjjob.time_step * offset

	node = BurstInjNode(binjjob)
	node.set_time_slide_file(time_slides_cache_entry.path)
	node.set_start(start)
	node.set_end(seg[1])
	if flow is not None:
		node.set_name("lalapps_binj_%s_%d_%d" % (tag, int(start), int(flow)))
	else:
		node.set_name("lalapps_binj_%s_%d" % (tag, int(start)))
	node.set_user_tag(tag)
	if flow is not None:
		node.add_macro("macroflow", flow)
	if fhigh is not None:
		node.add_macro("macrofhigh", fhigh)
	node.add_macro("macroseed", int(time.time()%100 + start))
	dag.add_node(node)
	return set([node])


def make_binjfind_fragment(dag, parents, tag, verbose = False):
	input_cache = collect_output_caches(parents)
	nodes = set()
	while input_cache:
		node = BinjfindNode(binjfindjob)
		node.add_input_cache([cache_entry for (cache_entry, parent) in input_cache[:binjfindjob.files_per_binjfind]])
		for parent in set(parent for cache_entry, parent in input_cache[:binjfindjob.files_per_binjfind]):
			node.add_parent(parent)
		del input_cache[:binjfindjob.files_per_binjfind]
		seg = cache_span(node.get_input_cache())
		node.set_name("lalburst_injfind_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
		node.add_macro("macrocomment", tag)
		dag.add_node(node)
		nodes.add(node)
	return nodes


def make_bucluster_fragment(dag, parents, tag, verbose = False):
	input_cache = collect_output_caches(parents)
	nodes = set()
	while input_cache:
		node = BuclusterNode(buclusterjob)
		node.add_input_cache([cache_entry for (cache_entry, parent) in input_cache[:buclusterjob.files_per_bucluster]])
		for parent in set(parent for cache_entry, parent in input_cache[:buclusterjob.files_per_bucluster]):
			node.add_parent(parent)
		del input_cache[:buclusterjob.files_per_bucluster]
		seg = cache_span(node.get_input_cache())
		node.set_name("lalburst_cluster_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
		node.add_macro("macrocomment", tag)
		node.set_retry(3)
		dag.add_node(node)
		nodes.add(node)
	return nodes


def make_bucut_fragment(dag, parents, tag, verbose = False):
	input_cache = collect_output_caches(parents)
	nodes = set()
	while input_cache:
		node = BucutNode(bucutjob)
		node.add_input_cache([cache_entry for (cache_entry, parent) in input_cache[:bucutjob.files_per_bucut]])
		for parent in set(parent for cache_entry, parent in input_cache[:bucutjob.files_per_bucut]):
			node.add_parent(parent)
		del input_cache[:bucutjob.files_per_bucut]
		seg = cache_span(node.get_input_cache())
		node.set_name("lalburst_cut_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
		node.add_macro("macrocomment", tag)
		dag.add_node(node)
		nodes.add(node)
	return nodes


def make_burca_fragment(dag, parents, tag, coincidence_segments = None, verbose = False):
	input_cache = collect_output_caches(parents)
	if coincidence_segments is not None:
		# doesn't sense to supply this keyword argument for
		# more than one input file
		assert len(input_cache) == 1
	nodes = set()
	while input_cache:
		node = BurcaNode(burcajob)
		node.add_input_cache([cache_entry for (cache_entry, parent) in input_cache[:burcajob.files_per_burca]])
		for parent in set(parent for cache_entry, parent in input_cache[:burcajob.files_per_burca]):
			node.add_parent(parent)
		del input_cache[:burcajob.files_per_burca]
		seg = cache_span(node.get_input_cache())
		node.set_name("lalburst_coinc_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
		if coincidence_segments is not None:
			node.set_coincidence_segments(coincidence_segments)
		node.add_macro("macrocomment", tag)
		dag.add_node(node)
		nodes.add(node)
	return nodes


def make_sqlite_fragment(dag, parents, tag, verbose = False):
	input_cache = collect_output_caches(parents)
	nodes = set()
	for cache_entry, parent in input_cache:
		node = SQLiteNode(sqlitejob)
		node.add_input_cache([cache_entry])
		node.add_parent(parent)
		node.set_name("ligolw_sqlite_%s_%d" % (tag, len(nodes)))
		node.set_output(cache_entry.path.replace(".xml.gz", ".sqlite"))
		dag.add_node(node)
		nodes.add(node)
	return nodes


def make_burca_tailor_fragment(dag, input_cache, seg, tag):
	input_cache = list(input_cache)
	input_cache.sort(reverse = True)
	nodes = set()
	max_cost_per_job = 25	# 10000 s -equivalent files
	while input_cache:
		cache = []
		cost = 0
		while input_cache and cost <= max_cost_per_job:
			cache.append(input_cache.pop())
			# cost porportional to segment duration squared
			cost += (float(abs(cache[-1].segment)) / 10000.0)**2
		node = BurcaTailorNode(burcatailorjob)
		node.add_input_cache(cache)
		node.set_name("lalburst_power_meas_likelihood_%s_%d_%d_%d" % (tag, int(seg[0]), int(abs(seg)), len(nodes)))
		node.set_output("%s_%d" % (tag, len(nodes)))
		dag.add_node(node)
		nodes.add(node)
	node = BurcaTailorNode(burcatailorjob)
	node.set_name("lalburst_power_meas_likelihood_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
	for parent in nodes:
		node.add_parent(parent)
		node.add_input_cache(parent.get_output_cache())
	del node.get_args()[:]
	node.add_var_arg("--add-from-cache %s" % node.cache_name)
	node.set_output(tag)
	dag.add_node(node)
	delete_cache = set(node.get_input_cache()) - set(node.get_output_cache())
	if delete_cache:
		rmnode = RMNode(rmjob)
		rmnode.set_name("lalburst_power_meas_likelihood_rm_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
		rmnode.add_parent(node)
		rmnode.add_input_cache(delete_cache)
		dag.add_node(rmnode)
	return set([node])


def make_burca2_fragment(dag, coinc_cache, likelihood_parents, tag):
	# FIXME:  pass a node set instead of a cache
	#input_cache = collect_output_caches(coinc_parents)
	coinc_cache = list(coinc_cache)
	coinc_cache.sort(reverse = True)

	likelihood_data_cache_filename = os.path.join(burca2job.cache_dir, "burca2_%s.cache" % tag)
	likelihood_data_cache_file = file(likelihood_data_cache_filename, "w")
	for cache_entry in [cache_entry for node in likelihood_parents for cache_entry in node.get_output_cache()]:
		print(str(cache_entry), file=likelihood_data_cache_file)

	nodes = set()
	max_cost_per_job = 10	# 10000 s -equivalent files
	while coinc_cache:
		cache = []
		cost = 0
		while coinc_cache and cost <= max_cost_per_job:
			cache.append(coinc_cache.pop())
			# cost porportional to segment duration squared
			cost += (float(abs(cache[-1].segment)) / 10000.0)**2
		node = BurcaNode(burca2job)
		node.set_name("lalburst_coinc2_%s_%d" % (tag, len(nodes)))
		node.add_macro("macrocomment", tag)
		node.add_var_arg("--likelihood-data-cache %s" % likelihood_data_cache_filename)
		node.add_input_cache(cache)
		for parent in likelihood_parents:
			node.add_parent(parent)
		dag.add_node(node)
		nodes.add(node)
	return nodes


#
# =============================================================================
#
#                             ligo_data_find Stage
#
# =============================================================================
#


def make_datafind_stage(dag, seglists, verbose = False):
	if verbose:
		print("building ligo_data_find jobs ...", file=sys.stderr)

	#
	# Fill gaps smaller than the padding added to each datafind job.
	# Filling in the gaps ensures that exactly 1 datafind job is
	# suitable for each lalapps_power job, and also hugely reduces the
	# number of ligo_data_find nodes in the DAG.
	#

	filled = seglists.copy().protract(datafind_pad / 2).contract(datafind_pad / 2)

	#
	# Build the nodes.  Do this in time order to assist depth-first job
	# submission on clusters.
	#

	segs = [(seg, instrument) for instrument, seglist in filled.iteritems() for seg in seglist]
	segs.sort()

	nodes = set()
	for seg, instrument in segs:
		if verbose:
			print("making datafind job for %s spanning %s" % (instrument, seg), file=sys.stderr)
		new_nodes = make_datafind_fragment(dag, instrument, seg)
		nodes |= new_nodes

		# add a post script to check the file list
		#required_segs_string = ",".join(segmentsUtils.to_range_strings(seglists[instrument] & segments.segmentlist([seg])))
		#for node in new_nodes:
		#	node.set_post_script(datafindjob.get_config_file().get("condor", "LSCdataFindcheck") + " --dagman-return $RETURN --stat --gps-segment-list %s %s" % (required_segs_string, node.get_output()))

	return nodes


#
# =============================================================================
#
#        Analyze All Segments in a segmentlistdict Using lalapps_power
#
# =============================================================================
#


#
# one segment
#


def make_power_segment_fragment(dag, datafindnodes, instrument, segment, tag, timing_params, psds_per_job, binjnodes = set(), verbose = False):
	"""
	Construct a DAG fragment for an entire segment, splitting the
	segment into multiple trigger generator jobs.
	"""
	# only one frame cache file can be provided as input, and only one
	# injection description file can be provided as input
	# the unpacking indirectly tests that the file count is correct
	[framecache] = [node.get_output() for node in datafindnodes]
	if binjnodes:
		[simfile] = [cache_entry.path for node in binjnodes for cache_entry in node.get_output_cache()]
		injargs = {"injection-file": simfile}
	else:
		injargs = {}
	seglist = split_segment(timing_params, segment, psds_per_job)
	if verbose:
		print("Segment split: " + str(seglist), file=sys.stderr)
	nodes = set()
	for seg in seglist:
		nodes |= make_power_fragment(dag, datafindnodes | binjnodes, instrument, seg, tag, framecache, injargs = injargs)
	return nodes


#
# all segments
#


def make_single_instrument_stage(dag, datafinds, seglistdict, tag, timing_params, psds_per_job, binjnodes = set(), verbose = False):
	nodes = []
	for instrument, seglist in seglistdict.iteritems():
		for seg in seglist:
			if verbose:
				print("generating %s fragment %s" % (instrument, str(seg)), file=sys.stderr)

			# find the datafind job this job is going to need
			dfnodes = set([node for node in datafinds if (node.get_ifo() == instrument) and (seg in segments.segment(node.get_start(), node.get_end()))])
			if len(dfnodes) != 1:
				raise ValueError("error, not exactly 1 datafind is suitable for trigger generator job at %s in %s" % (str(seg), instrument))

			# trigger generator jobs
			nodes += make_power_segment_fragment(dag, dfnodes, instrument, seg, tag, timing_params, psds_per_job, binjnodes = binjnodes, verbose = verbose)

	# done
	return nodes


#
# =============================================================================
#
#                         Coincidence Post-Processing
#
# =============================================================================
#


def group_coinc_parents(parents, offset_vectors, extentlimit = None, verbose = False):
	if not offset_vectors:
		# no-op
		return []

	if verbose:
		print("Grouping jobs for coincidence analysis:", file=sys.stderr)

	#
	# use ligolw_cafe to group each output file according to how they
	# need to be combined to perform the coincidence analysis
	#

	seglists, bins = cafe.ligolw_cafe([cache_entry for parent in parents for cache_entry in parent.get_output_cache()], offset_vectors, extentlimit = extentlimit, verbose = verbose)

	#
	# retrieve the file caches and segments.  note that ligolw_cafe
	# returns the bins sorted by segment, so we do too
	#

	caches = [frozenset(bin.objects) for bin in bins]
	assert len(set(caches)) == len(caches)
	segs = [cache_span(bin.objects) for bin in bins]

	#
	# determine the clipping boundaries to use for each coincidence job
	# if an extentlimit has been imposed
	#

	clipsegs = [None] * len(bins)
	if extentlimit is not None:
		extents = [bin.extent for bin in bins]
		for i, extent in enumerate(extents):
			# FIXME:  when we can rely on Python >= 2.5,
			#lo = segments.NegInfinity if i == 0 or extents[i - 1].disjoint(extent) else extent[0]
			# etc.
			if i == 0 or extents[i - 1].disjoint(extent):
				lo = segments.NegInfinity
			else:
				lo = extent[0]
			if i >= len(extents) - 1 or extents[i + 1].disjoint(extent):
				hi = segments.PosInfinity
			else:
				hi = extent[1]
			if lo is not segments.NegInfinity or hi is not segments.PosInfinity:
				clipsegs[i] = segments.segment(lo, hi)

	#
	# match parents to caches
	#

	if verbose:
		print("Matching jobs to caches ...", file=sys.stderr)
	parent_groups, unused = match_nodes_to_caches(parents, caches)
	if verbose and unused:
		# there were parents that didn't match any caches.  this
		# happens if ligolw_cafe decides their outputs aren't
		# needed
		print("Notice:  %d jobs (of %d) produce output that will not be used by a coincidence job" % (unused, len(parents)), file=sys.stderr)

	#
	# done
	#

	return zip(segs, parent_groups, caches, clipsegs)
