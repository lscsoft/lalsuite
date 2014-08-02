#
# Copyright (C) 2005  Kipp C. Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Pipeline generation script for the excess power online analysis.
"""


import ConfigParser
from optparse import OptionParser
import os
import sys
import tempfile


from glue import pipeline
from glue import segments
from glue import segmentsUtils
from pylal.date import LIGOTimeGPS
from lalapps import power


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"


#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#


def parse_command_line():
	parser = OptionParser(
		version="%prog CVS $Id$"
	)
	parser.add_option("-s", "--data-start", metavar = "GPSSECONDS", help = "set data segment start time")
	parser.add_option("-e", "--data-end", metavar = "GPSSECONDS", help = "set data segment end time")
	parser.add_option("-a", "--trig-start", metavar = "GPSSECONDS", help = "set analysis segment start time")
	parser.add_option("-b", "--trig-end", metavar = "GPSSECONDS", help = "set analysis segment end time")
	parser.add_option("-f", "--dag-name", metavar = "FILENAME", help = "set output .dag file name")
	parser.add_option("-t", "--aux-dir", metavar = "PATH", help = "set auxiliary data directory")
	parser.add_option("--condor-log-dir", metavar = "PATH", help = "set directory for Condor log")
	parser.add_option("--config-file", metavar = "FILENAME", default = "online_power.ini", help = "set .ini config file name")
	parser.add_option("--instrument", metavar = "INSTRUMENT", help = "set instrument name (default = value of instrument variable in [pipeline] section of .ini file)")
	parser.add_option("--publish-dest", metavar = "PATH", help = "set directory for output triggers")
	parser.add_option("--dmtmon-dest", metavar = "PATH", help = "set directory for DMT monitor output")
	parser.add_option("--gsiscp-dest", metavar = "PATH", help = "set destination for gsiscp")
	parser.add_option("--user-tag", metavar = "TAG", help = "set user tag on jobs that need it")

	options, extra_args = parser.parse_args()

	# data segment
	options.data_seg = segments.segment(LIGOTimeGPS(options.data_start), LIGOTimeGPS(options.data_end))
	try:
		options.data_seglists = segments.segmentlistdict({options.instrument: segments.segmentlist([options.data_seg])})
	except:
		raise ValueError, "failure parsing -s and/or -e; try --help"

	# trigger segment
	options.trig_seg = segments.segment(LIGOTimeGPS(options.trig_start), LIGOTimeGPS(options.trig_end))
	try:
		options.trig_seglists = segments.segmentlistdict({options.instrument: segments.segmentlist([options.trig_seg])})
	except:
		raise ValueError, "failure parsing -a and/or -b; try --help"
	if True in map(bool, (options.trig_seglists - options.data_seglists).itervalues()):
		raise ValueError, "trigger segment not contained in data segment!"

	# .dag name
	try:
		options.dag_name = os.path.splitext(options.dag_name)[0]
	except:
		raise ValueError, "failure parsing -f; try --help"

	# config file
	options.config_file = os.path.join(options.aux_dir, options.config_file)

	return options


options = parse_command_line()


#
# =============================================================================
#
#                                 Config File
#
# =============================================================================
#


# .ini file config parser
config_parser = ConfigParser.SafeConfigParser()
config_parser.read(options.config_file)


# Condor log file
condor_log = tempfile.mkstemp(".log", "power_", options.condor_log_dir)[1]


# instrument
if options.instrument:
	config_parser.set('pipeline', 'instrument', options.instrument)
else:
	options.instrument = config_parser.get('pipeline', 'instrument')


# publish location
if options.publish_dest:
	config_parser.set('publish', 'destination', options.publish_dest)
else:
	options.publish_dest = config_parser.get('publish', 'destination')


# gsiscp location
if options.gsiscp_dest:
	config_parser.set('gsiscp', 'destination', options.gsiscp_dest)
else:
	options.gsiscp_dest = config_parser.get('gsiscp', 'destination')


# user tag
if options.user_tag:
	config_parser.set('pipeline', 'user_tag', options.user_tag)
else:
	options.user_tag = config_parser.get('pipeline', 'user_tag')


# timing params
timing_params = power.get_timing_parameters(config_parser)


#
# =============================================================================
#
#                              Static Parameters
#
# =============================================================================
#


if abs(options.data_seg) < 3600:
	psds_per_job = 4
	psds_per_injection = 4
else:
	psds_per_job = 16
	psds_per_injection = 16


#
# =============================================================================
#
#                               Custom Job Types
#
# =============================================================================
#


class Burst2MonJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", power.get_executable(config_parser, "ligolw_burst2mon"))
		self.set_sub_file("ligolw_burst2mon.sub")
		self.set_stdout_file(os.path.join(power.get_out_dir(config_parser), "ligolw_burst2mon-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(power.get_out_dir(config_parser), "ligolw_burst2mon-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_burst2mon")


class Burst2MonNode(pipeline.AnalysisNode):
	def __init__(self, job):
		# FIXME: this shouldn't be needed.
		pipeline.CondorDAGNode.__init__(self, job)
		pipeline.AnalysisNode.__init__(self)


class PublishJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", power.get_executable(config_parser, "publish"))
		self.set_sub_file("publish.sub")
		self.set_stdout_file(os.path.join(power.get_out_dir(config_parser), "publish-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(power.get_out_dir(config_parser), "publish-$(cluster)-$(process).err"))

	def write_sub_file(self):
		# this is a hack, but I can't be bothered...
		self.add_file_arg("$(macrodestination)")
		pipeline.CondorDAGJob.write_sub_file(self)


class PublishNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)
		for c in cache:
			filename = c.path
			pipeline.CondorDAGNode.add_file_arg(self, filename)
	
	def set_output(self, destination):
		self.add_macro("macrodestination", destination)


class GsiScpJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "local", power.get_executable(config_parser, "gsiscp"))
		self.set_sub_file("gsiscp.sub")
		self.set_stdout_file(os.path.join(power.get_out_dir(config_parser), "gsiscp-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(power.get_out_dir(config_parser), "gsiscp-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")

	def write_sub_file(self):
		# this is a hack, but I can't be bothered...
		self.add_file_arg("$(macrodestination)")
		pipeline.CondorDAGJob.write_sub_file(self)


class GsiScpNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)
		for c in cache:
			filename = c.path
			pipeline.CondorDAGNode.add_file_arg(self, filename)
	
	def set_output(self, destination):
		self.add_macro("macrodestination", destination)


#
# =============================================================================
#
#                              Define .sub Files
#
# =============================================================================
#


power.init_job_types(config_parser, types = ["datafind", "binj", "power", "lladd", "binjfind", "bucluster", "bucut"])


llb2mjob = Burst2MonJob(config_parser)
publishjob = PublishJob(config_parser)
gsiscpjob = GsiScpJob(config_parser)


#
# =============================================================================
#
#                             Publish DAG Fragment
#
# =============================================================================
#


def make_publish_fragment(dag, parents, segment, tag, destination):
	start = int(segment[0])
	destination = os.path.join(destination, "%d/%d" % (start/1000000, start/100000))
	try:
		os.makedirs(destination)
	except OSError, e:
		# errno 17 == "file exists"
		if e.errno != 17:
			raise e
	node = PublishNode(publishjob)
	node.set_name("publish_%s" % tag)
	for parent in parents:
		node.add_parent(parent)
		node.add_input_cache(parent.get_output_cache())
	node.set_output(destination)
	dag.add_node(node)

	return [node]


#
# =============================================================================
#
#                       DMT Monitor Output DAG Fragment
#
# =============================================================================
#


# FIXME: still broken
def make_burst2mon_fragment(dag, parent, instrument, seg, tag):
	cluster_output = "%s-POWERMON_%s-%s-%s.xml" % (instrument, tag, int(seg[0]), int(abs(seg)))
	cluster = power.make_bucluster_fragment(dag, [], instrument, seg, "POWERMON_%s" % tag)
	cluster.add_parent(parent)
	cluster.set_pre_script("/bin/cp %s %s" % (parent.get_output_files()[0], cluster_output))
	cluster.add_file_arg(cluster_output)

	node = Burst2MonNode(llb2mjob)
	node.set_name("ligolw_burst2mon")
	node.add_parent(cluster)
	node.set_input(cluster.get_output_files()[0])
	node.set_output(os.path.join(options.dmtmon_dest, cluster.get_output_files()[0]))
	node.add_macro("macrocomment", tag)
	dag.add_node(node)

	return node


#
# =============================================================================
#
#                             GSISCP DAG Fragment
#
# =============================================================================
#


def make_gsiscp_fragment(dag, parents, destination):
	node = GsiScpNode(gsiscpjob)
	node.set_name("gsiscp")
	for parent in parents:
		node.add_parent(parent)
		node.add_input_cache(parent.get_output_cache())
	node.set_output(destination)
	dag.add_node(node)
	return [node]


#
# =============================================================================
#
#                               DAG Construction
#
# =============================================================================
#


power.make_dag_directories(config_parser)


dag = pipeline.CondorDAG(condor_log)
dag.set_dag_file(options.dag_name)


datafinds = power.make_datafind_stage(dag, options.data_seglists, verbose = True)


nodes = power.make_single_instrument_stage(dag, datafinds, options.data_seglists, options.user_tag, timing_params, psds_per_job, verbose = True)
nodes = power.make_lladded_bucluster_fragment(dag, nodes, options.data_seg, options.user_tag)


make_publish_fragment(dag, nodes, options.data_seg, options.user_tag, options.publish_dest)
# FIXME: still broken
if options.dmtmon_dest:
	make_burst2mon_fragment(dag, darmpowerfrag, options.instrument, options.data_seg, options.user_tag)


# FIXME: still broken
nodes = power.make_single_instrument_injections_stage(dag, datafinds, binjnodes, options.data_seglists, "INJECTIONS_%s" % options.user_tag, timing_params, psds_per_injection, verbose = True)
nodes = power.make_lladded_bucluster_fragment(dag, nodes, options.data_seg, "INJECTIONS_%s" % options.user_tag)
nodes = power.make_bucut_fragment(dag, nodes, "INJECTIONS_%s" % options.user_tag)
nodes = power.make_binjfind_fragment(dag, nodes, "INJECTIONS_%s" % options.user_tag)


make_publish_fragment(dag, nodes, options.data_seg, "INJECTIONS_%s" % options.user_tag, options.publish_dest)


# FIXME: still broken
make_gsiscp_fragment(dag, [darmpowerfrag, binjfindfrag], options.gsiscp_dest)


dag.write_sub_files()
dag.write_dag()


print str(options.trig_seg[0]), str(options.trig_seg[1])
