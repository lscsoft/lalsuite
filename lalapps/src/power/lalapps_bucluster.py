#
# Copyright (C) 2006,2012  Kipp Cannon
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


from optparse import OptionParser


from lal.utils import CacheEntry


from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import process as ligolw_process
from lalburst import git_version
from lalburst import ligolw_bucluster


lsctables.use_in(ligolw.LIGOLWContentHandler)


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#


def parse_command_line():
	parser = OptionParser(
		version = "Name: %%prog\n%s" % git_version.verbose_msg,
		usage = "%prog [options] [file ...]",
		description = "Run a single-instrument burst clustering algorithm on the sngl_burst events contained in LIGO Light Weight XML files.  Files can be listed on the command line and/or in one or more LAL cache files.  If no files are named, then input is read from stdin and written to stdout."
	)
	parser.add_option("--comment", metavar = "text", help = "Set the comment string to be recorded in the process table (default = None).")
	parser.add_option("-c", "--cluster-algorithm", metavar = "[excesspower]", help = "Set clustering method (required).")
	parser.add_option("-i", "--input-cache", metavar = "filename", action = "append", default = [], help = "Process the files listed in this LAL cache.")
	parser.add_option("-p", "--program", metavar = "name", default = "lalapps_power", help = "Set the name of the program that generated the events as it appears in the process table (default = \"lalapps_power\").")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.cluster_algorithm is None:
		raise ValueError("missing required argument --cluster-algorithm")
	if options.cluster_algorithm not in ("excesspower",):
		raise ValueError("unrecognized --cluster-algorithm %s" % options.cluster_algorithm)

	for cache in options.input_cache:
		filenames += [CacheEntry(line).path for line in file(cache)]

	return options, (filenames or [None])


#
# =============================================================================
#
#                                    Input
#
# =============================================================================
#


#
# Use interning row builder to save memory.
#


lsctables.table.RowBuilder = lsctables.table.InterningRowBuilder


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


options, filenames = parse_command_line()


prefunc = {
	"excesspower": ligolw_bucluster.ExcessPowerPreFunc
}[options.cluster_algorithm]
postfunc = {
	"excesspower": ligolw_bucluster.ExcessPowerPostFunc
}[options.cluster_algorithm]
testfunc = {
	"excesspower": ligolw_bucluster.ExcessPowerTestFunc
}[options.cluster_algorithm]
sortfunc = {
	"excesspower": ligolw_bucluster.ExcessPowerSortFunc
}[options.cluster_algorithm]
bailoutfunc = {
	"excesspower": ligolw_bucluster.ExcessPowerBailoutFunc
}[options.cluster_algorithm]
clusterfunc = {
	"excesspower": ligolw_bucluster.ExcessPowerClusterFunc
}[options.cluster_algorithm]


for filename in filenames:
	#
	# Load document
	#

	xmldoc = utils.load_filename(filename, verbose = options.verbose, contenthandler = ligolw.LIGOLWContentHandler)
	lsctables.table.InterningRowBuilder.strings.clear()

	# FIXME:  don't do this:  fix lalapps_power's output
	if options.cluster_algorithm in ("excesspower",):
		ligolw_bucluster.add_ms_columns(xmldoc)

	#
	# Add process information
	#

	process = ligolw_bucluster.append_process(xmldoc, cluster_algorithm = options.cluster_algorithm, comment = options.comment)

	#
	# Call clustering library
	#

	xmldoc, changed = ligolw_bucluster.ligolw_bucluster(
		xmldoc,
		program = options.program,
		process = process,
		prefunc = prefunc,
		postfunc = postfunc,
		testfunc = testfunc,
		clusterfunc = clusterfunc,
		sortfunc = sortfunc,
		bailoutfunc = bailoutfunc,
		verbose = options.verbose
	)

	#
	# Finish process information
	#

	ligolw_process.set_process_end_time(process)

	#
	# Write document
	#

	if changed:
		utils.write_filename(xmldoc, filename, gz = (filename or "stdout").endswith(".gz"), verbose = options.verbose)
	xmldoc.unlink()
