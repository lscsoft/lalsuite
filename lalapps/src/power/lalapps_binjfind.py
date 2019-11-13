#
# Copyright (C) 2006  Kipp Cannon
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
Command-line interface to burst injection identification code.
"""


from __future__ import print_function


from optparse import OptionParser
import sys


from ligo.lw import ligolw
from ligo.lw import lsctables
from ligo.lw import utils as ligolw_utils
from ligo.lw.utils import process as ligolw_process
from lalburst import git_version
from lalburst import binjfind


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
		description = "Accepts as input one or more LIGO Light Weight XML files, each containing burst candidates and a list of injections, and adds entries to the coincidence tables indicating which burst events match which injections."
	)
	parser.add_option("-f", "--force", action = "store_true", help = "Process even if file has already been processed.")
	parser.add_option("--comment", metavar = "text", help = "Set the comment string to be written to the process table (default = None).")
	parser.add_option("-c", "--match-algorithm", metavar = "[stringcusp|excesspower|omega|waveburst]", default = None, help = "Set the algorithm used to match burst candidates with injections (required).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.match_algorithm is None:
		raise ValueError("missing required --match-algorithm option")
	if options.match_algorithm not in ("stringcusp", "excesspower", "omega", "waveburst"):
		raise ValueError("unrecognized match algorithm \"%s\"" % options.match_algorithm)

	return options, (filenames or [None])


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


#
# command line
#


options, filenames = parse_command_line()

# must match columns in sngl_burst table
search = {
	"stringcusp": "StringCusp",
	"excesspower": "excesspower",
	"omega": "omega",
	"waveburst": "waveburst"
}[options.match_algorithm]

snglcomparefunc = {
	"stringcusp": binjfind.StringCuspSnglCompare,
	"excesspower": binjfind.ExcessPowerSnglCompare,
	"omega": binjfind.OmegaSnglCompare,
	"waveburst": binjfind.CWBSnglCompare
}[options.match_algorithm]

nearcoinccomparefunc = {
	"stringcusp": binjfind.StringCuspNearCoincCompare,
	"excesspower": binjfind.ExcessPowerNearCoincCompare,
	"omega": binjfind.OmegaNearCoincCompare,
	"waveburst": binjfind.CWBNearCoincCompare
}[options.match_algorithm]


#
# loop over files
#


for n, filename in enumerate(filenames):
	#
	# load the document
	#

	if options.verbose:
		print("%d/%d:" % (n + 1, len(filenames)), end=' ', file=sys.stderr)
	xmldoc = ligolw_utils.load_filename(filename, verbose = options.verbose, contenthandler = ligolw.LIGOLWContentHandler)

	#
	# have we already procesed it?
	#

	if ligolw_process.doc_includes_process(xmldoc, binjfind.process_program_name):
		if options.verbose:
			print("warning: %s already processed," % (filename or "stdin"), end=' ', file=sys.stderr)
		if not options.force:
			if options.verbose:
				print("skipping (use --force to force)", file=sys.stderr)
			continue
		if options.verbose:
			print("continuing by --force", file=sys.stderr)

	#
	# add process metadata to document
	#

	process = binjfind.append_process(xmldoc, match_algorithm = options.match_algorithm, comment = options.comment)

	#
	# run binjfind algorithm
	#

	binjfind.binjfind(xmldoc, process, search, snglcomparefunc, nearcoinccomparefunc, verbose = options.verbose)

	#
	# close out the process metadata
	#

	ligolw_process.set_process_end_time(process)

	#
	# done
	#

	ligolw_utils.write_filename(xmldoc, filename, verbose = options.verbose, gz = (filename or "stdout").endswith(".gz"))
	xmldoc.unlink()
	lsctables.reset_next_ids(lsctables.TableByName.values())
