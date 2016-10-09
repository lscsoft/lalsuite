#
# Copyright (C) 2010  Kipp Cannon
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
import sqlite3
import string
import sys


from lal.utils import CacheEntry


from glue import segments
from glue.ligolw import dbtables
from glue.ligolw import utils
from glue.ligolw.utils import process as ligolwprocess
from lalburst import SnglBurstUtils
from lalburst import git_version
from lalburst import ligolw_burca_tailor
from lalburst import stringutils


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


# characters allowed to appear in the description string
T010150_letters = set(string.ascii_lowercase + string.ascii_uppercase + string.digits + "_+#")


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
		usage = "%prog [options] [filename ...]",
		description = "%prog analyzes a collection of sqlite3 database files containing ligolw_burca outputs of string-cusp coincidence events, and measures probability distributions for a variety of parameters computed from those coincidences.  The distributions are written to a likelihood data file in XML format, which can later be used by to assign likelihoods to the coincidences.  The files to be processed can be named on the command line and/or provided by a LAL cache file."
	)
	parser.add_option("-o", "--output", metavar = "filename", default = None, help = "Set the name of the likelihood data file to write (default = stdout).")
	parser.add_option("-c", "--input-cache", metavar = "filename", help = "Also process the files named in this LAL cache.  See lalapps_path2cache for information on how to produce a LAL cache file.")
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("--T010150", metavar = "description", default = None, help = "Write the output to a file whose name is compatible with the file name format described in LIGO-T010150-00-E, \"Naming Convention for Frame Files which are to be Processed by LDAS\".  The description string will be used to form the second field in the file name.")
	parser.add_option("--injection-reweight", metavar = "off|astrophysical", default = "off", help = "Set the weight function to be applied to the injections (default = \"off\").  When \"off\", the injections are all given equal weight and so the injection population is whatever was injected.  When set to \"astrophysical\", the injections are reweighted to simulate an amplitude^{-4} distribution.")
	parser.add_option("--injection-reweight-cutoff", metavar = "amplitude", default = 1e-20, type = "float", help = "When using the astrophysical injection reweighting, do not allow the weight assigned to arbitrarily low-amplitude injections to grow without bound, instead clip the weight assigned to injections to the weight given to injections with this amplitude (default = 1e-20, 0 = disabled).  This option is ignored when astrophysical reweighting is not being performed.")
	parser.add_option("--vetoes-name", metavar = "name", help = "Set the name of the segment lists to use as vetoes (default = do not apply vetoes).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.T010150 is not None:
		if options.output is not None:
			raise ValueError("cannot set both --T010150 and --output")
		if options.T010150 == "":
			options.T010150 = "STRING_LIKELIHOOD"
		elif set(options.T010150) - T010150_letters:
			raise ValueError("invalid characters in description \"%s\"" % options.T010150)

	if options.input_cache:
		filenames += [CacheEntry(line).path for line in file(options.input_cache)]

	if not filenames:
		raise ValueError("no input files!")

	if options.injection_reweight not in ("off", "astrophysical"):
		raise ValueError("--injection-reweight \"%s\" not recognized" % options.injections_reweight)

	return options, filenames


#
# =============================================================================
#
#                         Injection ReWeight Functions
#
# =============================================================================
#


def get_injection_weight_func(contents, reweight_type, amplitude_cutoff):
	if reweight_type == "off":
		def weight_func(sim, amplitude_cutoff = amplitude_cutoff):
			# amplitude cut-off is not used
			return 1.0
	elif reweight_type == "astrophysical":
		population, = ligolwprocess.get_process_params(contents.xmldoc, "lalapps_binj", "--population")
		if population != "string_cusp":
			raise ValueError("lalapps_binj was not run with --population=\"string_cusp\"")
		def weight_func(sim, amplitude_cutoff = amplitude_cutoff):
			# the "string_cusp" injection population is uniform
			# in log A, meaning the number of injections
			# between log A and log A + d log A is independent
			# of A.  that corresponds to a distribution density
			# in A of P(A) dA \propto A^{-1} dA.  we want P(A)
			# dA \propto A^{-4} dA, so each physical injection
			# needs to be treated as if it was A^{-3} virtual
			# injections, then the density of virtual
			# injections will be P(A) dA \propto A^{-4} dA.
			# the factor of 10^{21} simply renders the
			# amplitudes closer to 1;  since double-precision
			# arithmetic is used this should have no affect on
			# the results, but it might help make the numbers a
			# little easier for humans to look at should anyone
			# have occasion to do so.
			return (max(sim.amplitude, amplitude_cutoff) * 1e21)**-3
	else:
		raise ValueError(reweight_type)
	return weight_func


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


#
# Command line.
#


options, filenames = parse_command_line()


#
# Clear the statistics book-keeping object.
#


distributions = stringutils.StringCoincParamsDistributions()
segs = segments.segmentlistdict()


#
# Iterate over files
#


for n, filename in enumerate(filenames):
	#
	# Open the database file.
	#

	if options.verbose:
		print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)

	working_filename = dbtables.get_connection_filename(filename, tmp_path = options.tmp_space, verbose = options.verbose)
	connection = sqlite3.connect(working_filename)
	if options.tmp_space is not None:
		dbtables.set_temp_store_directory(connection, options.tmp_space, verbose = options.verbose)

	#
	# Summarize the database.
	#

	contents = SnglBurstUtils.CoincDatabase(connection, live_time_program = "StringSearch", search = "StringCusp", veto_segments_name = options.vetoes_name)
	if options.verbose:
		SnglBurstUtils.summarize_coinc_database(contents)
	if not contents.seglists and options.verbose:
		print >>sys.stderr, "\twarning:  no segments found"
	segs |= contents.seglists

	#
	# Build triangulators.  The timing uncertainty of +/- 8e-5 s was
	# measured with lalapps_string_plot_binj and is essentially
	# identical for H1, H2, L1, and V1.
	#

	triangulators = stringutils.triangulators(dict.fromkeys(contents.instruments, 8e-5))

	#
	# Record statistics.  Assume all files with sim_burst tables are
	# the outputs of injection runs, and others aren't.
	#

	if contents.sim_burst_table is None:
		# iterate over burst<-->burst coincs
		for is_background, events, offsetvector in ligolw_burca_tailor.get_noninjections(contents):
			params = distributions.coinc_params([event for event in events if event.ifo not in contents.vetoseglists or event.peak not in contents.vetoseglists[event.ifo]], offsetvector, triangulators)
			if params is not None:
				if is_background:
					distributions.add_background(params)
				else:
					distributions.add_zero_lag(params)
	else:
		weight_func = get_injection_weight_func(contents, options.injection_reweight, options.injection_reweight_cutoff)
		# iterate over burst<-->burst coincs matching injections
		# "exactly"
		for sim, events, offsetvector in ligolw_burca_tailor.get_injections(contents):
			params = distributions.coinc_params([event for event in events if event.ifo not in contents.vetoseglists or event.peak not in contents.vetoseglists[event.ifo]], offsetvector, triangulators)
			if params is not None:
				distributions.add_injection(params, weight = weight_func(sim))

	#
	# Clean up.
	#

	contents.xmldoc.unlink()
	connection.close()
	dbtables.discard_connection_filename(filename, working_filename, verbose = options.verbose)


#
# Output.
#


def T010150_basename(description, seglists):
	seg = seglists.extent_all()
	return "%s-%s-%s-%s" % ("+".join(sorted(seglists.keys())), description, str(seg[0]), str(abs(seg)))

if options.T010150:
	filename = "%s.xml.gz" % T010150_basename(options.T010150, segs)
else:
	filename = options.output

stringutils.write_likelihood_data(filename, distributions, segs, verbose = options.verbose)
