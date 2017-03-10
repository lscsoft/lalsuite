#
# Copyright (C) 2006--2014  Kipp Cannon
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
import sys


from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw.utils import time_slide as ligolw_time_slide
from glue import offsetvector
from lalapps import git_version
from lalburst import timeslides


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


def parse_normalize(normalize):
	return dict((name.strip(), float(offset)) for name, offset in map(lambda s: s.split("="), normalize))


def parse_command_line():
	parser = OptionParser(
		version = "Name: %%prog\n%s" % git_version.verbose_msg,
		usage = "%prog [options] [filename ...]",
		description = "%prog constructs time_slide tables, writing the result to one or more files.  The time slide table to be constructed is described by specifying one or more ranges of offsets for each instrument.  If more than one instrument and set of offsets is given, then the time slide table will contain entries corresponding to all combinations of offsets, one each from the different instuments.  If no file names are given on the command line, output is written to stdout.  If more than one file name is given on the command line, then the time slides are distributed uniformly between files with each file being given a disjoint subset of the time slides.  Output files whose names end in \".gz\" will be gzip compressed.\n\nExample:\n\n%prog --verbose --instrument H1=-100:+100:10 --instrument L1=-100:+100:+10 time_slides.xml.gz"
	)
	parser.add_option("-a", "--add-to", metavar = "filename", action = "append", default = [], help = "Add the time slides from this file to the newly-generated time slides.  If the name ends in \".gz\" it will be gzip-decompressed on input.")
	parser.add_option("--comment", metavar = "text", help = "Set comment string in process table (default = None).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	parser.add_option("-i", "--instrument", metavar = "name=first:last:step[,first:last:step[,...]]", action = "append", default = [], help = "Provide a description of the set of offsets to use for a particular instrument.  The set of offsets is (first + n * step) where n is an integer such that first <= offset <= last.  More than one set of offsets can be given for the same instrument, in which case the union is used.  As a short-hand, the sets can be combined into a single command line argument by separating the first:last:step triples with commas.")
	parser.add_option("--inspiral-num-slides", metavar = "count:instrument=offset[,instrument=offset...]", action = "append", default = [], help = "Generate a set of inspiral group style time slides.  The collection of instrument/offset pairs defines an offset vector, and the time slides produced are integer multiples of that vector n * {offsets} where n is a non-zero integer in [-count, +count] (so if count=50, you get 101 time slides).  If this option is given more than once, then multiple sets of inspiral-style time slides are generated.")
	parser.add_option("-n", "--normalize", metavar = "name=offset", default = [], action = "append", help = "Normalize the time slides so that this instrument has the specified offset in all.  The other offsets in each time slide are adjusted so that the relative offsets are preserved.  Time slides that do not involve this instrument are unaffected.  If this option is given multiple times, then for each time slide they are considered in alphabetical order by instrument until the first is found that affects the offsets of that time slide.")
	parser.add_option("--remove-zero-lag", action = "store_true", help = "Remove the time slide with offsets of 0 for all instruments.")
	options, filenames = parser.parse_args()

	try:
		parse_normalize(options.normalize)
	except Exception as e:
		raise ValueError("unable to parse --normalize arguments: %s" % str(e))

	return options, (filenames or [None])


#
# =============================================================================
#
#                                 Preparation
#
# =============================================================================
#


def new_doc(comment = None, **kwargs):
	doc = ligolw.Document()
	doc.appendChild(ligolw.LIGO_LW())
	process = ligolw_process.register_to_xmldoc(
		doc,
		program = u"lalapps_gen_timeslides",
		paramdict = kwargs,
		version = __version__,
		cvs_repository = u"lscsoft",
		cvs_entry_time = __date__,
		comment = comment
	)
	timeslidetable = lsctables.New(lsctables.TimeSlideTable)
	doc.childNodes[0].appendChild(timeslidetable)

	return doc, process


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
# Load initial time slides.
#


@lsctables.use_in
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	pass

time_slides = {}
for filename in options.add_to:
	time_slide_table = lsctables.TimeSlideTable.get_table(ligolw_utils.load_filename(filename, verbose = options.verbose, contenthandler = LIGOLWContentHandler))
	extra_time_slides = time_slide_table.as_dict().values()
	if options.verbose:
		print >>sys.stderr, "Loaded %d time slides." % len(extra_time_slides)
	for offsetvect in extra_time_slides:
		time_slides[lsctables.TimeSlideTable.get_next_id()] = offsetvect


#
# Make new time slides.
#


if options.verbose:
	print >>sys.stderr, "Computing new time slides ..."

# dictionary mapping time_slide_id --> (dictionary mapping insrument --> offset)

for offsetvect in timeslides.SlidesIter(timeslides.parse_slides(options.instrument)):
	time_slides[lsctables.TimeSlideTable.get_next_id()] = offsetvect
for inspiral_slidespec in options.inspiral_num_slides:
	for offsetvect in timeslides.Inspiral_Num_Slides_Iter(*timeslides.parse_inspiral_num_slides_slidespec(inspiral_slidespec)):
		time_slides[lsctables.TimeSlideTable.get_next_id()] = offsetvect

if options.verbose:
	print >>sys.stderr, "Total of %d time slides." % len(time_slides)


#
# Remove duplicates.
#


if options.verbose:
	print >>sys.stderr, "Identifying and removing duplicates ..."

map(time_slides.pop, ligolw_time_slide.time_slides_vacuum(time_slides, verbose = options.verbose).keys())

if options.verbose:
	print >>sys.stderr, "%d time slides remain." % len(time_slides)


#
# Remove zero-lag
#


if options.remove_zero_lag:
	if options.verbose:
		print >>sys.stderr, "Identifying and removing zero-lag ..."

	null_ids = [time_slide_id for time_slide_id, offsetvect in time_slides.items() if not any(offsetvect.deltas.values())]
	for time_slide_id in null_ids:
		del time_slides[time_slide_id]

	if options.verbose:
		print >>sys.stderr, "%d time slides remain." % len(time_slides)


#
# Convert to list of offset dictionaries.  We no longer require the IDs,
# new ones will be assigned as they are put into the output XML document.
#


time_slides = time_slides.values()


#
# Normalize the time slides
#


if options.normalize:
	if options.verbose:
		print >>sys.stderr, "Normalizing the time slides ..."
	constraints = parse_normalize(options.normalize)
	time_slides = [offsetvect.normalize(**constraints) for offsetvect in time_slides]


#
# Make documents.
#


lsctables.TimeSlideTable.reset_next_id()
while filenames:
	#
	# Create an empty document, populate the process information.
	#

	xmldoc, process = new_doc(**options.__dict__)
	timeslidetable = lsctables.TimeSlideTable.get_table(xmldoc)

	#
	# How many slides will go into this file?
	#

	N = int(round(float(len(time_slides)) / len(filenames)))

	#
	# Put them in.
	#

	for offsetvect in time_slides[:N]:
		timeslidetable.append_offsetvector(offsetvect, process)
	del time_slides[:N]

	#
	# Finish off the document.
	#

	ligolw_process.set_process_end_time(process)

	#
	# Write.
	#

	filename = filenames.pop(0)
	ligolw_utils.write_filename(xmldoc, filename, verbose = options.verbose, gz = (filename or "stdout").endswith(".gz"))

assert not time_slides
