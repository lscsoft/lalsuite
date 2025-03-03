##python
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


import os.path
import sys
from optparse import OptionParser


import igwn_segments as segments


from lal.utils import CacheEntry


from lal import git_version


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
		#version = "Name: %%prog\n%s" % git_version.verbose_msg,
		usage = "usage: %prog [options]\n\nExample:\n\tls *.xml | %prog"
	)
	parser.add_option("-a", "--include-all", action = "store_true", help = "Include all files in output.  Unparseable file names are assigned empty metadata.")
	parser.add_option("-f", "--force", action = "store_true", help = "Ignore errors.  Unparseable file names are removed from the output.  This has no effect if --include-all is given.")
	parser.add_option("-i", "--input", metavar="filename", help="Read input from this file (default = stdin).")
	parser.add_option("-o", "--output", metavar="filename", help="Write output to this file (default = stdout).")
	parser.add_option("-v", "--verbose", action="store_true", help="Be verbose.")
	return parser.parse_args()[0]


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


#
# Parse command line
#


options = parse_command_line()


#
# Open input and output streams
#


if options.input is not None:
	src = open(options.input)
else:
	src = sys.stdin

if options.output is not None:
	dst = open(options.output, "w")
else:
	dst = sys.stdout

#
# Other initializations
#


path_count = 0
seglists = segments.segmentlistdict()


#
# Filter input one line at a time
#


for line in src:
	path, filename = os.path.split(line.strip())
	url = "file://localhost%s" % os.path.abspath(os.path.join(path, filename))
	try:
		cache_entry = CacheEntry.from_T050017(url)
	except ValueError as e:
		if options.include_all:
			cache_entry = CacheEntry(None, None, None, url)
		elif options.force:
			continue
		else:
			raise e
	print(str(cache_entry), file=dst)
	path_count += 1
	if cache_entry.segment is not None:
		seglists |= cache_entry.segmentlistdict.coalesce()


#
# Summary
#


if options.verbose:
	print("Size of cache: %d URLs" % path_count, file=sys.stderr)
	for instrument, seglist in list(seglists.items()):
		ext = seglist.extent()
		dur = abs(seglist)
		print("Interval spanned by %s: %s (%s s total, %.4g%% duty cycle)" % (instrument, str(ext), str(dur), 100.0 * float(dur) / float(abs(ext))), file=sys.stderr)
	span = seglists.union(seglists)
	ext = span.extent()
	dur = abs(span)
	print("Interval spanned by union: %s (%s s total, %.4g%% duty cycle)" % (str(ext), str(dur), 100.0 * float(dur) / float(abs(ext))), file=sys.stderr)
