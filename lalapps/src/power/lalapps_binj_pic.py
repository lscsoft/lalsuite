#
# Copyright (C) 2010,2013  Kipp Cannon
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


import math
from optparse import OptionParser
import sys
from PIL import Image


from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
import lalburst
import lalmetaio
import lalsimulation
from lalburst import git_version


class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	pass
lsctables.use_in(LIGOLWContentHandler)


# FIXME:  the lalmetaio.SimBurst type needs to be subclassed to get the ID
# columns to behave like ilwd:char strings.
lsctables.SimBurst = lsctables.SimBurstTable.RowType = lalmetaio.SimBurst


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
		usage = "%prog [options] filename",
		description = "Convert an image into a LIGO Light-Weight XML file containing a list of sine-Gaussian burst injections.  When injected into data, the injections will cause a waterfall plot to display the image."
	)
	parser.add_option("-l", "--f-low", metavar = "Hz", type = "float", default = 64.0, help = "Set the low-frequency limit of the tiling (default = 64).")
	parser.add_option("-d", "--delta-f", metavar = "Hz", type = "float", default = 16.0, help = "Set the frequency spacing of the tiling (default = 16).")
	parser.add_option("-t", "--delta-t", metavar = "s", type = "float", default = 1.0 / 16, help = "Set the time spacing of the tiling (default = 1/16).")
	parser.add_option("-H", "--height", metavar = "pixels", type = "int", default = 64, help = "Set the number of tiles in the frequency domain (default = 64).")
	parser.add_option("-o", "--output", metavar = "filename", help = "Set the name of the output file (default = stdout).")
	parser.add_option("-s", "--gps-start-time", metavar = "seconds", help = "Set the start time of the tiling in GPS seconds (required).")
	parser.add_option("-f", "--overlap-fraction", metavar = "fraction", type = "float", default = 0.25, help = "Set the fraction by which adjacent tiles overlap (default = 0.25).  The value must be in the range [0, 1).")
	parser.add_option("--sample-rate", metavar = "Hz", type = "int", default = 16384, help = "Set the sample rate of the data into which the injections will be placed (default = 16384).  This information is required in order to normalize each pixel accurately.  If the wrong value is used, the result will be the addition of noise to the image.")
	parser.add_option("-n", "--hrss-scale", metavar = "hrss", type = "float", default = 1e-20, help = "Set the hrss corresponding to \"white\" (default = 1e-20).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	parser.add_option("-T", "--time-slide-xml", metavar = "filename", help = "Associate injections with the first time slide ID in this XML file (required).")
	options, filenames = parser.parse_args()

	if options.gps_start_time is None:
		raise ValueError("missing required option --gps-start-time")
	if not (0 <= options.overlap_fraction < 1.0):
		raise ValueError("--overlap-fraction must be in [0, 1)")
	if options.delta_f * options.delta_t < 2 / math.pi:
		raise ValueError("the product of --delta-t and --delta-f must be >= 2/pi")
	if options.time_slide_xml is None:
		raise ValueError("missing required option time-slide-xml")

	# for the process_params table
	options.options_dict = dict(options.__dict__)

	options.gps_start_time = lsctables.LIGOTimeGPS(options.gps_start_time)

	return options, filenames


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


options, filenames = parse_command_line()


if options.verbose:
	print >>sys.stderr, "time-frequency tiles have %g degrees of freedom" % (2 * options.delta_t * options.delta_f)


xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())
process = ligolw_process.register_to_xmldoc(xmldoc, u"lalapps_binj_pic", options.options_dict, version = git_version.version)
time_slide_table = xmldoc.childNodes[-1].appendChild(lsctables.TimeSlideTable.get_table(ligolw_utils.load_filename(options.time_slide_xml, verbose = options.verbose, contenthandler = LIGOLWContentHandler)))
time_slide_id = time_slide_table[0].time_slide_id
if options.verbose:
	print >>sys.stderr, "associating injections with time slide ID \"%s\":  %s" % (time_slide_id, time_slide_table.as_dict()[time_slide_id])
for row in time_slide_table:
	row.process_id = process.process_id
sim_burst_tbl = xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.SimBurstTable, ["process_id", "simulation_id", "time_slide_id", "waveform", "waveform_number", "ra", "dec", "psi", "time_geocent_gps", "time_geocent_gps_ns", "duration", "frequency", "bandwidth", "egw_over_rsquared", "pol_ellipse_angle", "pol_ellipse_e"]))


for filename in filenames:
	if options.verbose:
		print >>sys.stderr, "loading %s ..." % filename
	img = Image.open(filename)

	width, height = img.size
	width, height = int(round(width / float(height) * options.height)), options.height
	if options.verbose:
		print >>sys.stderr, "converting to %dx%d grayscale ... " % (width, height)
	img = img.resize((width, height)).convert("L")

	for i in xrange(width):
		for j in xrange(height):
			# new row
			row = lsctables.SimBurst()
			row.process_id = int(process.process_id)
			row.simulation_id = int(sim_burst_tbl.get_next_id())
			row.time_slide_id = int(time_slide_id)

			# source orientation
			row.ra = row.dec = row.psi = 0

			# band- and time-limited white-noise burst
			row.waveform = "BTLWNB"
			row.waveform_number = 0
			row.pol_ellipse_e = row.pol_ellipse_angle = 0

			# time-frequency co-ordinates
			row.frequency = options.f_low + j * options.delta_f * (1.0 - options.overlap_fraction)
			gps = options.gps_start_time + i * options.delta_t * (1.0 - options.overlap_fraction)
			row.time_geocent_gps, row.time_geocent_gps_ns = gps.seconds, gps.nanoseconds
			row.bandwidth = options.delta_f
			row.duration = options.delta_t

			# amplitude.  hrss column is ignored by waveform
			# generation code.  it is included for convenience.
			# the waveform is normalized by generating the
			# burst, measuring its hrss, then rescaling to
			# achieve the desired value.  this process requires
			# the sample rate to be known.
			row.hrss = options.hrss_scale * img.getpixel((i, height - 1 - j)) / 255.0
			if row.hrss != 0:
				row.egw_over_rsquared = 1
				row.egw_over_rsquared *= row.hrss / lalsimulation.MeasureHrss(*lalburst.GenerateSimBurst(row, 1.0 / 16384))
				sim_burst_tbl.append(row)
			if options.verbose:
				print >>sys.stderr, "generating sim_burst table ... %d injections\r" % len(sim_burst_tbl),
	if options.verbose:
		print >>sys.stderr


ligolw_utils.write_filename(xmldoc, options.output, gz = (options.output or "stdout").endswith(".gz"), verbose = options.verbose)
