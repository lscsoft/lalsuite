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


from __future__ import print_function


import math
from optparse import OptionParser
import sys
from PIL import Image


from glue.text_progress_bar import ProgressBar
from ligo.lw import ligolw
from ligo.lw import lsctables
from ligo.lw import utils as ligolw_utils
from ligo.lw.utils import process as ligolw_process
import lalburst
import lalmetaio
import lalsimulation
from lalburst import git_version


@lsctables.use_in
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	pass


class SimBurst(lalmetaio.SimBurst):
	def __init__(self, **kwargs):
		super(SimBurst, self).__init__()
		for key, value in kwargs.items():
			setattr(self, key, value)
	@property
	def time_geocent_gps_ns(self):
		return self.time_geocent_gps.gpsNanoSeconds
lsctables.SimBurst = lsctables.SimBurstTable.RowType = SimBurst


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
	parser.add_option("-l", "--f-low", metavar = "Hz", type = "float", default = 64.0, help = "Set the low-frequency limit of the tiling in Hertz (default = 64).")
	parser.add_option("-d", "--delta-f", metavar = "Hz", type = "float", default = 16.0, help = "Set the bandwidth of the pixels in Hertz (default = 16).  Must be > 0.  The product of --delta-f and --delta-t must be at least 2/pi.")
	parser.add_option("-t", "--delta-t", metavar = "s", type = "float", default = 1.0 / 16, help = "Set the duration of the pixels in seconds (default = 1/16).  Must be > 0.  The product of --delta-f and --delta-t must be at least 2/pi.")
	parser.add_option("-f", "--overlap-fraction", metavar = "fraction", type = "float", default = 0.25, help = "Set the fraction by which adjacent tiles overlap (default = 0.25).  The pixels centres are spaced (1 - --overlap-fraction) * --delta-f apart in the frequency domain.  The value must be in the range [0, 1).")
	parser.add_option("-H", "--height", metavar = "pixels", type = "int", default = 64, help = "Set the number of pixles in the frequency domain (default = 64).  The image will be scaled to this vertical size, and the number of pixels in the time domain (horizontal size) will be fixed by the image aspect ratio.")
	parser.add_option("-o", "--output", metavar = "filename", help = "Set the name of the output file (default = stdout).")
	parser.add_option("-s", "--gps-start-time", metavar = "seconds", help = "Set the start time of the tiling in GPS seconds (required).")
	parser.add_option("--sample-rate", metavar = "Hz", type = "int", default = 16384, help = "Set the sample rate in Hertz of the data into which the injections will be placed (default = 16384).  This information is required in order to normalize each pixel accurately.  If the wrong value is used, the result will be the addition of noise to the image.  The highest frequency pixel must have a centre frequency < 1/2 this frequency.")
	parser.add_option("-n", "--hrss-scale", metavar = "hrss", type = "float", default = 1e-20, help = "Set the single-pixel hrss (root-sum-square strain) corresponding to \"white\" (default = 1e-20).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	parser.add_option("-T", "--time-slide-xml", metavar = "filename", help = "Associate injections with the first time slide ID in this XML file (required).")
	options, filenames = parser.parse_args()

	# save for the process_params table
	options.options_dict = dict(options.__dict__)

	if options.gps_start_time is None:
		raise ValueError("missing required option --gps-start-time")
	if not 0. < options.delta_f:
		raise ValueError("--delta-f must be > 0")
	if not 0. < options.delta_t:
		raise ValueError("--delta-t must be > 0")
	if options.delta_f * options.delta_t < 2. / math.pi:
		raise ValueError("the product of --delta-t and --delta-f must be >= 2/pi")
	if not (0. <= options.overlap_fraction < 1.):
		raise ValueError("--overlap-fraction must be in [0, 1)")
	if options.f_low + options.height * (1. - options.overlap_fraction) * options.delta_f >= options.sample_rate / 2.:
		raise ValueError("total waveform bandwidth exceeds Nyquist frequency")
	if options.time_slide_xml is None:
		raise ValueError("missing required option time-slide-xml")

	# type-cast
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


#
# use the time-slide file to start the output document
#


xmldoc = ligolw_utils.load_filename(options.time_slide_xml, verbose = options.verbose, contenthandler = LIGOLWContentHandler)


#
# add our metadata
#


process = ligolw_process.register_to_xmldoc(xmldoc, u"lalapps_binj_pic", options.options_dict, version = git_version.version)


#
# use whatever time slide vector comes first in the table (lazy)
#


time_slide_table = lsctables.TimeSlideTable.get_table(xmldoc)
time_slide_id = time_slide_table[0].time_slide_id
if options.verbose:
	print("associating injections with time slide (%d) %s" % (time_slide_id, time_slide_table.as_dict()[time_slide_id]), file = sys.stderr)


#
# find or add a sim_burst table
#


try:
	lsctables.SimBurstTable.get_table(xmldoc)
except ValueError:
	# no sim_burst table in document
	pass
else:
	raise ValueError("%s contains a sim_burst table.  this program isn't smart enough to deal with that." % options.time_slide_xml)

sim_burst_tbl = xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.SimBurstTable, ["process:process_id", "simulation_id", "time_slide:time_slide_id", "waveform", "waveform_number", "ra", "dec", "psi", "time_geocent_gps", "time_geocent_gps_ns", "duration", "frequency", "bandwidth", "egw_over_rsquared", "pol_ellipse_angle", "pol_ellipse_e"]))


#
# populate the sim_burst table with pixels
#


if options.verbose:
	print("time-frequency tiles have %g degrees of freedom" % (2. * options.delta_t * options.delta_f), file = sys.stderr)

for n, filename in enumerate(filenames, 1):
	if options.verbose:
		print("%d/%d: loading %s ..." % (n, len(filenames), filename), file = sys.stderr)
	img = Image.open(filename)

	width, height = img.size
	width, height = int(round(width * options.height / float(height))), options.height
	if options.verbose:
		print("converting to %dx%d grayscale ... " % (width, height), file = sys.stderr)
	img = img.resize((width, height)).convert("L")

	progress = ProgressBar("computing pixels", max = width * height) if options.verbose else None
	for i in range(width):
		for j in range(height):
			if progress is not None:
				progress.increment()
			# amplitude.  hrss column is ignored by waveform
			# generation code.  it is included for convenience,
			# to record the desired pixel brightness.  because
			# band- and time-limited white-noise burst
			# waveforms are random, the waveform's final
			# ampltiude (in the egw_over_rsquared column) is
			# determined by generating the burst at a canonical
			# amplitude, measuring its hrss, then rescaling to
			# achieve the desired value.  this process requires
			# the final sample rate to be known.
			hrss = options.hrss_scale * img.getpixel((i, height - 1 - j)) / 255.0
			if hrss == 0.:
				# don't generate injections for black
				# pixels
				continue

			# create and initialize table row
			row = lsctables.SimBurst(
				# metadata
				process_id = process.process_id,
				simulation_id = sim_burst_tbl.get_next_id(),
				time_slide_id = time_slide_id,

				# source orientation
				ra = 0.,
				dec = 0.,
				psi = 0.,

				# band- and time-limited white-noise burst
				waveform = "BTLWNB",
				waveform_number = 0,
				pol_ellipse_e = 0.,
				pol_ellipse_angle = 0.,

				# time-frequency co-ordinates
				time_geocent_gps = options.gps_start_time + i * options.delta_t * (1.0 - options.overlap_fraction),
				frequency = options.f_low + j * options.delta_f * (1.0 - options.overlap_fraction),
				bandwidth = options.delta_f,
				duration = options.delta_t,

				# brightness.  hrss is ignored, set
				# egw_over_rsquared to 1, to be re-scaled
				# after measuring waveform's real
				# brightness
				hrss = hrss,
				egw_over_rsquared = 1.
			)

			# generate waveform, then scale egw_over_rsquared
			# to get desired brightness
			row.egw_over_rsquared *= hrss / lalsimulation.MeasureHrss(*lalburst.GenerateSimBurst(row, 1.0 / options.sample_rate))

			# put row into table
			sim_burst_tbl.append(row)
	del progress


#
# write output
#


ligolw_utils.write_filename(xmldoc, options.output, gz = (options.output or "stdout").endswith(".gz"), verbose = options.verbose)
