#
# Copyright (C) 2013  Branson Stephens and Chris Pankow
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

import pygtk
pygtk.require("2.0")
import pygst
pygst.require("0.10")
import gobject
import gst

from gstlal import simplehandler
from gstlal import pipeparts
from laldetchar.idq import vetosrc
from optparse import OptionParser, Option
import os

#
# =============================================================================
#
# Pad probes
# XXX This is an ugly hack.
# The caps filter is getting rid of my gap buffer flag. 
# I want to get it back.
#
# =============================================================================
#

class handlerClass():
    def __init__(self):
        self.gap_list = []

    def earlyBufferHandler(self, pad, gst_buffer):
        if gst_buffer.flag_is_set(gst.BUFFER_FLAG_GAP):
            self.gap_list.append(gst_buffer.timestamp) 
        return True

    def lateBufferHandler(self, pad, gst_buffer):
        if gst_buffer.timestamp in self.gap_list:
            gst_buffer.flag_set(gst.BUFFER_FLAG_GAP)
            self.gap_list.remove(gst_buffer.timestamp)
        return True

#
# =============================================================================
#
#                                   Options
#
# =============================================================================
#

parser = OptionParser(description = __doc__)

# General options
parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose (optional).")
parser.add_option("--log-file", metavar = "name", default = "/tmp/idq_sender.log", help = "Full path to log for sent frames.")
parser.add_option("--output-type", metavar = "name", help = "Method of output. Valid choices are files 'files', tcp (default), and fake.", default = "tcp")

# Options for accessing iDQ results
parser.add_option("--init-time", metavar = "s", default = 0, type = "int", help = "GPS time of first file to ingest. Note: default of zero will be replaced with 'now'")
parser.add_option("--wait-time", metavar = "s", default = 32, type = "int", help = "Time to wait for next input file before giving up.")
parser.add_option("--dir-digits", metavar = "s", default = 5, type = "int", help = "The number of digits of the gps time in the directory name")
parser.add_option("--input-path", metavar = "name", help = "Path to input numpy files.")
parser.add_option("--input-prefix", metavar = "name", help = "Prefix for numpy files.")
parser.add_option("--input-ext", metavar = "name", help = "Extension for numpy files.", default = ".npy.gz")
parser.add_option("--idq-log-file", metavar = "name", help = "Full path to log for for the iDQ realtime process.")

# Options for the muxer
parser.add_option("--frame-duration", metavar = "s", default = 32, type = "int", help = "Set the duration of the output frames")
parser.add_option("--frames-per-file", metavar = "s", default = 1, type = "int", help = "Output frames per file")
parser.add_option("--channel-name", metavar = "name", help = "Output channel name.")

# Options for file output
parser.add_option("--frame-type", metavar = "name", help = "Specify the non-instrumental part of the frame type. The full frame type will be constructed by prepending the instrument.")
parser.add_option("--instrument", metavar = "name", help = "Specify the instrumental part of the frame type.")
parser.add_option("--output-path", metavar = "name", help = "Path to output frame files.")

# Options for the tcpclientsink
parser.add_option("--port", metavar = "s", default = 4953, type = "int", help = "Port for tcpclientsink")

options, filenames = parser.parse_args()

#
# =============================================================================
#
#                                    Main
#
# =============================================================================
#

#
# Set the pipeline up
#

pipeline = gst.Pipeline("veto_source")
mainloop = gobject.MainLoop()
handler = simplehandler.Handler(mainloop,pipeline)

# Check the instrument from the input_prefix
# XXX Used to be something like 'L-...' and now its 'L1_....'
#obsStr = options.input_prefix.split('-')[0]
obsStr = options.input_prefix[0]
if not options.instrument.startswith(obsStr):
    raise ValueError("Output channel instrument clashes with input prefix.")

# Setup the source class
if options.init_time > 0:
    init_time = options.init_time
else:
    init_time = None
vsrc = vetosrc.vetoSource(options.input_path, options.input_prefix,
            options.input_ext, options.wait_time, init_time, options.dir_digits, 
            options.log_file, options.idq_log_file)

# Create the appsrc with accoutrements
appsrc = pipeparts.mkgeneric(pipeline, None, "appsrc", caps=gst.Caps(vsrc.caps), 
    format="time")
appsrc.connect('need-data', vsrc.need_data)

# Set debug level for logging purposes
gst.debug_set_threshold_for_name('python', gst.LEVEL_INFO)

# Define the muxer.
mux = pipeparts.mkframecppchannelmux(pipeline, None, 
    frames_per_file = options.frames_per_file, 
    frame_duration = options.frame_duration)

# Link the source to the muxer. 
appsrc.get_pad("src").link(mux.get_pad(options.instrument + ':' + options.channel_name))

# XXX Hacking. Attach probe to the muxer.
hc = handlerClass()
mux_sink = mux.get_pad(options.instrument + ':' + options.channel_name)
mux_sink.add_buffer_probe(hc.earlyBufferHandler)

# Final destination.
if options.output_type == "files":
    try:
        os.makedirs(options.output_path)
    except Exception as e:
        print "Failed with %s" % e

    # Inject tags.  The framecpp_filesink element uses the tags to figure
    # out the output filename.
    print "Setting tag instrument: %s" % options.instrument
    tagInj = pipeparts.mktaginject(pipeline, mux, 
        "instrument=%s,channel-name=%s" % (options.instrument, options.channel_name))

    path = options.output_path
    if path:
        fs = pipeparts.mkframecppfilesink(pipeline, tagInj, 
            frame_type = options.frame_type, path = options.output_path)
    else:
        fs = pipeparts.mkframecppfilesink(pipeline, mux, 
            frame_type = options.frame_type)
elif options.output_type == "tcp":
    # NB: like the appsrc, the tcp client sink does not use the 'blocksize' property. 
    # It is just hanging around because it is inherited
    # gstdataprotocol serializes and sends the caps first. Then serializes and 
    # sends the buffers. Less work for us.

    # XXX Sadness. The capsfilter is not "gap-aware". So the gap flag will be unset.
    mux = pipeparts.mkcapsfilter(pipeline,mux,caps="application/x-igwd-frame,framed=true")
    tcpclientsink = pipeparts.mkgeneric(pipeline, mux, "tcpclientsink")
    tcpclientsink.set_property("protocol", "GST_TCP_PROTOCOL_GDP")
    tcpclientsink.set_property("port", options.port)
    tcpclientsink.set_property("sync", False)
    # XXX Hack. Attach pad probe to sink pad of tcpclientsink
    tcs_pad = tcpclientsink.get_pad("sink")
    tcs_pad.add_buffer_probe(hc.lateBufferHandler)
else:
    raise ValueError("Invalid output type.")

#
# Start the thing going.
#

pipeline.set_state(gst.STATE_PLAYING)
mainloop.run()

#
# done
#
