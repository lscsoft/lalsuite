#!/usr/bin/python
#
# Copyright (C) 2013 Branson Stephens and Chris Pankow
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

import glob
import gzip
import os
import time
import logging
import datetime
import json
from collections import deque
import numpy as np
from glue import gpstime

from laldetchar.idq import idq

#
# =============================================================================
#
#                                Utility Functions
#
# =============================================================================
#


# Older versions of numpy.load choke on .gz files.
# Not sure what version fixes this, but the numpy package 
# for debian squeeze (1.4.1) doesn't. See
# https://github.com/numpy/numpy/issues/1593

def wrapNpLoad(filePath):
    if filePath.endswith('gz'):
        fileObj = gzip.open(filePath)
        payload = np.load(fileObj)
        fileObj.close()
        return payload
    else:
        return np.load(filePath)


#
# =============================================================================
#
#                            Pipeline Utility Functions
#
# =============================================================================
#

# Pad probe handler.  Useful for debugging.
def probeBufferHandler(pad,gst_buffer):
    print 'gpsstart  = %s' % gst_buffer.timestamp
    print 'offset    = %s' % gst_buffer.offset
    print 'is gap    = %s' % gst_buffer.flag_is_set(gst.BUFFER_FLAG_GAP)
    return True

# A class for handling the veto timeseries files.

class vetoSource:
    def __init__(self, inputPath, inputPre, inputExt, waitTime, initTime=None, dirDigits=5, 
        log_file=None, idq_log_file=None):
        self.inputPath  = inputPath
        self.inputPre   = inputPre
        self.inputExt   = inputExt
        self.waitTime   = waitTime
        self.dirDigits  = dirDigits
        self.fullCurrentPrefix = ""

        self.logger = None
        if log_file:
            self.logger = logging.getLogger('idq_sender')
            self.logger.setLevel(logging.INFO)

            fh = logging.FileHandler(log_file)
            self.logger.addHandler(fh)
            fh.setFormatter(logging.Formatter('%(message)s'))

        # This is the log file of the running iDQ realtime process.
        # We will use it to see which stride iDQ is currently working on.
        self.idq_log = None
        if idq_log_file:
            self.idq_log = open(idq_log_file, "r")

        # Initialize the list of files.
        if not initTime:
            # FIXME:  Badness.  Sutracting max iDQ latency by hand.
            initTime = int(gpstime.GpsSecondsFromPyUTC(time.time())) - (2*64 + 2*128)
        self.check_for_new_files(initTime*gst.SECOND)
        # This is the offset of a given buffer with respect to the global stream in 
        # units of *samples*.
        self.current_offset = 0

        # Determine the rate by looking at the first file in the list.
        # XXX If the files with a given prefix have different rates, this will break.
        # Assume file names are of the form:
        # (old:)
        # L-KW_TRIGGERS_ovl_1058313600_1058400000-1058411456-64.npy.gz
        # <input_prefix><training_start>_<training_end>-<gpsstart>-<duration><input_ext>
        # (new:)
        # L1_idq_ovl_1117544400_1117548000-1117498800_1117542600_rank-1117552416-32.npy.gz
        # <input_prefix><something number>_<training_start>_<training_end>_rank-<gpsstart>-<duration><input_ext>
        # Ex:
        # get the first one.
        filePath = self.fileQueue.popleft()
        # put it back, since this is just for informational purposes.
        self.fileQueue.extendleft([filePath])
        firstVals = wrapNpLoad(filePath)
        #training_period, gps_start, rest = filePath[len(self.fullCurrentPrefix):].split('-')
        something, something_else, gps_start, rest = filePath[len(self.fullCurrentPrefix):].split('-')

        self.duration = int(rest[:rest.find(inputExt)])
        self.next_output_timestamp = int(gps_start) * gst.SECOND

        # The rate determination assumes the file contains N values with
        # t_i = gpsstart + (i-1)*Delta_t ,
        # such that the last data point is at time t = (gpsstart+duration)-Delta_t.
        # NOTE gstreamer expects an integer rate in Hz.
        # XXX This cast makes me uncomfortable.  That's why I added 0.1.
        self.rate = int(float(firstVals.size)/float(self.duration) + 0.1)
        # Now that we have the rate, we can set the caps for the appsrc
        self.caps = "audio/x-raw-float,width=32,channels=1,rate=%d" % self.rate

    def check_for_new_files(self, timestamp=0):
        # Figure out which directory the files should be in.
        timeStr = str(timestamp/gst.SECOND)
        if len(timeStr)>=self.dirDigits:
            dirSuffix = timeStr[:self.dirDigits]
        else:
            # What is this? The 1980s? Let's just pad with zeros.
            dirSuffix = timeStr.zfill(self.dirDigits)
        inputPath = self.inputPath + dirSuffix
        self.fullCurrentPrefix = os.path.join(inputPath,self.inputPre)
        # Excluding anything without 'rank' in the middle
        pattern = self.fullCurrentPrefix + '*rank*' + self.inputExt

        def is_current_file(path):
            filePath = os.path.basename(path)
            rest = filePath[len(self.inputPre):]
            if len(rest)>0:
                return int(rest.split('-')[2])*gst.SECOND>=timestamp
            else:
                return None
        filePathList = filter(is_current_file, glob.glob(pattern))
        filePathList.sort()
        self.fileQueue = deque(filePathList)
        return

    # 
    # Push a gap into the given appsrc element, given a duration in seconds.
    # The beginning of the gap period is assumed to be self.next_output_timestamp
    # Notice that next_output_timestamp is reset according to the duration of the
    # gap. We also write a line of JSON in a log indicating that a gap has been pushed.
    #
    def push_gap(self, src, gap_duration_s):
        gap_start = self.next_output_timestamp / gst.SECOND
        gap_duration = gap_duration_s * gst.SECOND
        gap_samples = gap_duration_s * self.rate
        gap_end = gap_start + self.waitTime
        # FIXME
        # For a real gap buffer, you would want zero samples. (I think.)
        #gap_samples = 0
        gap_vals = np.zeros(gap_samples)
        gap_vals = gap_vals.astype(np.float32)
        buffer_len = gap_vals.nbytes
        buf = gst.buffer_new_and_alloc(buffer_len)
        buf[:buffer_len-1] = np.getbuffer(gap_vals)
        # Again, for a real gap buffer, you'd want to set this flag.
        #buf.flag_set(gst.BUFFER_FLAG_GAP)
        buf.timestamp = self.next_output_timestamp
        buf.duration = gap_duration
        buf.offset = self.current_offset
        buf.offset_end = self.current_offset + gap_samples
        src.emit("push-buffer", buf)
        src.info("gst clock = %d" % int(src.get_clock().get_time()))
        src.info("No files! Pushed gap with start=%d, duration=%d latency=%d" %
            (buf.timestamp/gst.SECOND,gap_duration/gst.SECOND,(gpstime.GpsSecondsFromPyUTC(time.time()) - gap_end)))

        if self.logger:
            outDict = {}
            outDict['type'] = 'buffer'
            outDict['time'] = datetime.datetime.now().isoformat()
            outDict['is_gap'] = True
            outDict['gpsstart'] = buf.timestamp/gst.SECOND
            outDict['duration'] = gap_duration/gst.SECOND
            outDict['latency'] = gpstime.GpsSecondsFromPyUTC(time.time()) - gap_end
            self.logger.info(json.dumps(outDict))
        self.next_output_timestamp += buf.duration
        self.current_offset = buf.offset_end
        return

    #
    # Respond to the need_data signal from the appsrc. 
    # 
    def need_data(self, src, need_bytes=None):
        src.info("----------------------------------------------------")
        src.info("Received need-data signal, %s." % time.asctime())

        self.check_for_new_files(self.next_output_timestamp)
        try:
            filePath = self.fileQueue.popleft()
        except IndexError:
            if self.idq_log:
                # Check to see what is the latest stride that iDQ is working on.
                stride_start, stride_end = idq.most_recent_realtime_stride(self.idq_log)

                if stride_start * gst.SECOND > self.next_output_timestamp:
                    # The iDQ process has moved on. Push a gap.
                    gap_duration_s = stride_start - self.next_output_timestamp / gst.SECOND
                    self.push_gap(src, gap_duration_s)
            
            # If no iDQ log file is available, or if the stride we want is still being 
            # worked on by iDQ, we wait.
            time.sleep(self.waitTime)
            # Try it again.
            self.check_for_new_files(self.next_output_timestamp)
            try:
                filePath = self.fileQueue.popleft()
            except IndexError: 
                # Push gap equivalent to the wait time and return.
                self.push_gap(src, self.waitTime)
                return True
                    
        # Ah, we have a file.
        # Get the gpsstart time from the filename.
        rest = filePath[len(self.fullCurrentPrefix):]
        gpsstart = int(rest.split('-')[2])
        # Let's re-derive the duration.  maybe it changed?
        rest = rest.split('-')[3]
        duration = int(rest[:rest.find(self.inputExt)])

        # Is this file later than the one we were expecting?  
        # If so, we can be sure that the one we wanted will never show up.
        # We might as well push a gap and then continue to process
        # the file that we have available.
        if gpsstart * gst.SECOND > self.next_output_timestamp:
            gap_duration_s = gpsstart - self.next_output_timestamp / gst.SECOND
            self.push_gap(src, gap_duration_s)

        # Down to business. We must read in data from the npy file in order to
        # create the buffer. First load the numpy array.
        src.info("processing %s" % filePath)
        veto_vals = wrapNpLoad(filePath)
        veto_vals = veto_vals.astype(np.float32)

        # Build the buffer.
        buffer_len = veto_vals.nbytes
        buf = gst.buffer_new_and_alloc(buffer_len)
        buf[:buffer_len-1] = np.getbuffer(veto_vals)
        buf.timestamp = gpsstart * gst.SECOND
        # gst buffers require:
        # buffer_duration * rate / gst.SECOND = (offset_end - offset)
        # The offset is zero since our data begin at the beginning 
        # of the buffer.
        buf_end = gpsstart + duration
        buf.duration = duration * gst.SECOND
        buf.offset = self.current_offset
        buf.offset_end = self.current_offset + duration * self.rate
        buf.caps = self.caps
        # Push the buffer into the stream (a side effect of 
        # emitting this signal).
        src.emit("push-buffer", buf)
        # XXX FIXME Use a real gstreamer clock for latency instead? Make sure to fix it in 
        # the other places too. Okay?
        #src.info("pushed buffer with start=%d, duration=%d, latency=%d" % 
        #    (gpsstart,duration, (int(src.get_clock().get_time())-buf.timestamp)/gst.SECOND))
        src.info("pushed buffer with start=%d, duration=%d, latency=%d" % 
            (gpsstart,duration,(gpstime.GpsSecondsFromPyUTC(time.time()) - buf_end)))
        if self.logger:
            outDict = {}
            outDict['type'] = 'buffer'
            outDict['time'] = datetime.datetime.now().isoformat()
            outDict['is_gap'] = False
            outDict['gpsstart'] = gpsstart
            outDict['duration'] = duration
            outDict['latency'] = gpstime.GpsSecondsFromPyUTC(time.time()) - buf_end
            self.logger.info(json.dumps(outDict))        
    
        self.next_output_timestamp += buf.duration
        self.current_offset = buf.offset_end
        return True



