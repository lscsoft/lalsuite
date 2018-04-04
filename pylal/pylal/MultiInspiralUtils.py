# Copyright (C) 2006  Sukanta Bose
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

import numpy
import itertools

from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from glue import segments
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import ilwd
from glue.ligolw import ligolw

#
# =============================================================================
#
#                                   Input
#
# =============================================================================
#

def ReadMultiInspiralFromFiles(fileList):
  """
  Read the multiInspiral tables from a list of files
  @param fileList: list of input files
  """
  if not fileList:
    return multiInspiralTable(), None

  multis = None

  for thisFile in fileList:
    doc = utils.load_filename(thisFile,
        gz=(thisFile or "stdin").endswith(".gz"), contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))
    # extract the multi inspiral table
    try:
      multiInspiralTable = table.get_table(doc,
          lsctables.MultiInspiralTable.tableName)
      if multis: multis.extend(multiInspiralTable)
      else: multis = multiInspiralTable
    except: multiInspiralTable = None
  return multis

def ReadMultiInspiralTimeSlidesFromFiles(fileList,generate_output_tables=False):
  """
  Read time-slid multiInspiral tables from a list of files
  @param fileList: list of input files
  """
  if not fileList:
    return multiInspiralTable(), None

  multis = None
  timeSlides = []

  segmentDict = {}
  for thisFile in fileList:

    doc = utils.load_filename(thisFile,
        gz=(thisFile or "stdin").endswith(".gz"), contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))
    # Extract the time slide table
    timeSlideTable = table.get_table(doc,
          lsctables.TimeSlideTable.tableName)
    slideMapping = {}
    currSlides = {}
    # NOTE: I think some of this is duplicated in the glue definition of the
    # time slide table. Probably should move over to that
    for slide in timeSlideTable:
      currID = int(slide.time_slide_id)
      if currID not in currSlides.keys():
        currSlides[currID] = {}
        currSlides[currID][slide.instrument] = slide.offset
      elif slide.instrument not in currSlides[currID].keys():
        currSlides[currID][slide.instrument] = slide.offset

    for slideID,offsetDict in currSlides.items():
      try:
        # Is the slide already in the list and where?
        offsetIndex = timeSlides.index(offsetDict)
        slideMapping[slideID] = offsetIndex
      except ValueError:
        # If not then add it
        timeSlides.append(offsetDict)
        slideMapping[slideID] = len(timeSlides) - 1

    # Get the mapping table
    segmentMap = {}
    timeSlideMapTable = table.get_table(doc,
        lsctables.TimeSlideSegmentMapTable.tableName)
    for entry in timeSlideMapTable:
      segmentMap[int(entry.segment_def_id)] = int(entry.time_slide_id)

    # Extract the segment table
    segmentTable = table.get_table(doc,
        lsctables.SegmentTable.tableName)
    for entry in segmentTable:
      currSlidId = segmentMap[int(entry.segment_def_id)]
      currSeg = entry.get()
      if not segmentDict.has_key(slideMapping[currSlidId]):
        segmentDict[slideMapping[currSlidId]] = segments.segmentlist()
      segmentDict[slideMapping[currSlidId]].append(currSeg)
      segmentDict[slideMapping[currSlidId]].coalesce()
    
    # extract the multi inspiral table
    try:
      multiInspiralTable = table.get_table(doc,
          lsctables.MultiInspiralTable.tableName)
      # Remap the time slide IDs
      for multi in multiInspiralTable:
        newID = slideMapping[int(multi.time_slide_id)]
        multi.time_slide_id = ilwd.ilwdchar(\
                              "time_slide:time_slide_id:%d" % (newID))
      if multis: multis.extend(multiInspiralTable)
      else: multis = multiInspiralTable
#    except: multiInspiralTable = None
    except: raise

  if not generate_output_tables:
    return multis,timeSlides,segmentDict
  else:
    # Make a new time slide table
    timeSlideTab = lsctables.New(lsctables.TimeSlideTable)

    for slideID,offsetDict in enumerate(timeSlides):
      for instrument in offsetDict.keys():
        currTimeSlide = lsctables.TimeSlide()
        currTimeSlide.instrument = instrument
        currTimeSlide.offset = offsetDict[instrument]
        currTimeSlide.time_slide_id = ilwd.ilwdchar(\
                                "time_slide:time_slide_id:%d" % (slideID))
        currTimeSlide.process_id = ilwd.ilwdchar(\
                                "process:process_id:%d" % (0))
        timeSlideTab.append(currTimeSlide)

    # Make a new mapping table
    timeSlideSegMapTab = lsctables.New(lsctables.TimeSlideSegmentMapTable)
    
    for i in range(len(timeSlides)):
      currMapEntry = lsctables.TimeSlideSegmentMap()
      currMapEntry.time_slide_id = ilwd.ilwdchar(\
                                "time_slide:time_slide_id:%d" % (i))
      currMapEntry.segment_def_id = ilwd.ilwdchar(\
                                "segment_def:segment_def_id:%d" % (i))
      timeSlideSegMapTab.append(currMapEntry)

    # Make a new segment table
    newSegmentTable = lsctables.New(lsctables.SegmentTable)

    segmentIDCount = 0
    for i in range(len(timeSlides)):
      currSegList = segmentDict[i]
      for seg in currSegList:
        currSegment = lsctables.Segment()
        currSegment.segment_id = ilwd.ilwdchar(\
                              "segment:segment_id:%d" %(segmentIDCount))
        segmentIDCount += 1
        currSegment.segment_def_id = ilwd.ilwdchar(\
                                "segment_def:segment_def_id:%d" % (i))
        currSegment.process_id = ilwd.ilwdchar(\
                                "process:process_id:%d" % (0))
        currSegment.set(seg)
        currSegment.creator_db = -1
        currSegment.segment_def_cdb = -1
        newSegmentTable.append(currSegment)
    return multis,timeSlides,segmentDict,timeSlideTab,newSegmentTable,\
           timeSlideSegMapTab


#
# =============================================================================
#
#                                  Clustering
#
# =============================================================================
#

def CompareMultiInspiralByEndTime(a, b):
  """
  Orders a and b by peak time.
  """
  return cmp(a.get_end(), b.get_end())


def CompareMultiInspiralBySnr(a, b):
  """
  Orders a and b by peak time.
  """
  return cmp(a.snr, b.snr)


def CompareMultiInspiral(a, b, twindow = LIGOTimeGPS(0)):
  """
  Returns 0 if a and b are less than twindow appart.
  """
  tdiff = abs(a.get_end() - b.get_end())
  if tdiff < twindow:
    return 0
  else:
    return cmp(a.get_end(), b.get_end())


def cluster_multi_inspirals(mi_table, dt, loudest_by="snr"):
    """Cluster a MultiInspiralTable with a given ranking statistic and
    clustering window.

    This method separates the table rows into time bins, returning those
    row that are
        * loudest in their own bin, and
        * louder than those events in the preceeding and following bins
          that are within the clustering time window

    @return: a new MultiInspiralTable containing those clustered events

    @param mi_table:
        MultiInspiralTable to cluster
    @param dt:
        width (seconds) of clustering window
    @keyword loudest_by:
        column by which to rank events, default: 'snr'

    @type mi_table: glue.ligolw.lsctables.MultiInspiralTable
    @type dt: float
    @type loudest_by: string
    @rtype: glue.ligolw.lsctables.MultiInspiralTable
    """
    cluster_table = table.new_from_template(mi_table)
    if not len(mi_table):
        return cluster_table

    # get data
    end_time = numpy.asarray(mi_table.get_end()).astype(float)
    if hasattr(mi_table, "get_%s" % loudest_by):
        stat = numpy.asarray(getattr(mi_table, "get_%s" % loudest_by)())
    else:
        stat = numpy.asarray(mi_table.get_column(loudest_by))

    # get times
    start = round(end_time.min())
    end = round(end_time.max()+1)

    # generate bins
    num_bins  = int((end-start)//dt + 1)
    time_bins = []
    loudest_stat = numpy.zeros(num_bins)
    loudest_time = numpy.zeros(num_bins)
    for n in range(num_bins):
        time_bins.append([])

    # bin triggers
    for i,(t,s) in enumerate(itertools.izip(end_time, stat)):
        bin_ = int(float(t-start)//dt)
        time_bins[bin_].append(i)
        if s > loudest_stat[bin_]:
            loudest_stat[bin_] = s
            loudest_time[bin_] = t

    # loop over all bins
    for i,bin_ in enumerate(time_bins):
        if len(bin_)<1:
            continue
        first = i==0
        last = i==(num_bins-1)

        prev = i-1
        next_ = i+1
        check_prev = (not first and len(time_bins[prev]) > 0)
        check_next = (not last and len(time_bins[next_]) > 0)

        # pick loudest event in bin
        idx = bin_[stat[bin_].argmax()]
        s = stat[idx]
        t = end_time[idx]

        # trigger was loudest in it's bin, search loudest event
        # in previous bin
        if (check_prev and (t - loudest_time[prev]) < dt and
                s < loudest_stat[prev]):
            continue

        # Same for the next bin
        if (check_next and (loudest_time[next_] - t) < dt and
                s < loudest_stat[next_]):
            continue

        loudest=True

        # trigger was loudest in it's bin, search previous bin
        if check_prev and not (t - loudest_time[prev]) < dt:
            for idx2 in time_bins[prev]:
                t2 = end_time[idx2]
                if (t - end_time[idx2]) < dt and s < stat[idx2]:
                    loudest = False
                    break
        if not loudest:
            continue

        # if still loudest, check the next bin
        if check_next and not (loudest_time[next_] - t) < dt:
            for idx2 in time_bins[next_]:
                if (end_time[idx2] - t) < dt and s < stat[idx2]:
                    loudest = False
                    break
        if not loudest:
            continue

        # this was the loudest trigger in its vicinity,
        # keep it and move to the next bin
        cluster_table.append(mi_table[idx])

    return cluster_table
