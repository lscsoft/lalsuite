# Copyright (C) 2006  Duncan A. Brown
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

import copy

from pylal import SearchSummaryUtils
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw.utils import ligolw_add
from glue import iterutils
from glue import segments


#
# =============================================================================
#
#                                   Input
#
# =============================================================================
#

class ExtractSnglInspiralTableLIGOLWContentHandler(ligolw.PartialLIGOLWContentHandler):
  """
  LIGOLWContentHandler that will extract only the SnglInspiralTable from a document.
  See glue.ligolw.LIGOLWContentHandler help for more info.
  """
  def __init__(self,document):
    def filterfunc(name,attrs):
      if name==ligolw.Table.tagName and attrs.has_key('Name'):
        return 0==table.CompareTableNames(attrs.get('Name'), lsctables.SnglInspiralTable.tableName)
      else:
        return False
    ligolw.PartialLIGOLWContentHandler.__init__(self,document,filterfunc)


class SnglInspiralID_old(object):
  """
  Custom row ID thing for sngl_inspiral tables with int_8s event IDs.
  """
  column_name = "event_id"

  def __init__(self, n = 0):
    self.n = n

  def new(self, row):
    self.n += 1
    a = self.n // 100000
    b = self.n % 100000
    return lsctables.SnglInspiralID(a * 1000000000 + row.get_id_parts()[1] * 100000 + b)


def ReadSnglInspiralFromFiles(fileList, verbose=False, filterFunc=None):
  """
  Read the SnglInspiralTables from a list of files.
  If filterFunc is not None, only keep triggers for which filterFunc
  evaluates to True.  Ex.: filterFunc=lambda sng: sng.snr >= 6.0

  @param fileList: list of input files
  @param verbose: print progress
  """
  # NOTE: this function no longer carries out event ID mangling (AKA
  # reassignment). Please adjust calling codes accordingly!
  # This means that identical event IDs produced by lalapps_thinca in
  # non-slide files having the same GPS start time will stay identical,
  # affecting zerolag and injection runs made over the same data.
  #
  # In consequence, if the calling code is going to reconstruct coincs
  # from the sngl event IDs, and if these include multiple injection
  # runs, coinc finding should be done one file at a time - see the
  # readCoincInspiralFromFiles function in CoincInspiralUtils.py

  sngls = lsctables.New(lsctables.SnglInspiralTable, \
      columns=lsctables.SnglInspiralTable.loadcolumns)

  lsctables.use_in(ExtractSnglInspiralTableLIGOLWContentHandler)
  for i,file in enumerate(fileList):
    if verbose: print str(i+1)+"/"+str(len(fileList))+": "
    xmldoc = utils.load_filename(file, verbose=verbose, contenthandler=ExtractSnglInspiralTableLIGOLWContentHandler)
    try:
      sngl_table = table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
      if filterFunc is not None:
        iterutils.inplace_filter(filterFunc, sngl_table)
    except ValueError: #some xml files have no sngl table, that's OK
      sngl_table = None
    if sngl_table: sngls.extend(sngl_table)
    xmldoc.unlink()    #free memory

  return sngls


def ReadSnglInspiralSlidesFromFiles(fileList, shiftVector, vetoFile=None,
  verbose=False):
  """
  Function for reading time-slided single inspiral triggers
  with automatic resliding the times, given a list of input files.

  @param fileList: List of files containing single inspiral time-slided
                   triggers
  @param shiftVector: Dictionary of time shifts to apply to triggers
                      keyed by IFO
  @param vetoFile: segwizard formatted file used to veto all triggers
  @param verbose: print progress
  """
  # NOTE: This function also will not mangle (reassign) lalapps_thinca
  # style event IDs (see the previous function).
  # Hence it will fail if fed event ID-related options.

  # read raw triggers
  inspTriggers = ReadSnglInspiralFromFiles(fileList, verbose=verbose)
  if inspTriggers:
    # get the rings
    segDict = SearchSummaryUtils.GetSegListFromSearchSummaries(fileList)
    rings = segments.segmentlist(iterutils.flatten(segDict.values()))
    rings.sort()

    # perform the veto
    if vetoFile is not None:
      segList = segmentsUtils.fromsegwizard(open(vetoFile))
      inspTriggers = inspTriggers.veto(segList)

    # now slide all the triggers within their appropriate ring
    slideTriggersOnRingWithVector(inspTriggers, shiftVector, rings)

  # return the re-slided triggers
  return inspTriggers


def ReadSnglInspiralsForPipelineStage(xmldoc, slideDict, stage):
  """
  Collect the sngl_inspiral rows from a desired stage in the pipeline.
  -- For the INSPIRAL and zerolag THINCA stages, the entire sngl_inspiral table is returned.
  -- For the slide THINCA stages, return only the rows in the sngl_inspiral that
     compose a coincident event from the desired time-slide
  @param xmldoc:    ligolw_xml doc
  @param slideDict: dictionary of the desired time-slide (eg. {u'H1': 0.0, u'L1': 100.0})
  @param stage:    the name of the stage (INSPIRAL_FIRST, THINCA_0_CAT_2, etc)
  """

  sngls_tbl = table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
  if 'THINCA' in stage and slideDict:
    # get the time-slides as a dictionary
    time_slide_tbl = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName)
    time_slide_dict = time_slide_tbl.as_dict()

    coinc_event_tbl = table.get_table(xmldoc, lsctables.CoincTable.tableName)
    coinc_map_tbl = table.get_table(xmldoc, lsctables.CoincMapTable.tableName)

    # get the time_slide_ids that have 
    time_slide_ids = set()
    for tsid, offsets in time_slide_dict.items():
      if offsets.values() == slideDict.values():
        time_slide_ids.add( tsid )
    # get the coinc_event_ids for coincs from the desired slides
    coinc_event_ids = set()
    for coinc in coinc_event_tbl:
      if coinc.time_slide_id in time_slide_ids:
        coinc_event_ids.add( coinc.coinc_event_id )
    # get the inspiral event_ids associated with the above coinc events
    event_ids = set()
    for row in coinc_map_tbl:
      if row.coinc_event_id in coinc_event_ids:
        event_ids.add( row.event_id )

    sngls_tbl_eid = sngls_tbl.getColumnByName("event_id")
    coinc_sngls_tbl = xmldoc.childNodes[0].insertBefore( lsctables.New(lsctables.SnglInspiralTable), sngls_tbl)
    for idx, event_id in enumerate(event_ids):
      coinc_sngls_tbl.insert( idx, sngls_tbl[sngls_tbl_eid.index(event_id)] )

    return coinc_sngls_tbl
  else:
    return sngls_tbl


#
# =============================================================================
#
#                                  Clustering
#
# =============================================================================
#

def CompareSnglInspiralByEndTime(a, b):
  """
  Orders a and b by peak time.
  """
  return cmp(a.get_end(), b.get_end())


def CompareSnglInspiralBySnr(a, b):
  """
  Orders a and b by peak time.
  """
  return cmp(a.snr, b.snr)


def CompareSnglInspiral(a, b, twindow = LIGOTimeGPS(0)):
  """
  Returns 0 if a and b are less than twindow appart.
  """
  tdiff = abs(a.get_end() - b.get_end())
  if tdiff < twindow:
    return 0
  else:
    return cmp(a.get_end(), b.get_end())

#
# =============================================================================
#
#                              Time Sliding on Ring
#
# =============================================================================
#

def slideTriggersOnLines(triggerList, shifts):
  """
  In-place modify trigger_list so that triggers are slid by appropriate value
  of shifts.

  @param triggerList: a SnglInspiralTable
  @param shifts:      a dictionary of the time-shifts keyed by IFO
  """
  for trigger in triggerList:
    end_time = trigger.get_end()
    trigger.set_end( end_time + shifts[trigger.ifo] )

def slideTimeOnRing(time, shift, ring):
  """
  Return time after adding shift, but constrained to lie along the ring
  """
  assert time in ring
  # use ring[0] as an epoch, do arithmetic using floats relative to epoch
  return ring[0] + (float(time - ring[0]) + shift) % float(abs(ring))


def slideTriggersOnRings(triggerList, rings, shifts):
  """
   In-place modify trigger_list so that triggers are slid by appropriate value
   of shifts along their enclosing ring segment by the algorithm given in XXX.
   This function calls the function slideTimeOnRing

   @param triggerList: a SnglInspiralTable
   @param rings:       sorted segment list of possible rings
   @param shifts:      a dictionary of the time-shifts keyed by IFO
  """
  for trigger in triggerList:
    end_time = trigger.get_end()
    trigger.set_end(slideTimeOnRing(end_time, shifts[trigger.ifo], rings[rings.find(end_time)]))

def unslideTriggersOnRings(triggerList, rings, shifts):
  """
   In-place modify trigger_list so that triggers are unslid by appropriate
   value of shifts along their enclosing ring segment by the algorithm given in
   XXX.
   This function calls the function slideTriggersOnRing

   @param triggerList: a SnglInspiralTable
   @param rings:       sorted segment list of possible rings
   @param shifts:      a dictionary of the time-shifts keyed by IFO
  """
  slideTriggersOnRings(triggerList, rings, dict((ifo, -shift) for ifo, shift in shifts.items()))

def slideTriggersOnRingWithVector(triggerList, shiftVector, rings):
   """
   In-place modify trigger_list so that triggers are slid by
   along their enclosing ring segment by the algorithm given in XXX.
   Slide numbers are extracted from the event_id of each trigger,
   and multiplied by the corresponding (ifo-keyed) entry in shift_vector
   to get the total slide amount.
   This function is called by ReadSnglInspiralSlidesFromFiles and
   calls the function slideTimeOnRing

   @param triggerList: a SnglInspiralTable
   @param shiftVector: a dictionary of the unit time-shift vector,
                       keyed by IFO
   @param rings:       sorted segment list of possible rings
   """
   for trigger in triggerList:
     end_time = trigger.get_end()
     trigger.set_end(slideTimeOnRing(end_time, trigger.get_slide_number() * shiftVector[trigger.ifo], rings[rings.find(end_time)]))

def slideSegListDictOnRing(ring, seglistdict, shifts):
  """
   Return seglistdict with segments that are slid by appropriate values of
   shifts along the ring segment by the algorithm given in XXX.

   @param ring:        segment on which to cyclicly slide segments in
                       seglistdict
   @param seglistdict: segments to be slid on ring
   @param shifts:      a dictionary of the time-shifts keyed by IFO
  """
  # don't do this in loops
  ring_duration = float(abs(ring))

  # automate multi-list arithmetic
  ring = segments.segmentlistdict.fromkeys(seglistdict.keys(), segments.segmentlist([ring]))

  # make a copy, and extract the segments that are in the ring.
  seglistdict = seglistdict & ring

  # apply the shift vector.  the shift vector is first normalized so that
  # all the shifts are non-negative and less than the duration of the ring.
  seglistdict.offsets.update(dict((instrument, shift % ring_duration) for instrument, shift in shifts.items()))

  # split the result into the pieces that are still in the ring, and the
  # pieces that have fallen off the edge.  both retain the shift vector in
  # the offsets attribute.
  extra = seglistdict - ring
  seglistdict &= ring

  # wrap the fallen-off pieces around the ring.
  for instrument in extra.keys():
    extra.offsets[instrument] -= ring_duration

  # return the union of the results.  the shifts vector is retained in the
  # offsets attribute of the result
  return seglistdict | extra


def compute_thinca_livetime(on_instruments, off_instruments, rings, vetoseglistdict, offsetvectors):
  """
  @on_instruments is an iterable of the instruments that must be on.

  @off_instruments is an iterable of the instruments that must be off.

  on_instruments and off_instruments must be disjoint.

  @rings is a list of segments defining the analysis ring boundaries.  They
  can overlap, and do not need to be ordered.

  @vetoseglistdict is a coalesced glue.segments.segmentlistdict object
  providing the veto segments for whatever instruments have vetoes defined
  for them.  This can include veto lists for instruments other than those
  listed in on_ and off_instruments (extra veto lists will be ignored), and
  it need not provide lists for all instruments (instruments for which
  there are no veto segment lists are assumed to be on at all times).

  @offsetvectors is an iterable of dictionaries of instrument-offset pairs.
  Each dictionary must contain entries for all instruments in the union of
  on_instruments and off_instruments (it is allowed to name others as well,
  but they will be ignored).  An example of one dictionary of
  instrument-offset pairs:  {"H1": 0.0, "H2": 5.0, "L1": 10.0}.

  The return value is a float giving the livetime in seconds.
  """
  # local copies so they can be modified and iterated over more than once
  # (in case generator expressions have been passed in)
  on_instruments = set(on_instruments)
  off_instruments = set(off_instruments)

  # check that the on and off instruments are disjoint
  if on_instruments & off_instruments:
    raise ValueError, "on_instruments and off_instruments not disjoint"

  # instruments that are not vetoed are assumed to be on
  on_instruments &= set(vetoseglistdict.keys())

  # performance aid:  only need offsets for instruments whose state is
  # important
  all_instruments = on_instruments | off_instruments
  offsetvectors = tuple(dict((key, value) for key, value in offsetvector.items() if key in all_instruments) for offsetvector in offsetvectors)

  # performance aid:  if there are no offset vectors to consider, the
  # livetime is trivial
  if not offsetvectors:
    return []

  # check that each offset vector provides values for all instruments of
  # interest
  for offsetvector in offsetvectors:
    if not set(offsetvector.keys()).issuperset(all_instruments):
      raise ValueError, "incomplete offset vector %s;  missing instrument(s) %s" % (repr(offsetvector), ", ".join(all_instruments - set(offsetvector.keys())))

  # initialize the livetime sums
  live_time = [0.0] * len(offsetvectors)

  # the livetime is trivial if an instrument that must be off is never
  # vetoed
  if not set(vetoseglistdict.keys()).issuperset(off_instruments):
    return live_time

  # performance aid:  don't need veto segment lists for instruments whose
  # state is unimportant, nor veto segments that don't intersect the rings
  coalesced_rings = segments.segmentlist(rings).coalesce()
  vetoseglistdict = segments.segmentlistdict((key, segments.segmentlist(seg for seg in seglist if coalesced_rings.intersects_segment(seg))) for key, seglist in vetoseglistdict.items() if key in all_instruments)

  # tot up the time when exactly the instruments that must be on are on
  for ring in rings:
    # don't do this in loops
    ring = segments.segmentlist([ring])

    # performance aid:  this is done in the loop, inside
    # slideSegListDictOnRing(), but we can make that go faster by doing it
    # here first
    clipped_vetoseglistdict = segments.segmentlistdict((key, seglist & ring) for key, seglist in vetoseglistdict.items())

    # performance aid:  if an instrument that must be vetoed is never
    # vetoed in this ring, the livetime is zero
    if not all(clipped_vetoseglistdict[key] for key in off_instruments):
      continue

    # iterate over offset vectors
    for n, offsetvector in enumerate(offsetvectors):
      # apply the offset vector to the vetoes, wrapping around the ring
      slidvetoes = slideSegListDictOnRing(ring[0], clipped_vetoseglistdict, offsetvector)

      # slidvetoes = times when instruments are vetoed,
      # slidvetoes.union(on_instruments) = times when an instrument that
      # must be on is vetoed
      #
      # ~slidvetoes = times when instruments are not vetoed,
      # (~slidvetoes).union(off_instruments) = times when an instrument
      # that must be off is not vetoed
      live_time[n] += float(abs(ring - slidvetoes.union(on_instruments) - (~slidvetoes).union(off_instruments)))

  # done
  return live_time
