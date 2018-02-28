#!/usr/bin/python
#
# Copyright (C) 2011 Chris Pankow
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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

# TODO:
# 2. Find a way to include other tables which require the IDs to be linked
# (i. e. coinc_events and such).

"""
Convert a ROOT file produced by cWB into a ligolw document.
"""

import sys
import os
import platform
import re

# ligo_lw xml handling modules
import glue

from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw import types
from glue.ligolw import utils
from glue.ligolw import ilwd
from glue.ligolw.utils import time_slide

from glue.lal import LIGOTimeGPS
from glue.lal import CacheEntry
from lalburst import timeslides as ligolw_tisi
from pylal import llwapp
from pylal import git_version

# TODO: Reimplement this when the cwb_table module is included
#try:
  # Import cwb module
  #import cwb_table
#except:
  #sys.exit("ERROR: Unable to import modules from cwb_table.")

try:
  from ROOT import TFile, TChain, TTree
  from ROOT import gDirectory
except ImportError:
  print >>sys.stderr, "WARNING: Unable to import modules from ROOT. ROOT file processing will not be possible."

#CONVERT_ROOT_VERSION = "0.0.8"
__author__ = "Chris Pankow <chris.pankow@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

#
# =============================================================================
#
#                                  Helper Functions
#
# =============================================================================
#

def branch_array_to_list(branch, len):
  """
  Turn a ROOT TBranch into a python list
  """

  list = []
  for i in range(0, len):
    list.append( branch[i] )

  return list

def get_ifos_from_index(indx):
  """
  Index map copied from detector definitions in WaveBurst. Constructs a string with all relevant detector strings from a set of indexes given.
  """

  if not isinstance(indx, list):
    indx = [indx]

  ifolist = []
  for i in indx :
    if i == 1 : ifolist.append("L1")
    if i == 2 : ifolist.append("H1")
    if i == 3 : ifolist.append("H2")
    if i == 4 : ifolist.append("G1")
    if i == 5 : ifolist.append("T1")
    if i == 6 : ifolist.append("V1")
    if i == 7 : ifolist.append("A1")

  if len(ifolist) == 0 : return None 
  if len(ifolist) == 1 : return ifolist[0] 
  return ifolist

#
# =============================================================================
#
#                                  XML Table Functions
#
# =============================================================================
#

class CWB2Coinc(object):
  """
  Class to convert a set of rootfiles to a ligolw_document.
  """

  def __init__(self, joblist=None, start=None, end=None, instruments=None,
          #cwbtables=False,
          waveoffset=0, verbose=False):
	self.job_list = joblist
	self.start = start
	self.end = end
	#self.cwbtable = cwbtables
	self.verbose = verbose
	self.instruments = lsctables.ifos_from_instrument_set( instruments.split(",") )
	self.waveoffset = waveoffset

  def create_tables(self, xmldoc, rootfiles):
    """
    Sets up table structures and calls populating methods.
    """

    if os.path.splitext(rootfiles[0])[1] == ".root":
      sim_tree = TChain("waveburst")
    else: # If the file is (for example) text, use a proxy class
      sim_tree = CWBTextConverter()

    for rootfile in rootfiles:
	  sim_tree.Add(rootfile)

    # Define tables
    sngl_burst_table = lsctables.New(lsctables.SnglBurstTable,
	  ["peak_time_ns", "start_time_ns", "stop_time_ns",
	  "process_id", "ifo", "peak_time", "start_time", "stop_time", 
	  "duration", "time_lag", "peak_frequency", "search",
	  "flow", "fhigh", "bandwidth", "tfvolume", "hrss", "event_id", "snr"])
    xmldoc.childNodes[0].appendChild(sngl_burst_table)
    sngl_burst_table.sync_next_id()
  
    coinc_event_table = lsctables.New(lsctables.CoincTable,
	  [ "process_id", "coinc_event_id", "nevents", "instruments", "time_slide_id",
	  "coinc_def_id", "likelihood"])
    xmldoc.childNodes[0].appendChild(coinc_event_table)
    coinc_event_table.sync_next_id()
  
    # TODO: Reimplement this when the cwb_table module is included
    #if self.cwbtable:
      #cohwb_table = lsctables.New(cwb_table.CoherentWaveburstTable,
	    #["ellipticity", "correlated_energy", "eff_correlated_energy", 
      #"coherent_network_amp", "penalty", "network_correlation", 
      #"energy_disbalance", "ellip_energy_disbalance", "process_id", 
      #"coinc_event_id", "cwb_id"])
      #xmldoc.childNodes[0].appendChild(cohwb_table)
      #cohwb_table.sync_next_id()
  
    multi_burst_table = lsctables.New(lsctables.MultiBurstTable,
	  ["process_id", "peak_time", "peak_time_ns", "coinc_event_id", "snr", "ifos",
    # NOTE: Added to the table definition
    "false_alarm_rate", "ligo_axis_ra", "ligo_axis_dec"])
    xmldoc.childNodes[0].appendChild(multi_burst_table)
  
    coinc_event_map_table = lsctables.New(lsctables.CoincMapTable)
    xmldoc.childNodes[0].appendChild(coinc_event_map_table)
  
    jobsegment = None
    if self.start and self.end :
      jobsegment = segments.segment(LIGOTimeGPS(self.start), LIGOTimeGPS(self.end))
  
    if self.verbose:
      print "Creating Process Table...",
  
    if self.job_list:
      self.do_process_table(xmldoc, sim_tree)
    else :
      self.do_process_table_from_segment(xmldoc, sim_tree, jobsegment)
  
    process_index = dict((int(row.process_id), row) for row in table.get_table(xmldoc, lsctables.ProcessTable.tableName))
  
    if self.verbose:
      print " done."
  
    if self.verbose:
      print "Creating Summary Table...",
    # If we are provided a joblist, use it to generate the list
    if self.job_list :
      self.do_summary_table_from_joblist(xmldoc, sim_tree)
    elif self.job_list == None and self.start and self.end :
      self.do_summary_table_from_segment(xmldoc, jobsegment, sim_tree)
    else :
      self.do_summary_table(xmldoc, sim_tree)
    if self.verbose:
      print " done."
  
    # create coinc_definer row
    row = self.get_coinc_def_row(sim_tree)
    coinc_def_id = llwapp.get_coinc_def_id(xmldoc, row.search, row.search_coinc_type, description = row.description)
  
    entries = sim_tree.GetEntries()
    for i in range(0, entries):
      sim_tree.GetEntry(i)
      if self.start != None :
        if float(self.start) > sim_tree.start[0] : continue
      if self.end != None :
        if float(self.end) < sim_tree.stop[0] : continue
        
  
      offset_vector = dict((get_ifos_from_index(instrument_index), offset) for instrument_index, offset in zip(sim_tree.ifo, sim_tree.lag))
  
      coinc_event = coinc_event_table.RowType()
      coinc_event.process_id = process_index[sim_tree.run].process_id
      coinc_event.coinc_event_id = coinc_event_table.get_next_id()
      coinc_event.coinc_def_id = coinc_def_id
      coinc_event.nevents = sim_tree.ndim
      coinc_event.instruments = lsctables.ifos_from_instrument_set( get_ifos_from_index( branch_array_to_list ( sim_tree.ifo, sim_tree.ndim ) ) )
      coinc_event.time_slide_id = time_slide.get_time_slide_id(xmldoc, offset_vector, process_index[sim_tree.run])
      coinc_event.likelihood = sim_tree.likelihood
      coinc_event_table.append(coinc_event)
  
      for d in range(0, sim_tree.ndim):
        sngl_burst = self.get_sngl_burst_row(sngl_burst_table, sim_tree, d)
        sngl_burst.process_id = coinc_event.process_id
        sngl_burst.event_id = sngl_burst_table.get_next_id()
        sngl_burst_table.append(sngl_burst)
  
        coinc_event_map = coinc_event_map_table.RowType()
        coinc_event_map.event_id = sngl_burst.event_id
        coinc_event_map.table_name = sngl_burst.event_id.table_name
        coinc_event_map.coinc_event_id = coinc_event.coinc_event_id
        coinc_event_map_table.append(coinc_event_map)
  
      # TODO: Reimplement when cwb_table module is included
      #if self.cwbtable:
        #cwb = get_cwb_row(cohwb_table, sim_tree)
        #cwb.process_id = coinc_event.process_id
        #cwb.coinc_event_id = coinc_event.coinc_event_id
        #cwb.cwb_id = cohwb_table.get_next_id()
        #cohwb_table.append(cwb)
  
      multi_burst = self.get_multi_burst_row(multi_burst_table, sim_tree)
      multi_burst.process_id = coinc_event.process_id
      multi_burst.coinc_event_id = coinc_event.coinc_event_id
      # NOTE: Until we have an embedded cwb table definition, this will be
      # copied here so official tools have a ranking statistic to use
      try:
        multi_burst.snr = sim_tree.rho[1]
      except TypeError: # difference in definition between ROOT and text
        multi_burst.snr = sim_tree.rho
      # NOTE: To be filled in later by farburst
      multi_burst.false_alarm_rate = -1.0 
      # Reconstructed right ascension and declination
      multi_burst.ligo_axis_ra = sim_tree.phi[2]
      multi_burst.ligo_axis_dec = sim_tree.theta[2]
      multi_burst.ifos = lsctables.ifos_from_instrument_set( get_ifos_from_index( branch_array_to_list ( sim_tree.ifo, sim_tree.ndim ) ) )
  
      multi_burst_table.append(multi_burst)
  
  def get_coinc_def_row(self, sim_tree):
    """
    Fill in a coinc definer entry for waveburst
    """
  
    return lsctables.CoincDef(search = u"waveburst", search_coinc_type = 0, description = u"fully coherent search (burst events)")
  
  def get_cwb_row(self, cohwb_table, sim_tree):
    """
    Fill in a custom cwb event table
    """
  
    row = cohwb_table.RowType()
    # ellipticity of reconstructed event
    row.ellipticity = sim_tree.norm
    # sum of correlated energy
    row.correlated_energy = sim_tree.ecor
    # sum of correlated energy weighted by likelihood matrix
    row.eff_correlated_energy = sim_tree.ECOR
	  # Coherent Network Amplitude
	  # The currently used definition uses effective correlated energy which
	  # is the second entry in the array (the first is normal corrrelated energy)
    row.coherent_network_amp = sim_tree.rho[1]
	  # Penalty factor
	  # is an expression of the mismatch between signal energy in the h(t)
	  # and the energy in the reconstructed signal, if the recon. signal energy is
	  # greater than the data signal energy -- otherwise this is identically unity.
    row.penalty = sim_tree.penalty
	  # Network Correlation Coefficient
	  # Correlated energy divided by sum of null and correlated energy
	  # Note again, the first entry is used, corresponding to ecor rather than 
	  # eff_ecor in this case
    row.network_correlation = sim_tree.netcc[0]
	  # Energy Disbalance
	  # Summed difference of the reconstructed and data signal energies
	  # See Penaly factor for description of above
    row.energy_disbalance = sim_tree.neted[0]/sim_tree.ecor
	  # ELliptical Energy Disbalance
	  # Same as above, but between the signal and phase shifted signal
    row.ellip_energy_disbalance = sim_tree.neted[1]/sim_tree.ecor
    
    return row
  
  def get_sngl_burst_row(self, sngl_burst_table, sim_tree, d):
    """
    Fill in a sngl_burst row for a cwb event 
    """
    row = sngl_burst_table.RowType()
    row.search = u"waveburst"
    # Interferometer name -> ifo
    row.ifo = get_ifos_from_index(sim_tree.ifo[d] )
    # Timing
    peak = LIGOTimeGPS(sim_tree.time[d])
    seg = segments.segment(LIGOTimeGPS(sim_tree.start[d]), LIGOTimeGPS(sim_tree.stop[d]))
    # Central time in the detector -> cent_time
    row.set_peak(peak)
    # Start time in the detector -> end_time
    row.set_start(seg[0])
    # Stop time in the detector -> star_time
    row.set_stop(seg[1])
    # Event duration
    row.duration = abs(seg)
    # TODO: Make sure this is right = Time lag used to shift detector -> lag
    row.time_lag = sim_tree.lag[d]
    # Frequency
    # Central frequency in the detector -> frequency
    row.peak_frequency = sim_tree.frequency[d]
    try:
        # Low frequency of the event in the detector -> flow
        row.flow = sim_tree.low[d]
        # High frequency of the event in the detector ->  fhigh
        row.fhigh = sim_tree.high[d]
    except TypeError:
        row.flow = sim_tree.low
        row.fhigh = sim_tree.high
    # Bandwidth
    row.bandwidth = sim_tree.bandwidth[d]
    # Shape
    # number of pizels on the TF plane -> tfsize
    row.tfvolume = sim_tree.size[d]
    # Energy
    # energy / noise variance -> snr
    row.snr = sim_tree.snr[d]**(1./2.)
    # TODO: What to do with this? GW strain
    #row.strain = sim_tree.strain[d])
    # h _ root square sum
    row.hrss = sim_tree.hrss[d]
    
    return row
  
  def get_multi_burst_row(self, multi_burst_table, sim_tree):
    """
    Fill in a multi_burst row for a cwb event  -- these should exactly correlate to the sngl_burst events within a single coherent event
    """
    row = multi_burst_table.RowType()
	  # TODO : Script chokes on the line below: we need a better way to handle this anyway
    #row.set_peak(sum(LIGOTimeGPS(t) for t in sim_tree.time[0:3]) / int(sim_tree.ndim))
    row.set_peak(LIGOTimeGPS(sim_tree.time[0]))
    return row
  
  def do_process_table_from_segment(self, xmldoc, sim_tree, segment, jobid=-1):
    """
    Create the process_table for the cWB job(s) from a job segment.
    """
  
    try: 
      process_table = table.get_table(xmldoc, lsctables.ProcessTable.tableName)
    except ValueError:
      process_table = lsctables.New(lsctables.ProcessTable,
      ["process_id", "ifos", "comment", "program", "start_time", "jobid", 
		  "end_time",
      # Unused, but required for column synchronicity		
      "version", "cvs_repository", "cvs_entry_time", "is_online", "node",
      "username", "unix_procid", "domain"])
      xmldoc.childNodes[0].appendChild(process_table)

  
    # get run IDs
    if(jobid < 0):
      N = sim_tree.GetEntries()
      if N == 0:
        exit("There is no information available to parse from the ROOT file, and no external information was provided.")
      runs = set()
      for i in xrange(N):
        sim_tree.GetEntry(i)
        run = sim_tree.run
        runs.add(int(run))
    else:
      run = jobid
    seg = segment
  
    # Imstruments involved in the search
    ifos = lsctables.ifos_from_instrument_set( get_ifos_from_index( branch_array_to_list ( sim_tree.ifo, sim_tree.ndim ) ) )

    if( ifos == None or len(ifos) == 0 ): 
        if( self.instruments ):
            ifos = self.instruments
        else: # Not enough information to completely fill out the table
            sys.exit("Found a job with no IFOs on, or not enough to determine IFOs. Try specifying instruments directly.")

    for run in runs:
      row = process_table.RowType()
      row.process_id = type(process_table.next_id)(run)
      row.ifos = ifos
      row.comment = u"waveburst"
      row.program = u"waveburst"
      row.start_time = None
      row.end_time = None
      row.jobid = run
      row.version = __version__
      row.cvs_repository = u"lscsoft"
      # TODO: Need a handle on where this comes from
      row.cvs_entry_time = 0
      row.is_online = 0
      row.username = os.environ['USER']
      row.node = platform.node()
      row.unix_procid = os.getpid()
      row.domain = u"pylal"
      process_table.append(row)
  
  def do_process_table(self, xmldoc, sim_tree):
    """
    Create the process_table for the cWB job(s)
    """
  
    try: 
      process_table = table.get_table(xmldoc, lsctables.ProcessTable.tableName)
    except ValueError:
      process_table = lsctables.New(lsctables.ProcessTable,
      ["process_id", "ifos", "comment", "program", "start_time", "jobid", 
		  "end_time",
      # Unused, but required for column synchronicity		
      "version", "cvs_repository", "cvs_entry_time", "is_online", "node",
      "username", "unix_procid", "domain"])
      xmldoc.childNodes[0].appendChild(process_table)
  
    runids = set()
    runid = dict()
  
    itr = 0 
    for line in file(self.job_list) :
      if line[0] == '#' :
        continue # skip comments
  
      # Determine the file type
      # WARNING: Mixing the two types could cause overwrite behavior
      line = line.split()
      if len(line) == 2 :  # start and stop
        seg = segments.segment( map( LIGOTimeGPS, map(float, line[0:2]) ) )
        runid[itr] = seg
        itr += 1
      elif len(line) == 3 :  # index start and stop
        seg = segments.segment( map( LIGOTimeGPS, map(float, line[1:3]) ) )
        runid[int(line[0])] = seg
      elif len(line) == 4 :  # index start and stop and length
        seg = segments.segment( map( LIGOTimeGPS, map(float, line[1:3]) ) )
        runid[int(line[0])] = seg
      else : # dunno!
        sys.exit("Unable to understand job list segment format. Why do you confused me so?")
  
  
    sim_tree.GetEntry(0)
    for run, seg in runid.iteritems() :
      if self.start != None :
        if float(self.start) > seg[0] : continue
      if self.end != None :
        if float(self.end) < seg[1] : continue
  
      row = process_table.RowType()
      row.process_id = type(process_table.next_id)(run)
      #print "Run id %d was given process_id %d" % (sim_tree.run, row.process_id)
      runids.add(sim_tree.run)
  
      # Imstruments involved in the search
      ifos = lsctables.ifos_from_instrument_set( get_ifos_from_index( branch_array_to_list ( sim_tree.ifo, sim_tree.ndim ) ) )

      if( ifos == None or len(ifos) == 0 ): 
          if( self.instruments ):
              ifos = self.instruments
          else: # Not enough information to completely fill out the table
              sys.exit("Found a job with no IFOs on, or not enough to determine IFOs. Try specifying instruments directly.")
  
      row.ifos = ifos
      row.comment = u"waveburst"
      row.program = u"waveburst"
  
      row.start_time = None
      row.end_time = None
      row.jobid = run
  
      row.version = __version__
      row.cvs_repository = u"lscsoft"
      # TODO: Need a handle on where this comes from
      row.cvs_entry_time = 0
      row.is_online = 0
      row.username = os.environ['USER']
      row.node = platform.node()
      row.unix_procid = os.getpid()
      row.domain = "pylal"
  
      process_table.append(row)
  
  def do_summary_table_from_segment(self, xmldoc, segment, sim_tree, jobid=-1):
    """
    Create the search_summary table for a cWB from a segment specified from the command line. The function will try to determine the proper job intervals from the waveoffset, if specified.
    """
  
    try: 
      search_summary = table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName)
    except ValueError:
      search_summary = lsctables.New(lsctables.SearchSummaryTable,
      ["process_id", "nevents", "ifos", "comment", "in_start_time",
      "in_start_time_ns", "out_start_time", "out_start_time_ns",
      "in_end_time", "in_end_time_ns", "out_end_time", "out_end_time_ns"])
      xmldoc.childNodes[0].appendChild(search_summary)
  
    process_id_type = type(table.get_table(xmldoc, lsctables.ProcessTable.tableName).next_id)
  
    sim_tree.GetEntry(0)

    if(jobid < 0):
      run = sim_tree.run
    else: run = jobid
    seg = segment
  
    # Search Summary Table
    # events found in the run -> nevents
    row = search_summary.RowType()
    row.process_id = process_id_type(run)
    row.nevents = sim_tree.GetEntries()

    ifos = lsctables.ifos_from_instrument_set( get_ifos_from_index( branch_array_to_list ( sim_tree.ifo, sim_tree.ndim ) ) )
    # Imstruments involved in the search
    if( ifos == None or len(ifos) == 0 ): 
        if( self.instruments ):
            ifos = self.instruments
        else: # Not enough information to completely fill out the table
            sys.exit("Found a job with no IFOs on, or not enough to determine IFOs. Try specifying instruments directly.")

    row.ifos = ifos
    row.comment = "waveburst"
  
    # Begin and end time of the segment
    waveoffset = self.waveoffset
    if waveoffset == None: waveoffset = 0
  
    # in -- with waveoffset
    row.set_in(seg)
    # out -- without waveoffset
    waveoffset = LIGOTimeGPS(waveoffset)
    row.set_out(segments.segment(seg[0]+waveoffset, seg[1]-waveoffset))
    search_summary.append(row)
  
  def do_summary_table_from_joblist(self, xmldoc, sim_tree):
    """
    Create the search_summary table for the cWB job(s) and a provided cWB joblist. The function will try to determine the proper job intervals from the waveoffset, if specified.
    """
  
    try: 
      search_summary = table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName)
    except ValueError:
      search_summary = lsctables.New(lsctables.SearchSummaryTable,
      ["process_id", "nevents", "ifos", "comment", "in_start_time",
      "in_start_time_ns", "out_start_time", "out_start_time_ns",
      "in_end_time", "in_end_time_ns", "out_end_time", "out_end_time_ns"])
      xmldoc.childNodes[0].appendChild(search_summary)
  
    process_id_type = type(table.get_table(xmldoc, lsctables.ProcessTable.tableName).next_id)
  
    runid = dict()
    itr = 0 
    for line in file(self.job_list) :
      if line[0] == '#' :
        continue # skip comments
  
      # Determine the file type
      # WARNING: Mixing the two types could cause overwrite behavior
      line = line.split()
	  
      if len(line) == 2 :  # start and stop
        seg = segments.segment( map( LIGOTimeGPS, map(float, line[0:2]) ) )
        runid[itr] = seg
        itr += 1
      elif len(line) == 3 :  # index start and stop
        seg = segments.segment( map( LIGOTimeGPS, map(float, line[1:3]) ) )
        runid[int(line[0])] = seg
      elif len(line) == 4 :  # index start and stop and length
        seg = segments.segment( map( LIGOTimeGPS, map(float, line[1:3]) ) )
        runid[int(line[0])] = seg
      else : # dunno!
        sys.exit("Unable to understand job list segment format.")
  
    for run, seg in runid.iteritems() :
      if self.start != None :
        if float(self.start) > seg[0] : continue
      if self.end != None :
        if float(self.end) < seg[1] : continue
  
      row = search_summary.RowType()
      row.process_id = process_id_type(run)
  
      # Search Summary Table
      # events found in the run -> nevents
      row.nevents = 0
      #entries = sim_tree.GetEntries()
  
      # Imstruments involved in the search
      sim_tree.GetEntry(0)
      ifos = lsctables.ifos_from_instrument_set( get_ifos_from_index( branch_array_to_list ( sim_tree.ifo, sim_tree.ndim ) ) )
      # Imstruments involved in the search
      if( ifos == None or len(ifos) == 0 ): 
          if( self.instruments ):
              ifos = self.instruments
          else: # Not enough information to completely fill out the table
              sys.exit("Found a job with no IFOs on, or not enough to determine IFOs. Try specifying instruments directly.")

      row.ifos = ifos
      row.comment = "waveburst"
  
      # Begin and end time of the segment
      # TODO: This is a typical offset on either side of the job for artifacts
      # It can, and probably will change in the future, and shouldn't be hardcoded
      #waveoffset, livetime = 8, -1
      waveoffset, livetime = self.waveoffset, -1
      if waveoffset == None: waveoffset = 0
      livetime = abs(seg)
  
      if livetime < 0:
        print >>sys.stderr, "WARNING: Run %d will have zero livetime because no events are recorded in this run, and therefore livetime cannot be calculated." % run
        # in -- with waveoffset
        row.set_in(segments.segment(LIGOTimeGPS(0), LIGOTimeGPS(1)))
        # out -- without waveoffset
        row.set_out(segments.segment(LIGOTimeGPS(0), LIGOTimeGPS(1)))
      else:
        # in -- with waveoffset
        row.set_in(seg)
        # out -- without waveoffset
        row.set_out(segments.segment(seg[0]+waveoffset, seg[1]-waveoffset))
  
      search_summary.append(row)
  
  def do_summary_table(self, xmldoc, sim_tree):
    """
    Create the search_summary table for the cWB job(s). This function exists as a backup in case no job list exists. It will try and reconstruct the job segments from the event data, but this list will be incomplete in the case where no events were found during a job.
    """
  
    try: 
      search_summary = table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName)
    except ValueError:
      search_summary = lsctables.New(lsctables.SearchSummaryTable,
      ["process_id", "nevents", "ifos", "comment", "in_start_time",
      "in_start_time_ns", "out_start_time", "out_start_time_ns",
      "in_end_time", "in_end_time_ns", "out_end_time", "out_end_time_ns"])
      xmldoc.childNodes[0].appendChild(search_summary)
  
    process_id_type = type(table.get_table(xmldoc, lsctables.ProcessTable.tableName).next_id)
  
    runids = set()
    entries = sim_tree.GetEntries()
    for i in range(0, entries) :
      sim_tree.GetEntry(i)
  
      if self.start != None :
        if float(self.start) > sim_tree.start[0] : continue
      if self.end != None :
        if float(self.end) < sim_tree.stop[0] : continue
  
      # Id for the run processed by WaveBurst -> process ID
      run = sim_tree.run
      if run in runids :
        continue
  
      row = search_summary.RowType()
      row.process_id = process_id_type(run)
      runids.add(run)
  
      # Search Summary Table
      # events found in the run -> nevents
	  # TODO: Destroy ROOT, because it hates me
      #row.nevents = sim_tree.nevent
      row.nevents = 0
  
      # Imstruments involved in the search
      ifos = lsctables.ifos_from_instrument_set( get_ifos_from_index( branch_array_to_list ( sim_tree.ifo, sim_tree.ndim ) ) )
      # Imstruments involved in the search
      if( ifos == None or len(ifos) == 0 ): 
          if( self.instruments ):
              ifos = self.instruments
          else: # Not enough information to completely fill out the table
              sys.exit("Found a job with no IFOs on, or not enough to determine IFOs. Try specifying instruments directly.")

      row.ifos = ifos
      row.comment = "waveburst"
  
      # Begin and end time of the segment
      if self.waveoffset != None :
          waveoffset = self.waveoffset
      else: 
          waveoffset = 0
  
      livetime = -1
	  # WARNING: This will fail miserably if there are no events in analyzed 
	  # see do_search_summary_from_joblist for a better solution
      livetime = sim_tree.left[0] + sim_tree.duration[0] + sim_tree.right[0]
  
      if livetime < 0:
        print >>sys.stderr, "WARNING: Run %d will have zero livetime because no events are recorded in this run, and therefore livetime cannot be calculated." % run
        # in -- with waveoffset
        row.set_in(segments.segment(LIGOTimeGPS(0), LIGOTimeGPS(1)))
        # out -- without waveoffset
        row.set_out(segments.segment(LIGOTimeGPS(0), LIGOTimeGPS(1)))
      else:

        seg_start_with_offset = sim_tree.gps - waveoffset
        seg_start_without_offset = sim_tree.gps
        seg_end_with_offset = sim_tree.gps + waveoffset + livetime
        seg_end_without_offset = sim_tree.gps + livetime 
        # in -- with waveoffset
        row.set_in(segments.segment(LIGOTimeGPS(seg_start_with_offset), LIGOTimeGPS(seg_end_with_offset)))
        # out -- without waveoffset
        row.set_out(segments.segment(LIGOTimeGPS(seg_start_without_offset), LIGOTimeGPS(seg_end_without_offset)))
  
      search_summary.append(row)
  
class CWBTextConverter(object):
	"""
	Class to transparently emit data as if the input were a ROOT file. The input is an unstructured text file made by cWB at runtime.
	"""

	CWB_SEARCH_TYPES = { 'r': "unmodelled", 'i': "elliptical", 's': "linear", 'g': "circular", 'R': "unmodelled", 'I': "elliptical", 'S': "linear", 'G': "circular" }

	def __init__(self):
		self.data = []
		self.this_data = None
		self.version = None
		self.online_version = None
		self.search_type = None
		self._skymap_reg = re.compile( "^map_lenght" ) # NOT A TYPO, THIS IS ACTUALLY IN THE FILE.
		self._far_info_reg = re.compile( "^##" )

	def __getattr__(self, name):
		"""
		This overrides the attribute getter for this class so that when the proxy class is required to produce information based on an attribute (say 'rho') it responds just like a TBranch would. If no event is 'selected' via GetEntry, then None is returned, regardless of the attribute, because there's no way to know whether we're asking for a TBranch leaf or something else.
		"""
		if self.this_data is None:
			return None
		try:
			return self.this_data[name]
		except KeyError:
			return getattr(self, name)

	def GetEntry(self, i):
		"""
		Emulate a TTree GetEntry() call by switching internal pointer to a row in the table.
		"""
		self.this_data = self.data[i]

	def GetEntries(self):
		"""
		Emulate a TTree GetEntries() call by returning the length of the internal data table.
		"""
		return len(self.data)

	def Add(self, filen):
		"""
		Add a file which will be parsed by the class and added to the internal table.
		"""
		self.parse_transcript( open(filen) )

	def get_far_segments( self ):
		"""
		Gets the segments for which various FARs have been defined.
		"""
		if self.this_data["far_data"] is None:
			raise ValueError( "No FAR info for this event." )

		return self.this_data["far_data"].keys()


	def get_far_information( self, segment=None ):
		"""
		Returns the false alarm rate for a given segment, usually one of the segments retrieved from 'get_far_segments'. If segment is None (default), then the FAR from the biggest segment is returned.
		"""

		if self.this_data["far_data"] is None:
			raise ValueError( "No FAR info for this event." )

		# return FAR from largest segment of data
		if segment is None:
			segs = sorted( self.get_far_segments(), key=lambda s: abs(s) )
			return self.this_data["far_data"][segs[-1]]

		return self.this_data["far_data"][segment]


	def get_skymap_coords( self ):
		"""
		Get a list of the skymap coordinates (R.A. and dec) which cWB has probability values for.
		"""
		if self.this_data["skymap"] is None:
			raise ValueError( "No skymap set for this event." )

		return self.this_data["skymap"].keys()

	def get_skymap_prob( self, ra, dec ):
		"""
		Return the probability for a given R.A. and dec from the skymap data. If there is no entry for this, then 0 is returned. This is the implied value for no R.A. and dec entry in the skymap data.
		"""
		if self.this_data["skymap"] is None:
			raise ValueError( "No skymap set for this event." )

		try:
			return self.this_data["skymap"][(ra, dec)]
		except KeyError:
			return 0

	def parse_with_type(self, val):
		"""
		Attempt to coerce data to the proper type. Note that this depends on python's finickiness with int / float parsing. It won't do string -> float -> int truncation with one call to float(). So, this check whether or not int() works, then tries float(). Then it gives up and gives the value (presumed a string) back.
		"""
		try:
			return int(val)
		except ValueError:
			pass
		try:
			return float(val)
		except ValueError:
			pass
		return val

	def parse_transcript(self, text):
		"""
		Parses a textual event transcript and adds the data to the internal table.
		"""
		row, mode = {}, "var"
		row["skymap"] = {}
		row["far_data"] = {}

		for line in text.readlines():

			# ALl the unstructured stuff is at the end...
			if "wat version" in line:
				self.version = line.split("=")[1].strip()
				continue
			elif "online version" in line:
				self.online_version = line.split("=")[1].strip()
				continue
			elif "search" in line:
				self.search_type = self.CWB_SEARCH_TYPES[line.split("=")[1].strip()]
				continue
			elif re.match( self._far_info_reg, line ) is not None:
				mode = "far"
				continue
			elif re.match( self._skymap_reg, line ) is not None:
				mode = "skymap"
				continue
			elif len(line.strip()) == 0:
				# a blank line is assumed to be a new event entry
				row["gps"] = row["segment"][0]
				self.data.append( row )
				row = {}
				row["skymap"] = {}
				row["far_data"] = {}
				mode = "var"
				continue

			if line[0] == "#": continue

			# Modal parsing logic
			if mode == "var":

				varn, val = line.split(":")
				val = [ self.parse_with_type(v) for v in val.split() ]
				if len(val) == 1: val = val[0]

				row[varn] = val
			elif mode == "skymap":
				sid, theta, dec, step, phi, ra, step2, prob, cum = line.split()
				row["skymap"][(float(ra), float(dec))] = float(prob)
			elif mode == "far":
				try:
					rank, rate, segstart, segend, segdur = line.split()
				except ValueError:
					# There is now two lines at the end of the file which
					# contain extraneous paths and URLs
					continue
				far_seg = segments.segment( float(segstart), float(segend) )
				row["far_data"][far_seg] = float(rate)

		# Not currently in the file, but an alias is needed for compatibility
		row["gps"] = row["segment"][0]
		self.data.append( row )
		return row
