#!/usr/bin/env python

# Author: Chris Pankow

# TODO:
# 1. Allow for multiple ROOT files to be placed on the command line
# 2. Find a way to include other tables which require the IDs to be linked
# (i. e. coinc_events and such).
# 3. Convert coherent information (netcc, rho, etc) into a new table

import sys
from optparse import OptionParser

try:
  # ligo_lw xml handling modules
  import glue

  from glue import segments
  from glue.ligolw import ligolw
  from glue.ligolw import lsctables
  from glue.ligolw import table
  from glue.ligolw import types
  from glue.ligolw import utils
  from glue.ligolw import ilwd

  from glue.lal import LIGOTimeGPS
  from pylal import llwapp
except:
  sys.exit("ERROR: Unable to import modules from GLUE.")

try:
  from ROOT import TFile, TChain, TTree
  from ROOT import gDirectory
except:
  sys.exit("ERROR: Unable to import modules from ROOT.")

def branch_array_to_list(branch, len):
  list = []
  for i in range(0, len):
    list.append( branch[i] )

  return list

# Index map copied from detector definitions in WaveBurst
# Constructs a string with all relevant detector strings
# frm a set of indexes given
def get_ifos_from_index(indx):
  ifostr = ""
  if not isinstance(indx, list):
    indx=[indx]

  for i in indx :
    if i == 1 : ifostr = "%s%s" % (ifostr, "L1")
    if i == 2 : ifostr = "%s%s" % (ifostr, "H1")
    if i == 3 : ifostr = "%s%s" % (ifostr, "H2")
    if i == 4 : ifostr = "%s%s" % (ifostr, "G1")
    if i == 5 : ifostr = "%s%s" % (ifostr, "T1")
    if i == 6 : ifostr = "%s%s" % (ifostr, "V1")
    if i == 7 : ifostr = "%s%s" % (ifostr, "A1")
	
  return ifostr

# Sets up table structures and calls populating methods
def create_tables(xmldoc, rootfiles):

  sim_tree = TChain("waveburst")
  liv_tree = TChain("liveTime")
  for rootfile in rootfiles :
    sim_tree.Add(rootfile)
    liv_tree.Add(rootfile)

  # Define tables
  sngl_burst_table = lsctables.New(lsctables.SnglBurstTable,
	["peak_time_ns", "start_time_ns", "stop_time_ns",
        "process_id", "ifo", "peak_time", "start_time", "stop_time", 
	"duration", "time_lag", "peak_frequency", "search",
	"flow", "fhigh", "bandwidth", "tfvolume", "hrss", "event_id"])
  xmldoc.childNodes[0].appendChild(sngl_burst_table)
  sngl_burst_table.sync_next_id()

  coinc_event_table = lsctables.New(lsctables.CoincTable,
	[ "process_id", "coinc_event_id", "nevents", "instruments", "time_slide_id",
	"coinc_def_id"])
  xmldoc.childNodes[0].appendChild(coinc_event_table)
  coinc_event_table.sync_next_id()

  multi_burst_table = lsctables.New(lsctables.MultiBurstTable,
	["process_id", "peak_time", "peak_time_ns", "coinc_event_id"])
  xmldoc.childNodes[0].appendChild(multi_burst_table)

  coinc_event_map_table = lsctables.New(lsctables.CoincMapTable)
  xmldoc.childNodes[0].appendChild(coinc_event_map_table)

  do_process_table(xmldoc, sim_tree, liv_tree)
  process_index = dict((int(row.process_id), row) for row in table.get_table(xmldoc, lsctables.ProcessTable.tableName))

  do_summary_table(xmldoc, sim_tree, liv_tree)

  # create coinc_definer row
  row = get_coinc_def_row(sim_tree)
  coinc_def_id = llwapp.get_coinc_def_id(xmldoc, row.search, row.search_coinc_type, description = row.description)

  for i in range(0, sim_tree.GetEntries()):
    sim_tree.GetEntry(i)

    offset_vector = dict((get_ifos_from_index(instrument_index), offset) for instrument_index, offset in zip(sim_tree.ifo, sim_tree.lag))

    coinc_event = coinc_event_table.RowType()
    coinc_event.process_id = process_index[sim_tree.run].process_id
    coinc_event.coinc_event_id = coinc_event_table.get_next_id()
    coinc_event.coinc_def_id = coinc_def_id
    coinc_event.nevents = sim_tree.ndim
    coinc_event.instruments = get_ifos_from_index( branch_array_to_list ( sim_tree.ifo, sim_tree.ndim ) )
    coinc_event.time_slide_id = llwapp.get_time_slide_id(xmldoc, offset_vector, process_index[sim_tree.run])
    coinc_event_table.append(coinc_event)

    for d in range(0, sim_tree.ndim):
      sngl_burst = get_sngl_burst_row(sngl_burst_table, sim_tree, d)
      sngl_burst.process_id = coinc_event.process_id
      sngl_burst.event_id = sngl_burst_table.get_next_id()
      sngl_burst_table.append(sngl_burst)

      coinc_event_map = coinc_event_map_table.RowType()
      coinc_event_map.event_id = sngl_burst.event_id
      coinc_event_map.table_name = sngl_burst.event_id.table_name
      coinc_event_map.coinc_event_id = coinc_event.coinc_event_id
      coinc_event_map_table.append(coinc_event_map)

    multi_burst = get_multi_burst_row(multi_burst_table, sim_tree)
    multi_burst.process_id = coinc_event.process_id
    multi_burst.coinc_event_id = coinc_event.coinc_event_id
    multi_burst_table.append(multi_burst)

def get_coinc_def_row(sim_tree):
  return lsctables.CoincDef(search = u"waveburst", search_coinc_type = 0, description = u"fully coherent search (burst events)")

def get_sngl_burst_row(sngl_burst_table, sim_tree, d):
  row = sngl_burst_table.RowType()
  setattr(row, "search", "waveburst")
  # Interferometer name -> ifo
  setattr(row, "ifo", get_ifos_from_index(sim_tree.ifo[d]) )
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
  setattr(row, "time_lag", sim_tree.lag[d])
  # Frequency
  # Central frequency in the detector -> frequency
  setattr(row, "peak_frequency", sim_tree.frequency[d])
  # Low frequency of the event in the detector -> flow
  setattr(row, "flow", sim_tree.low[d])
  # High frequency of the event in the detector ->  fhigh
  setattr(row, "fhigh", sim_tree.high[d])
  # Bandwidth
  setattr(row, "bandwidth", sim_tree.bandwidth[d])
  # Shape
  # number of pizels on the TF plane -> tfsize
  setattr(row, "tfvolume", sim_tree.size[d])
  # Energy
  # energy / noise variance -> snr
  setattr(row, "snr", sim_tree.snr[d])
  # TODO: What to do with this? GW strain
  #setattr(row, "strain", sim_tree.strain[d])
  # h _ root square sum
  setattr(row, "hrss", sim_tree.hrss[d])
  
  return row

def get_multi_burst_row(multi_burst_table, sim_tree):
  row = multi_burst_table.RowType()
  row.set_peak(sum(LIGOTimeGPS(t) for t in sim_tree.time) / sim_tree.ndim)
  return row

def do_process_table(xmldoc, sim_tree, liv_tree):
  try: 
    process_table = table.get_table(xmldoc, lsctables.ProcessTable.tableName)
  except ValueError:
    process_table = lsctables.New(lsctables.ProcessTable,
    ["process_id", "ifos", "comment", "program", "start_time", "jobid", 
		"end_time"])
    xmldoc.childNodes[0].appendChild(process_table)

  runids = set()
  for i in range(0, sim_tree.GetEntries()) :
    sim_tree.GetEntry(i)

    # Id for the run processed by WaveBurst -> process ID
    if sim_tree.run in runids :
      continue

    row = process_table.RowType()
    row.process_id = type(process_table.next_id)(sim_tree.run)
    runids.add(sim_tree.run)

    # Imstruments involved in the search
    row.ifos = lsctables.ifos_from_instrument_set( get_ifos_from_index( branch_array_to_list ( sim_tree.ifo, sim_tree.ndim ) ) )
    row.comment = "waveburst"
    row.program = "waveburst"

    # Begin and end time of the segment
    # TODO: This is a typical offset on either side of the job for artifacts
    # It can, and probably will change in the future, and should not be hardcoded
    setattr(row, "start_time", None)
    setattr(row, "end_time", None)
    setattr(row, "jobid", sim_tree.run)

    process_table.append(row)

def do_summary_table(xmldoc, sim_tree, liv_tree):
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
  for i in range(0, sim_tree.GetEntries()) :
    sim_tree.GetEntry(i)

    # Id for the run processed by WaveBurst -> process ID
    if sim_tree.run in runids :
      continue

    row = search_summary.RowType()
    row.process_id = process_id_type(sim_tree.run)
    runids.add(sim_tree.run)

    # Search Summary Table
    # events found in the run -> nevents
    setattr(row, "nevents", sim_tree.GetEntries())

    # Imstruments involved in the search
    row.ifos = lsctables.ifos_from_instrument_set( get_ifos_from_index( branch_array_to_list ( sim_tree.ifo, sim_tree.ndim ) ) )
    setattr(row, "comment", "waveburst")

    # Begin and end time of the segment
    # TODO: This is a typical offset on either side of the job for artifacts
    # It can, and probably will change in the future, and should not be hardcoded
		# TODO: Make this work properly. We need a gps end from the livetime
    waveoffset = 8
    livetime = 600
    #live_entries = live_tree.GetEntries()
    # This is WAAAAAAAAAAAAAY too slow
    #for l in range(0, live_entries):
      #liv_tree.GetEntry(l)
      #livetime = max(livetime, liv_tree.live)

    #if livetime < 0:
      #sys.exit("Could not find livetime, cannot fill all of summary table.")
    # in -- with waveoffset
    # out -- without waveoffset
    row.set_in(segments.segment(LIGOTimeGPS(sim_tree.gps - waveoffset), LIGOTimeGPS(sim_tree.gps + livetime + waveoffset)))
    row.set_out(segments.segment(LIGOTimeGPS(sim_tree.gps), LIGOTimeGPS(sim_tree.gps + livetime)))

    search_summary.append(row)


# Begin main driver script

parser = OptionParser()
parser.add_option("-o", "--output-name", action="store", dest="outputname", 
  help="Name of XML file to output.")
#parser.add_option("-s", "--segment-file", action="store", dest="segfile", 
  #help="Name of segment file to integrate.")

opts, rootfiles = parser.parse_args()

if rootfiles == None:
  sys.exit("Must specify input Coherent WaveBurst ROOT file with -r option")

xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW());

create_tables(xmldoc, rootfiles)

if opts.outputname == None :
  print "Assigning name to output xml"
  output = "convertROOT.xml.gz"
else :
  output = "%s.xml.gz" % opts.outputname

utils.write_filename(xmldoc, output, verbose = True, 
  gz = (output or "stdout").endswith(".gz"))

xmldoc.unlink()
