#!/usr/bin/python
import scipy
from scipy import interpolate
import numpy
import pylab

try:
        import sqlite3
except ImportError:
        # pre 2.5.x
        from pysqlite2 import dbapi2 as sqlite3
from math import *
import sys
import glob
import copy
from optparse import OptionParser

from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue.ligolw import utils
from pylal import db_thinca_rings
from pylal import llwapp
from pylal import rate
from pylal import SimInspiralUtils
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

lsctables.LIGOTimeGPS = LIGOTimeGPS

def get_far_threshold_and_segments(zerofname, instruments, live_time_program, veto_seg_name="vetoes", verbose = True):
  """
  return the false alarm rate of the most rare zero-lag coinc, and a
  dictionary of the thinca segments indexed by instrument.
  """
  # open database
  working_filename = dbtables.get_connection_filename(zerofname, verbose = verbose)
  connection = sqlite3.connect(working_filename)
  dbtables.DBTable_set_connection(connection)

  # extract false alarm rate threshold
  query = 'SELECT MIN(coinc_inspiral.false_alarm_rate) FROM coinc_inspiral JOIN coinc_event ON (coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id) WHERE ( coinc_event.instruments = "' + instruments + '" AND NOT EXISTS(SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset != 0) );'
  print "\n", query
  far, = connection.cursor().execute(query).fetchone()

  # extract segments.
  seglists = db_thinca_rings.get_thinca_zero_lag_segments(connection, program_name = live_time_program)

  # done
  connection.close()
  dbtables.discard_connection_filename(zerofname, working_filename, verbose = verbose)
  dbtables.DBTable_set_connection(None)

  return far, seglists

def get_injections(injfnames, zero_lag_segments, ifos="*", FAR=1000.0, verbose = True):
  """
  The LV stat uses simulations above threshold, not some IFAR of the loudest event, so the default should be "inf"
  """
  def injection_was_made(geocent_end_time, geocent_end_time_ns, zero_lag_segments = zero_lag_segments):
    """
    return True if injection was made in the given segmentlist
    """
    return lsctables.LIGOTimeGPS(geocent_end_time, geocent_end_time_ns) in zero_lag_segments

  found = []
  missed = []
  print >>sys.stderr, ""
  for cnt, f in enumerate(injfnames):
    print >>sys.stderr, "getting injections: " + ifos + ":\t%.1f%%\r" % (100.0 * cnt / len(injfnames),),
    working_filename = dbtables.get_connection_filename(f, tmp_path = None, verbose = verbose)
    connection = sqlite3.connect(working_filename)
    connection.create_function("injection_was_made", 2, injection_was_made)

    make_sim_inspiral = lsctables.table.get_table(dbtables.get_xml(connection), lsctables.SimInspiralTable.tableName)._row_from_cols

    # FIXME may not be done correctly if injections are done in timeslides
    # FIXME may not be done correctly if injections aren't logarithmic in d
    # do we really want d^3 waiting?
    # FIXME do we really want injections independent of their FAR

    for values in connection.cursor().execute("""
SELECT
  sim_inspiral.*,
  -- true if injection matched a coinc below the false alarm rate threshold
  EXISTS (
    SELECT
      *
    FROM
      coinc_event_map AS mapa
      JOIN coinc_event_map AS mapb ON (
        mapa.coinc_event_id == mapb.coinc_event_id
      )
      JOIN coinc_inspiral ON (
        mapb.table_name == "coinc_event"
        AND mapb.event_id == coinc_inspiral.coinc_event_id
      )
    WHERE
      mapa.table_name == "sim_inspiral"
      AND mapa.event_id == sim_inspiral.simulation_id
      AND coinc_inspiral.false_alarm_rate < ?
      --AND glob(?, coinc_inspiral.ifos)
  )
FROM
  sim_inspiral
WHERE
  -- only interested in injections that were injected
  injection_was_made(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
  ORDER BY sim_inspiral.eff_dist_h
    """, (FAR,)):
      sim = make_sim_inspiral(values)
      if values[-1]:
        found.append(sim)
      else:
        missed.append(sim)

    # done
    connection.close()
    dbtables.discard_connection_filename(f, working_filename, verbose = verbose)
    dbtables.DBTable_set_connection(None)

  print >>sys.stderr, "\nFound = %d Missed = %d" % (len(found), len(missed))
  return found, missed


def get_on_instruments(zerofname, verbose=True):
  # open database
  working_filename = dbtables.get_connection_filename(zerofname, verbose = verbose)
  connection = sqlite3.connect(working_filename)
  dbtables.DBTable_set_connection(connection)

  # extract false alarm rate threshold
  # FIXME This may not be robust if triggers are missing from a given category, for
  # example no triples in zero lag or time slides.
  query = 'SELECT distinct(instruments) FROM coinc_event'
  inst_list = []
  for i in connection.cursor().execute(query): inst_list.append(i)

  # done
  connection.close()
  dbtables.discard_connection_filename(zerofname, working_filename, verbose = verbose)
  dbtables.DBTable_set_connection(None)
  return inst_list

def get_ifos(zerofname, verbose=True):
  # open database
  working_filename = dbtables.get_connection_filename(zerofname, verbose = verbose)
  connection = sqlite3.connect(working_filename)
  dbtables.DBTable_set_connection(connection)

  # extract false alarm rate threshold
  # FIXME This may not be robust if triggers are missing from a given category, for
  # example no triples in zero lag or time slides.
  query = 'SELECT distinct(ifos) FROM coinc_inspiral'
  ifo_list = []
  for i in connection.cursor().execute(query): ifo_list.append(i)

  # done
  connection.close()
  dbtables.discard_connection_filename(zerofname, working_filename, verbose = verbose)
  dbtables.DBTable_set_connection(None)
  return ifo_list

def get_vetoes(fname, veto_segments_name = "vetoes", verbose=True):
  working_filename = dbtables.get_connection_filename(fname, verbose = verbose)
  connection = sqlite3.connect(working_filename)
  dbtables.DBTable_set_connection(connection)
  veto_segments = db_thinca_rings.get_veto_segments(connection, veto_segments_name)
  connection.close()
  dbtables.discard_connection_filename(fname, working_filename, verbose = verbose)
  dbtables.DBTable_set_connection(None)
  return veto_segments

###############################################################################
################### MAIN PROGRAM ##############################################
###############################################################################
# This program assigns lvstats to coinc table databases
#
# ./lvstat.py FULL_DATA.sqlite INJ1.sqlite INJ2.sqlite .....


# FIXME this assumes that the first argument is the zero lag / time slide database
zerofname = sys.argv[1]
# FIXME this assumes that subsequent files are injections
injfnames = sys.argv[2:]
vetosegs = get_vetoes(zerofname)
# FIXME  get the on instruments in the zero lag files,
# this assumes that the instruments are same for injection runs
# they probably should be

pylab.figure(1)
#FAR, segs = get_far_threshold_and_segments(zerofname, on_ifos[0], "thinca")
#segs -= vetosegs
#zero_lag_segments = segs.intersection(instruments) - segs.union(set(segs.keys()) - instruments)
#f, m = get_injections(sys.argv[2:], zero_lag_segments, trigger_ifos[0])

for on_ifos in get_on_instruments(sys.argv[1]):
  pylab.figure(1)
  #get the correct segment lists and remove vetos
  instruments = lsctables.instrument_set_from_ifos(on_ifos[0])
  FAR, segs = get_far_threshold_and_segments(zerofname, on_ifos[0], "thinca")
  segs -= vetosegs
  print on_ifos[0]
  zero_lag_segments = segs.intersection(instruments) - segs.union(set(segs.keys()) - instruments)
  #f, m = get_injections(sys.argv[2:], zero_lag_segments)
  #FIXME okay, really we need the decisive effective distance, this is right for the current search,
  # but will be wrong when including virgo...
  #pylab.semilogy([s.mass1+s.mass2 for s in m], [max([s.eff_dist_h, s.eff_dist_l]) for s in m], '.', color=[0,0,0])
  pylab.hold(1)
  legend = []
  for trigger_ifos in get_ifos(sys.argv[1]):
    #get the ifos that participate in coincidence
    key = (on_ifos[0], trigger_ifos[0])
    f, tmp = get_injections(sys.argv[2:], zero_lag_segments, trigger_ifos[0])
    #pylab.semilogy([s.mass1+s.mass2 for s in f], [max([s.eff_dist_h, s.eff_dist_l]) for s in f], '.')
    #pylab.semilogy([s.mchirp for s in f], [max([s.eff_dist_h, s.eff_dist_l]) for s in f], '.')
    #legend.append(trigger_ifos[0])
    #pylab.semilogy([s.mass1+s.mass2 for s in m], [s.eff_dist_h for s in m], '.k')
    #pylab.hold(1)
  f, m = get_injections(sys.argv[2:], zero_lag_segments)
  #FIXME okay, really we need the decisive effective distance, this is right for the current search,
  # but will be wrong when including virgo...
  #pylab.semilogy([s.mass1+s.mass2 for s in m], [max([s.eff_dist_h, s.eff_dist_l]) for s in m], '.', color=[0,0,0])
  pylab.semilogy([s.mchirp for s in m], [max([s.eff_dist_h, s.eff_dist_l]) for s in m], '.', color=[0,0,0])
  #pylab.semilogy([s.h_end_time for s in m], [max([s.eff_dist_h, s.eff_dist_l]) for s in m], '.', color=[0.5,0.5,0.5])

  legend.append('Missed')
  pylab.hold(0)
  pylab.legend(legend)
  pylab.xlabel('chirp mass')
  pylab.ylabel('max(Hanford effective distance, Livingston effective distance)')
  pylab.savefig(on_ifos[0].replace(",","") + "found_missed.png")
  pylab.clf()

  #print [s.h_end_time + s.h_end_time_ns/10**9 for s in m], [min([s.eff_dist_h, s.eff_dist_l]) for s in m]
