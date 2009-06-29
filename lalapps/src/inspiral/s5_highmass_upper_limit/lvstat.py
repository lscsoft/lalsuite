#!/usr/bin/python
import scipy
from scipy import interpolate
import numpy
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
from pylal.xlal.date import LIGOTimeGPS

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

def get_injections(injfnames, zero_lag_segments, ifos="H1,H2,L1", FAR=1.0, verbose = True):
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
    print >>sys.stderr, "getting injections: " + str(FAR) + ":\t%.1f%%\r" % (100.0 * cnt / len(injfnames),),
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
      AND coinc_inspiral.ifos == ?
  )
FROM
  sim_inspiral
WHERE
  -- only interested in injections that were injected
  injection_was_made(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
    """, (FAR,ifos,)):
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

class LVstat(object):
  """
  LVstat stores the data necessary for computing the LVstat and
  updates the databases
  """
  def __init__(self):
    self.eff_factor = {}

  def calc_eff_factor(self, on_ifos, trigger_ifos, m, f):
    """
    The eff_factor depends on trigger type and on ifos
    it is the distance cubed weighted efficiency
    it might not be right for injections that are not log
    distributed in distance
    """
    key = (on_ifos, trigger_ifos)
    
    M = 0.0
    F = 0.0
    for inj in m:
      M += inj.distance**3
    for inj in f:
      F += inj.distance**3
    self.eff_factor[key] = F / (M+F)


  def lvstat(self, on_ifos, trigger_ifos, false_alarm_rate):
    """
    lvstat is a function created in the database to compute the effective likelihood
    """

    key = (on_ifos, trigger_ifos)
    return  numpy.log( self.eff_factor[key] / (false_alarm_rate * 31536000)  + 1.0 )

  def update_coincs(self, fnames):
    """
    This function iterates over the databases and updates the the likelihood 
    column with the proper lv stat
    """
    for f in fnames:
      working_filename = dbtables.get_connection_filename(zerofname, verbose = True)
      connection = sqlite3.connect(working_filename)
      dbtables.DBTable_set_connection(connection)
      connection.create_function("lvstat", 3, self.lvstat)
      query = "UPDATE coinc_event SET likelihood = (SELECT lvstat(coinc_event.instruments, coinc_inspiral.ifos, coinc_inspiral.false_alarm_rate) FROM coinc_inspiral WHERE coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id)"
      print query

      connection.cursor().execute(query)
      connection.commit()
      connection.close()
      dbtables.discard_connection_filename(zerofname, working_filename, verbose = True)
      dbtables.DBTable_set_connection(None)

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

# A class to compute the LVstat
lvs = LVstat()

# FIXME  get the on instruments in the zero lag files,
# this assumes that the instruments are same for injection runs
# they probably should be
for on_ifos in get_on_instruments(sys.argv[1]):
  #get the correct segment lists and remove vetos
  instruments = lsctables.instrument_set_from_ifos(on_ifos[0])
  FAR, segs = get_far_threshold_and_segments(zerofname, on_ifos[0], "thinca")
  segs -= vetosegs
  zero_lag_segments = segs.intersection(instruments) - segs.union(set(segs.keys()) - instruments)

  for trigger_ifos in get_ifos(sys.argv[1]):
    #get the ifos that participate in coincidence
    key = (on_ifos[0], trigger_ifos[0])
    f, m = get_injections(sys.argv[2:], zero_lag_segments, trigger_ifos[0])
    lvs.calc_eff_factor(on_ifos[0], trigger_ifos[0], m, f)
    print "Eff factor: ", lvs.eff_factor[key], " key: ", key

# update the coincs with the effective likelihood
lvs.update_coincs(sys.argv)
