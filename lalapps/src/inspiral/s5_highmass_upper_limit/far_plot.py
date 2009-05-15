import scipy
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab
from math import *
import sys

try:
        import sqlite3
except ImportError:
        # pre 2.5.x
        from pysqlite2 import dbapi2 as sqlite3

from glue import iterutils
from glue import segmentsUtils
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from pylal import db_thinca_rings

def fix_numbers(a, val=0.01):
  for i in range(len(a)):
    if a[i] <= 0: a[i] = val

def sigma_regions(m, s, nlist):
  out = {}
  for n in nlist:
    out[n] = numpy.concatenate((m-n*s, (m+n*s)[::-1]))
    fix_numbers(out[n])
  return out

class Summary(object):
  def __init__(self):
    self.background_livetime = {}
    self.playground_livetime = {}
    self.nonplayground_livetime = {}
    self.playground_ifar = {}
    self.nonplayground_ifar = {}
    self.background_ifar = {}

  def update_from_db(self, filename, live_time_program = "thinca", vetoes_name = "vetoes", verbose = False):
    working_filename = dbtables.get_connection_filename(filename, tmp_path = None, verbose = verbose)
    connection = sqlite3.connect(working_filename)

    #
    # extract veto segments
    #

    veto_segments = db_thinca_rings.get_veto_segments(connection, name = vetoes_name)

    #
    # compute the total background live time for each choice of on
    # instruments
    #

    for on_instruments, livetime in db_thinca_rings.get_thinca_livetimes(db_thinca_rings.get_thinca_rings_by_available_instruments(connection, program_name = live_time_program), veto_segments, db_thinca_rings.get_background_offset_vectors(connection), verbose = verbose).items():
      try:
        self.background_livetime[on_instruments] += livetime
      except KeyError:
        self.background_livetime[on_instruments] = livetime

    #
    # compute the playground and non-playground zero-lag live time for each
    # choice of on instruments
    #

    zero_lag_segs = db_thinca_rings.get_thinca_zero_lag_segments(connection, program_name = live_time_program) - veto_segments
    playground_segs = segmentUtils.S2playground(zero_lag_segs.extent_all())

    available_instruments = set(zero_lag_segs.keys())
    for on_instruments in (on_instruments for m in range(2, len(available_instruments) + 1) for on_instruments in iterutils.choices(sorted(available_instruments), m)):
      on_instruments = frozeset(on_instruments)

      # compute times when only exactly on_instruments were on
      segs = zero_lag_segs.intersection(on_instruments) - zero_lag_segs.union(available_instruments - on_instruments)
      playground_livetime = abs(segs & playground_segs)
      nonplayground_livetime = abs(seg - playground_segs)

      try:
        self.playground_livetime[on_instruments] += playground_livetime
        self.nonplayground_livetime[on_instruments] += nonplayground_livetime
      except KeyError:
        self.playground_livetime[on_instruments] = playground_livetime
        self.nonplayground_livetime[on_instruments] = nonplayground_livetime

    #
    # collect the false alarm rates from each choice of on instruments
    #

    for on_instruments in (on_instruments for on_instruments, in connection.cursor().execute("SELECT DISTINCT(instruments) FROM coinc_event;")):
      key = frozenset(lsctables.instrument_set_from_ifos(on_instruments))
      if key not in self.background_ifar:
        self.playground_ifar[key] = []
        self.nonplayground_ifar[key] = []
        self.background_ifar[key] = []

      for ifar, end_time, end_time_ns, is_background in connection.cursor().execute("""
SELECT
  1.0 / coinc_inspiral.false_alarm_rate,
  coinc_inspiral.end_time,
  coinc_inspiral.end_time_ns,
  EXISTS (
    SELECT
      *
    FROM
      time_slide
    WHERE
      time_slide.time_slide_id == coinc_event.time_slide_id
      AND time_slide.offset != 0
    )
FROM
  coinc_inspiral
  JOIN coinc_event ON (
    coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id
  )
WHERE
  coinc_event.instruments == ?
ORDER BY
  coinc_inspiral.false_alarm_rate
      """, (on_instruments,)):
        if is_background:
          self.background_ifar[key].append(ifar)
	elif lsctables.LIGOTimeGPS(end_time, end_time_ns) in playground_segs:
          self.playground_ifar[key].append(ifar)
        else:
          self.nonplayground_ifar[key].append(ifar)

    connection.close()
    dbtables.discard_connection_filename(filename, working_filename, verbose = verbose)


summary = Summary()
for file in sys.argv[1:]:
  summary.update_from_db(file, verbose = True)


figcnt = 0
#FIXME fails if there is no background in a given ifo combo
for on_instruments in summary.nonplayground_ifar.keys():
  figcnt += 1
  pylab.figure(figcnt)
  pylab.loglog(summary.nonplayground_ifar[on_instruments], numpy.arange(len(summary.nonplayground_ifar[on_instruments])) + 1, linewidth=2)
  pylab.hold(1)
  #FIXME compensation should be ratio of live times
  m = (numpy.arange(len(summary.background_ifar[on_instruments])) + 1) * (summary.playground_livetime[on_instruments] / summary.background_livetime[on_instruments])
  pylab.loglog(summary.background_ifar[on_instruments], m, 'k', linewidth=2)
  s = numpy.sqrt(m)
  x = numpy.array(summary.background_ifar[on_instruments])
  xc = numpy.concatenate((x, x[::-1]))
  y = sigma_regions(m, s, [1,2,3])
  pylab.fill(xc, y[1], alpha=0.25, facecolor=[0.25, 0.25, 0.25])
  pylab.fill(xc, y[2], alpha=0.25, facecolor=[0.5, 0.5, 0.5])
  pylab.fill(xc, y[3], alpha=0.25, facecolor=[0.75, 0.75, 0.75])
  pylab.xlim( x[-1], x[0] )
  pylab.ylim(0.01, len(summary.zero_ifar[on_instruments]) )
  pylab.ylabel('Number')
  pylab.xlabel('IFAR based Rank')
  pylab.legend(("zero lag", "expected background, N", "1*sqrt(N)", "2*sqrt(N)","3*sqrt(N)"), loc="lower left")
  pylab.hold(0)
  pylab.savefig(sys.argv[1]+'_'+''.join(sorted(on_instruments))+"_far.png")
