import scipy
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab
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
from pylal.xlal.date import LIGOTimeGPS

lsctables.LIGOTimeGPS = LIGOTimeGPS

def sigma_region(m, s, n):
  return numpy.concatenate((m-n*s, (m+n*s)[::-1]))

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

    if verbose:
      print >>sys.stderr, "computing background livetime:",
    for on_instruments, livetime in db_thinca_rings.get_thinca_livetimes(db_thinca_rings.get_thinca_rings_by_available_instruments(connection, program_name = live_time_program), veto_segments, db_thinca_rings.get_background_offset_vectors(connection), verbose = verbose).items():
      self.background_livetime[on_instruments] = self.background_livetime.get(on_instruments, 0.0) + float(livetime)
    if verbose:
      print >>sys.stderr

    #
    # compute the playground and non-playground zero-lag live time for each
    # choice of on instruments
    #

    zero_lag_segs = db_thinca_rings.get_thinca_zero_lag_segments(connection, program_name = live_time_program) - veto_segments
    playground_segs = segmentsUtils.S2playground(zero_lag_segs.extent_all())

    if verbose:
      print >>sys.stderr, "computing zero-lag livetime:",
    available_instruments = set(zero_lag_segs.keys())
    for on_instruments in (on_instruments for m in range(2, len(available_instruments) + 1) for on_instruments in iterutils.choices(sorted(available_instruments), m)):
      if verbose:
        print >>sys.stderr, ",".join(on_instruments),

      # compute times when only exactly on_instruments were on
      on_instruments = frozenset(on_instruments)
      segs = zero_lag_segs.intersection(on_instruments) - zero_lag_segs.union(available_instruments - on_instruments)

      self.playground_livetime[on_instruments] = self.playground_livetime.get(on_instruments, 0.0) + float(abs(segs & playground_segs))
      self.nonplayground_livetime[on_instruments] = self.nonplayground_livetime.get(on_instruments, 0.0) + float(abs(segs - playground_segs))
    if verbose:
      print >>sys.stderr

    #
    # collect the false alarm rates from each choice of on instruments
    #

    if verbose:
      print >>sys.stderr, "collecting false alarm rate ranks ..."

    for on_instruments in self.background_livetime:
      self.background_ifar.setdefault(on_instruments, [])
      self.playground_ifar.setdefault(on_instruments, [])
      self.nonplayground_ifar.setdefault(on_instruments, [])
    for ifar, end_time, end_time_ns, on_instruments, is_background in connection.cursor().execute("""
SELECT
  1.0 / coinc_inspiral.false_alarm_rate,
  coinc_inspiral.end_time,
  coinc_inspiral.end_time_ns,
  coinc_event.instruments,
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
    """):
      on_instruments = frozenset(lsctables.instrument_set_from_ifos(on_instruments))
      if is_background:
        self.background_ifar[on_instruments].append(ifar)
      elif lsctables.LIGOTimeGPS(end_time, end_time_ns) in playground_segs:
        self.playground_ifar[on_instruments].append(ifar)
      else:
        self.nonplayground_ifar[on_instruments].append(ifar)

    if verbose:
      print >>sys.stderr, "done"
    connection.close()
    dbtables.discard_connection_filename(filename, working_filename, verbose = verbose)


verbose = True

summary = Summary()
for file in sys.argv[1:]:
  summary.update_from_db(file, verbose = verbose)

if verbose:
  print >>sys.stderr, "generating plots:"
for figcnt, on_instruments in enumerate(summary.background_livetime):
  if verbose:
    print >>sys.stderr, "\t%s:" % ",".join(sorted(on_instruments))
  summary.background_ifar.setdefault(on_instruments, []).sort(reverse = True)
  summary.playground_ifar.setdefault(on_instruments, []).sort(reverse = True)
  summary.nonplayground_ifar.setdefault(on_instruments, []).sort(reverse = True)

  if verbose:
    print >>sys.stderr, "\t%f s background livetime, %f s non-playground zero-lag livetime (%f to 1)" % (summary.background_livetime[on_instruments], summary.nonplayground_livetime[on_instruments], summary.background_livetime[on_instruments] / summary.nonplayground_livetime[on_instruments])
  if not (summary.background_ifar[on_instruments] or summary.playground_ifar[on_instruments] or summary.nonplayground_ifar[on_instruments]):
    if verbose:
      print >>sys.stderr, "\tno events, skipping"
    continue

  pylab.figure(figcnt)
  N = numpy.arange(len(summary.nonplayground_ifar[on_instruments])) + 1
  if verbose:
    print >>sys.stderr, "\t%d non-playground zero-lag events" % max(N)
  minN, maxN = min(N), max(N)
  pylab.loglog(summary.nonplayground_ifar[on_instruments], N, 'k', linewidth=2)
  pylab.hold(1)
  N = numpy.arange(len(summary.background_ifar[on_instruments])) + 1
  if verbose:
    print >>sys.stderr, "\t%d background events" % max(N)
  N = N * summary.nonplayground_livetime[on_instruments] / summary.background_livetime[on_instruments]
  minN, maxN = min(minN, min(N)), max(maxN, max(N))
  pylab.loglog(summary.background_ifar[on_instruments], N, 'k--', linewidth=1)
  x = numpy.array(summary.background_ifar[on_instruments])
  minX, maxX = min(x), max(x)
  xc = numpy.concatenate((x, x[::-1]))
  pylab.fill(xc, sigma_region(N, numpy.sqrt(N), 1).clip(minN, maxN), alpha=0.25, facecolor=[0.25, 0.25, 0.25])
  pylab.fill(xc, sigma_region(N, numpy.sqrt(N), 2).clip(minN, maxN), alpha=0.25, facecolor=[0.5, 0.5, 0.5])
  pylab.fill(xc, sigma_region(N, numpy.sqrt(N), 3).clip(minN, maxN), alpha=0.25, facecolor=[0.75, 0.75, 0.75])
  pylab.xlim(minX, maxX)
  pylab.ylim(minN, maxN)
  pylab.ylabel('Number')
  pylab.xlabel('IFAR based Rank')
  pylab.legend(("zero lag", "expected background, N", r"$\pm\sqrt{N}$", r"$\pm 2\sqrt{N}$","$\pm 3\sqrt{N}$"), loc="lower left")
  pylab.hold(0)
  pylab.savefig(sys.argv[1]+'_'+''.join(sorted(on_instruments))+"_far.png")
  if verbose:
    print >>sys.stderr
