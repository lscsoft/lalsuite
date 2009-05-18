import scipy
import numpy
import matplotlib
matplotlib.rcParams.update({
  "font.size": 8.0,
  "axes.titlesize": 10.0,
  "axes.labelsize": 10.0,
  "xtick.labelsize": 8.0,
  "ytick.labelsize": 8.0,
  "legend.fontsize": 8.0,
  "figure.dpi": 600,
  "savefig.dpi": 600,
  "text.usetex": True
})
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import sys
try:
        import sqlite3
except ImportError:
        # pre 2.5.x
        from pysqlite2 import dbapi2 as sqlite3

from glue import segmentsUtils
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from pylal import db_thinca_rings
from pylal.xlal.date import LIGOTimeGPS

lsctables.LIGOTimeGPS = LIGOTimeGPS

def sigma_region(mean, nsigma):
  return numpy.concatenate((mean - nsigma * numpy.sqrt(mean), (mean + nsigma * numpy.sqrt(mean))[::-1]))

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
      self.background_livetime[on_instruments] = self.background_livetime.get(on_instruments, 0.0) + livetime
      self.background_ifar.setdefault(on_instruments, [])
      self.playground_ifar.setdefault(on_instruments, [])
      self.nonplayground_ifar.setdefault(on_instruments, [])
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
    for on_instruments in self.background_livetime:
      if verbose:
        print >>sys.stderr, ",".join(on_instruments),

      # compute times when only exactly on_instruments were on
      selected_segs = zero_lag_segs.intersection(on_instruments) - zero_lag_segs.union(available_instruments - on_instruments)

      self.playground_livetime[on_instruments] = self.playground_livetime.get(on_instruments, 0.0) + float(abs(selected_segs & playground_segs))
      self.nonplayground_livetime[on_instruments] = self.nonplayground_livetime.get(on_instruments, 0.0) + float(abs(selected_segs - playground_segs))
    if verbose:
      print >>sys.stderr

    #
    # collect the false alarm rates
    #

    if verbose:
      print >>sys.stderr, "collecting false alarm rate ranks ..."

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
for on_instruments in summary.background_livetime:
  if verbose:
    print >>sys.stderr, "\t%s:" % ",".join(sorted(on_instruments))
  summary.background_ifar[on_instruments].sort(reverse = True)
  summary.playground_ifar[on_instruments].sort(reverse = True)
  summary.nonplayground_ifar[on_instruments].sort(reverse = True)

  if verbose:
    print >>sys.stderr, "\t%f s background livetime, %f s non-playground zero-lag livetime (%f to 1)" % (summary.background_livetime[on_instruments], summary.nonplayground_livetime[on_instruments], summary.background_livetime[on_instruments] / summary.nonplayground_livetime[on_instruments])
  if not (summary.background_ifar[on_instruments] or summary.playground_ifar[on_instruments] or summary.nonplayground_ifar[on_instruments]):
    if verbose:
      print >>sys.stderr, "\tno events, skipping"
    continue

  fig = figure.Figure()
  FigureCanvas(fig)
  axes = fig.gca()
  axes.loglog()
  N = numpy.arange(len(summary.nonplayground_ifar[on_instruments]), dtype = "double") + 1.0
  if verbose:
    print >>sys.stderr, "\t%d non-playground zero-lag events" % max(N)
  x = numpy.array(summary.nonplayground_ifar[on_instruments], dtype = "double")
  axes.plot(x.repeat(2)[1:], N.repeat(2)[:-1], 'k', linewidth=2)
  minX, maxX = min(x), max(x)
  minN, maxN = min(N), max(N)
  N = numpy.arange(len(summary.background_ifar[on_instruments]), dtype = "double") + 1.0
  if verbose:
    print >>sys.stderr, "\t%d background events" % max(N)
  N *= summary.nonplayground_livetime[on_instruments] / summary.background_livetime[on_instruments]
  x = numpy.array(summary.background_ifar[on_instruments], dtype = "double")
  axes.plot(x.repeat(2)[1:], N.repeat(2)[:-1], 'k--', linewidth=1)
  minX, maxX = min(minX, min(x)), max(maxX, max(x))
  minN, maxN = min(minN, min(N)), max(maxN, max(N))
  x = x.repeat(2)[1:]
  x = numpy.concatenate((x, x[::-1]))
  N = N.repeat(2)[:-1]
  axes.fill(x, sigma_region(N, 3.0).clip(minN, maxN), alpha=0.25, facecolor=[0.75, 0.75, 0.75])
  axes.fill(x, sigma_region(N, 2.0).clip(minN, maxN), alpha=0.25, facecolor=[0.5, 0.5, 0.5])
  axes.fill(x, sigma_region(N, 1.0).clip(minN, maxN), alpha=0.25, facecolor=[0.25, 0.25, 0.25])
  axes.set_xlim(minX, maxX)
  axes.set_ylim(minN, maxN)
  axes.set_title('Zero-lag Events Observed Compared to Background')
  axes.set_xlabel('IFAR based Rank')
  axes.set_ylabel(r'Number of Events $\geq$ Rank')
  axes.legend(("zero lag", "expected background, N", r"$\pm 3\sqrt{N}$", r"$\pm 2\sqrt{N}$","$\pm\sqrt{N}$"), loc="lower left")
  fig.savefig(sys.argv[1]+'_'+''.join(sorted(on_instruments))+"_far.png")
  if verbose:
    print >>sys.stderr
