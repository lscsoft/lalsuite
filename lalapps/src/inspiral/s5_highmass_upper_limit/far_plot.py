import scipy
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab
from math import *
import sys
from glue.ligolw import ligolw
try:
        import sqlite3
except ImportError:
        # pre 2.5.x
        from pysqlite2 import dbapi2 as sqlite3
from glue.ligolw import dbtables

def fix_numbers(a, val=0.01):
  for i in range(len(a)):
    if a[i] <= 0: a[i] = val

def sigma_regions(m, s, nlist):
  out = {}
  for n in nlist:
    out[n] = numpy.concatenate((m-n*s, (m+n*s)[::-1]))
    fix_numbers(out[n])
  return out

def get_ifos(cursor):
  query = "SELECT distinct(instruments) FROM coinc_event;"
  return [ifo for ifo, in cursor.execute(query )]
  
def connect(fname):
  connection = sqlite3.connect(dbtables.get_connection_filename('FULL_DATACAT_3.sqlite', tmp_path = None, verbose = True))
  cursor = connection.cursor()
  return cursor

def extract_from_db(file, zero_ifar, slide_ifar):
  cursor = connect(sys.argv[1])
  ifos = get_ifos(cursor)
  for ifo in ifos:
    try: zero_ifar[ifo]
    except: zero_ifar[ifo] = []
    try: slide_ifar[ifo]
    except: slide_ifar[ifo] = []
    zero_ifar[ifo] += [ifar for ifar, in cursor.execute('SELECT 1.0 / coinc_inspiral.false_alarm_rate FROM coinc_inspiral JOIN coinc_event ON (coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id) WHERE (coinc_event.instruments == ? AND NOT EXISTS (SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset !=0 ) )ORDER BY coinc_inspiral.false_alarm_rate;', (ifo,))]
    slide_ifar[ifo] += [ifar for ifar, in cursor.execute('SELECT 1.0 / coinc_inspiral.false_alarm_rate FROM coinc_inspiral JOIN coinc_event ON (coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id) WHERE (coinc_event.instruments == ? AND EXISTS (SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset !=0 ) )ORDER BY coinc_inspiral.false_alarm_rate;', (ifo,) )]
  cursor.close()

zero_ifar = {}
slide_ifar = {}
for file in sys.argv[1:]:
  extract_from_db(file, zero_ifar, slide_ifar)

figcnt = 0
#FIXME fails if there is no background in a given ifo combo
for ifos in zero_ifar.keys():
  figcnt += 1
  pylab.figure(figcnt)
  pylab.loglog(zero_ifar[ifos], numpy.arange(len(zero_ifar[ifos])) + 1, linewidth=2)
  pylab.hold(1)
  #FIXME compensation should be ratio of live times
  m = (numpy.arange(len(slide_ifar[ifos])) + 1) / 100.00
  pylab.loglog(slide_ifar[ifos], m, 'k', linewidth=2)
  s = numpy.sqrt(m)
  x = numpy.array(slide_ifar[ifos])
  xc = numpy.concatenate((x, x[::-1]))
  y = sigma_regions(m, s, [1,2,3])
  pylab.fill(xc, y[1], alpha=0.25, facecolor=[0.25, 0.25, 0.25])
  pylab.fill(xc, y[2], alpha=0.25, facecolor=[0.5, 0.5, 0.5])
  pylab.fill(xc, y[3], alpha=0.25, facecolor=[0.75, 0.75, 0.75])
  pylab.xlim( x[-1], x[0] )
  pylab.ylim(0.01, len(zero_ifar[ifos]) )
  pylab.ylabel('Number')
  pylab.xlabel('IFAR based Rank')
  pylab.legend(("zero lag", "expected background, N", "1*sqrt(N)", "2*sqrt(N)","3*sqrt(N)"), loc="lower left")
  pylab.hold(0)
  pylab.savefig(sys.argv[1]+'_'+''.join(ifos.split(','))+"_far.png")
