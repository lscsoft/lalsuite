import sys 
from optparse import *
from glue import segments
from glue import segmentsUtils
from pylal import readMeta

def cleanlist(seglist, min_length):
  removals = segments.segmentlist()
  for seg in seglist:
    if seg.duration() < min_length:
      removals.append(seg)
  seglist = seglist - removals

  return seglist

# Chop a simInspiralTable into segments that we care about.
def getSegments ( sim_table, seglist, ifokey ):

  tmptable = readMeta.metaDataTable( None, "sim_inspiral" )

  for e in sim_table.table:
    end_time = e[ifokey]
    if seglist.__contains__(end_time):
      tmptable.table.append(e)

  return tmptable
# define usage and command line options and arguments - parse
usage = "usage: %prog ..."
parser = OptionParser( usage )

parser.add_option("-v", "--verbose", action="store_true",default=False,\
  help="make things verbose" )
parser.add_option("-H","--h1-segments",action="store",type="string",\
  default=None, metavar=" H1_SEGMENTS", help="H1 input segment to read" )
parser.add_option("-K","--h2-segments",action="store",type="string",\
  default=None, metavar=" H2_SEGMENTS", help="H2 input segment to read" )
parser.add_option("-L","--l1-segments",action="store",type="string",\
  default=None, metavar=" L1_SEGMENTS", help="L1 input segment to read" )
parser.add_option("-i","--injection-file",action="store",type="string",\
  default=None, metavar=" INJFILE", help="An injection file" )
parser.add_option("-g","--glitch-time",action="store",type="int",\
    default=None, metavar=" GLITCH_TIME",\
    help="produce plot of triggers around GLITCH_TIME")

( opts , args ) = parser.parse_args()

h1seglist = segmentsUtils.fromsegwizard(file(opts.h1_segments)).coalesce()
h1seglist = cleanlist(h1seglist, 2064)
h1seglist = h1seglist.contract(72)

h2seglist = segmentsUtils.fromsegwizard(file(opts.h2_segments)).coalesce()
h2seglist = cleanlist(h2seglist, 2064)
h2seglist = h2seglist.contract(72)

l1seglist = segmentsUtils.fromsegwizard(file(opts.l1_segments))
l1seglist = cleanlist(l1seglist, 2064)
l1seglist = l1seglist.contract(72)

triplelist = h1seglist & h2seglist & l1seglist
triplelist.coalesce()

h1h2doublelist = h1seglist & h2seglist
h1h2doublelist = h1h2doublelist.coalesce() - triplelist
h1l1doublelist = h1seglist & l1seglist
h1l1doublelist = h1l1doublelist - triplelist
h2l1doublelist = h2seglist & l1seglist
h2l1doublelist = h2l1doublelist - triplelist

if opts.injection_file:
  flist = [ opts.injection_file ]
  injections = readMeta.metaDataTable( flist, "sim_inspiral")
  print "No triple injections: %d" % \
      getSegments(injections, triplelist, "geocent_end_time").nevents()
  print "No H1H2 injections: %d" % \
      getSegments(injections, h1h2doublelist, "geocent_end_time").nevents()
  print "No H1L1 injections: %d" % \
      getSegments(injections, h1l1doublelist, "geocent_end_time").nevents()
  print "No H2L1 injections: %d" % \
      getSegments(injections, h2l1doublelist, "geocent_end_time").nevents()

if opts.glitch_time:
  if triplelist.__contains__(opts.glitch_time):
    print "Time " + str(opts.glitch_time) + " is in triple time"
  if h1h2doublelist.__contains__(opts.glitch_time):
    print "Time " + str(opts.glitch_time) + " is in h1h2 only time"
  if h1l1doublelist.__contains__(opts.glitch_time):
    print "Time " + str(opts.glitch_time) + " is in h1l1 only time"
  if h2l1doublelist.__contains__(opts.glitch_time):
    print "Time " + str(opts.glitch_time) + " is in h2l1 only time"
else:
  tmptime=triplelist.duration()
  print "Total triple time: %d s, %f yr" % (tmptime, tmptime/(365.25 * 24.0 * 3600.0))
  tmptime=h1h2doublelist.duration()
  print "Total H1H2-only time: %d, %f yr" % (tmptime, tmptime/(365.25 * 24.0 * 3600.0))
  tmptime=h1l1doublelist.duration()
  print "Total H1L1-only time: %d, %f yr" % (tmptime, tmptime/(365.25 * 24.0 * 3600.0))
  tmptime=h2l1doublelist.duration()
  print "Total H2L1-only time: %d, %f yr" % (tmptime, tmptime/(365.25 * 24.0 * 3600.0))
  tmptime=triplelist.duration()+ h1h2doublelist.duration() + h1l1doublelist.duration() + h2l1doublelist.duration()
  print "Total time: %d, %f yr" % (tmptime, tmptime/(365.25 * 24.0 * 3600.0))
