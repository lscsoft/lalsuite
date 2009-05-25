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


def get_far_threshold_and_segments(zerofname, live_time_program, verbose = False):
  """
  return the false alarm rate of the most rare zero-lag coinc, and a
  dictionary of the thinca segments indexed by instrument.
  """
  # open database
  working_filename = dbtables.get_connection_filename(zerofname, verbose = verbose)
  connection = sqlite3.connect(working_filename)
  dbtables.DBTable_set_connection(connection)

  # extract false alarm rate threshold
  query = 'SELECT MIN(coinc_inspiral.false_alarm_rate) FROM coinc_inspiral JOIN coinc_event ON (coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id) WHERE NOT EXISTS(SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset != 0);'
  far, = connection.cursor().execute(query).fetchone()

  # extract segments.
  seglists = db_thinca_rings.get_thinca_zero_lag_segments(connection, program_name = live_time_program)

  # done
  connection.close()
  dbtables.discard_connection_filename(zerofname, working_filename, verbose = verbose)
  dbtables.DBTable_set_connection(None)
  return far, seglists

def get_volume_derivative(injfnames, twoDMassBins, dBin, FAR, zero_lag_segments, gw):
  livetime = float(abs(zero_lag_segments))
  FARh = FAR*100000
  FARl = FAR*0.001
  nbins = 5
  FARS = rate.LogarithmicBins(FARl, FARh, nbins)
  #gw = rate.gaussian_window2d(3,3,10)
  vA = []
  vA2 = []
  for far in FARS.centres():
    m, f = get_injections(injfnames, far, zero_lag_segments)
    print >>sys.stderr, "computing volume at FAR " + str(far)
    vAt, vA2t = twoD_SearchVolume(f, m, twoDMassBins, dBin, gw, livetime, 1)  
    vAt.array = scipy.log10(vAt.array + 0.001)
    vA.append(vAt)
  # the derivitive is calcuated with respect to FAR * t
  FARS = rate.LogarithmicBins(FARl * livetime, FARh * livetime, nbins)
  return derivitave_fit(FARS, FAR * livetime, vA, twoDMassBins)
  
  
def derivitave_fit(farts, FARt, vAs, twodbin):
  '''
     Relies on scipy spline fits for each mass bin
     to find the derivitave of the volume at a given
     FAR.  See how this works for a simple case where
     I am clearly giving it a parabola.  To high precision it calculates
     the proper derivitave. 
     A = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]
     B = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
     C = interpolate.splrep(B,A,s=0, k=4)
     interpolate.splev(5,C,der=1) 
     10.000
  '''
  dA = rate.BinnedArray(twodbin)
  for m1 in range(dA.array.shape[0]):
    for m2 in range(dA.array.shape[1]):
      da = []
      for f in farts.centres():
        da.append(vAs[farts[f]].array[m1][m2])
      fit = interpolate.splrep(farts.centres(),da,k=4)
      val = 0.0 - interpolate.splev(FARt,fit,der=1)
      dA.array[m1][m2] = val # minus the derivitave
  return dA

def get_injections(injfnames, FAR, zero_lag_segments, verbose = False):
  """
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
    print >>sys.stderr, "getting injections below FAR: " + str(FAR) + ":\t%.1f%%\r" % (100.0 * cnt / len(injfnames),),
    working_filename = dbtables.get_connection_filename(f, tmp_path = None, verbose = verbose)
    connection = sqlite3.connect(working_filename)
    connection.create_function("injection_was_made", 2, injection_was_made)

    make_sim_inspiral = lsctables.table.get_table(dbtables.get_xml(connection), lsctables.SimInspiralTable.tableName)._row_from_cols

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
  )
FROM
  sim_inspiral
WHERE
  -- only interested in injections that were injected
  injection_was_made(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
    """, (FAR,)):
      sim = make_sim_inspiral(values)
      if values[-1]:
        found.append(sim)
      else:
        missed.append(sim)

    # done
    connection.close()
    dbtables.discard_connection_filename(zerofname, working_filename, verbose = verbose)
    dbtables.DBTable_set_connection(None)

  print >>sys.stderr, "\nFound = %d Missed = %d" % (len(found), len(missed))
  return found, missed


def trim_mass_space(eff, twodbin, minthresh=0.0, minM=25.0, maxM=100.0):
  """
  restricts array to only have data within the mass space and sets everything
  outside the mass space to some canonical value, minthresh
  """
  x = eff.shape[0]
  y = eff.shape[1]
  c1 = twodbin.centres()[0]
  c2 = twodbin.centres()[1]
  numbins = 0
  for i in range(x):
    for j in range(y):
      if c1[i] > c2[j] or (c1[i] + c2[j]) > maxM or (c1[i]+c2[j]) < minM: eff[i][j] = minthresh
      else: numbins+=1

def fix_masses(sims):
  """
  Function to duplicate the mass pairs to remove edge effects 
  on the equal mass line, takes a list of sim rows
  """
  sims2 = []
  for l in sims:
    l2 = copy.deepcopy(l)
    l2.mass1 = l.mass2
    l2.mass2 = l.mass1
    sims2.append(l2)
  sims.extend(sims2)

def get_2d_mass_bins(low, high, bins):
  """
  Given the component mass range low, high of the search it will
  return 2D bins with size bins in each direction
  """
  mass1Bin = rate.LinearBins(low,high,bins)
  mass2Bin = rate.LinearBins(low,high,bins)
  twoDMB=rate.NDBins( (mass1Bin,mass2Bin) )
  return twoDMB
    
def scramble_pop(m, f):
  """
  A function to draw a new injection sample in the "boot strap" method 
  http://en.wikipedia.org/wiki/Bootstrapping_(statistics) 
  and included refereneces.
  This was used in the stack-a-flare search to get MC errors etc. 
  """
  inj = m+f
  ix = scipy.random.randint(0,len(inj), (len(inj),))
  return [inj[i] for i in ix if i < len(m) ], [inj[i] for i in ix if i >=len(m) ]

def scramble_dist(dist,relerr):
  """
  function to handle random calibration error.  Individually srambles the distances
  of injection by an error.
  """
  #return dist * float( scipy.exp( relerr * scipy.random.standard_normal(1) ) )
  return dist * (1-relerr)

def twoD_SearchVolume(found, missed, twodbin, dbin, wnfunc, livetime, bootnum=1):
  """ 
  Compute the search volume in the mass/mass plane, bootstrap
  and measure the first and second moment (assumes the underlying 
  distribution can be characterized by those two parameters) 
  This is gonna be brutally slow
  """
  wnfunc /= wnfunc[(wnfunc.shape[0]-1) / 2, (wnfunc.shape[1]-1) / 2]
  x = twodbin.shape[0]
  y = twodbin.shape[1]
  z = dbin.n
  rArrays = []
  volArray=rate.BinnedArray(twodbin)
  volArray2=rate.BinnedArray(twodbin)
  #set up ratio arrays for each distance bin
  for k in range(z):
    rArrays.append(rate.BinnedRatios(twodbin))

  # Bootstrap to account for errors
  for n in range(bootnum):
    #initialize by setting these to zero
    for k in range(z):
      rArrays[k].numerator.array = numpy.zeros(rArrays[k].numerator.bins.shape)
      rArrays[k].denominator.array = numpy.zeros(rArrays[k].numerator.bins.shape)
    #Scramble the inj population
    if bootnum > 1: sm, sf = scramble_pop(missed, found)
    else: sm, sf = missed, found
    for l in sf:#found:
      tbin = rArrays[dbin[scramble_dist(l.distance,0.1)]]
      tbin.incnumerator( (l.mass1, l.mass2) )
    for l in sm:#missed:
      tbin = rArrays[dbin[scramble_dist(l.distance,0.1)]]
      tbin.incdenominator( (l.mass1, l.mass2) )
    
    tmpArray2=rate.BinnedArray(twodbin) #start with a zero array to compute the mean square
    for k in range(z): 
      tbins = rArrays[k]
      tbins.denominator.array += tbins.numerator.array
      rate.filter_array(tbins.denominator.array,wnfunc)
      rate.filter_array(tbins.numerator.array,wnfunc)
      tbins.regularize()
      # logarithmic(d)
      integrand = 4.0 * pi * tbins.ratio() * dbin.centres()[k]**3 * dbin.delta
      volArray.array += integrand
      tmpArray2.array += integrand #4.0 * pi * tbins.ratio() * dbin.centres()[k]**3 * dbin.delta
      print >>sys.stderr, "bootstrapping:\t%.1f%% and Calculating smoothed volume:\t%.1f%%\r" % ((100.0 * n / bootnum), (100.0 * k / z)),
    tmpArray2.array *= tmpArray2.array
    volArray2.array += tmpArray2.array
    
  print >>sys.stderr, "" 
  #Mean and variance
  volArray.array /= bootnum
  volArray2.array /= bootnum
  volArray2.array -= volArray.array**2 # Variance
  volArray.array *= livetime
  volArray2.array *= livetime*livetime # this gets two powers of live time
  return volArray, volArray2
 

def cut_distance(sims, mnd, mxd):
  """
  Exclude sims outside some distance range to avoid errors when binning
  """
  return [sim for sim in sims if mnd <= sim.distance <= mxd]
 

######################## ACTUAL PROGRAM #######################################
###############################################################################
###############################################################################


def parse_command_line():
  parser = OptionParser(version = "%prog CVS $Id$", usage = "%prog [options] [file ...]", description = "%prog computes mass/mass upperlimit")
  parser.add_option("-i", "--instruments", metavar = "name[,name,...]", help = "Set the list of instruments.  Required.  Example \"H1,H2,L1\"")
  parser.add_option("--live-time-program", default = "thinca", metavar = "name", help = "Set the name of the program whose rings will be extracted from the search_summary table.  Default = \"thinca\".")
  parser.add_option("-t", "--output-name-tag", default = "", metavar = "name", help = "Set the file output name tag, real name is 2Dsearchvolume-<tag>-<ifos>.xml")
  parser.add_option("-f", "--full-data-file", default = "FULL_DATACAT_3.sqlite", metavar = "pattern", help = "File in which to find the full data, example FULL_DATACAT_3.sqlite")
  parser.add_option("-s", "--inj-data-glob", default = "*INJCAT_3.sqlite", metavar = "pattern", help = "Glob for files to find the inj data, example *INJCAT_3.sqlite")
  parser.add_option("-b", "--bootstrap-iterations", default = 1, metavar = "integer", type = "int", help = "Number of iterations to compute mean and variance of volume MUST BE GREATER THAN 1 TO GET USABLE NUMBERS, a good number is 10000")
  parser.add_option("--veto-segments-name", help = "Set the name of the veto segments to use from the XML document.")
  parser.add_option("--verbose", action = "store_true", help = "Be verbose.")

  opts, filenames = parser.parse_args()

  if opts.instruments is None:
    raise ValueError, "missing required argument --instruments"
  opts.instruments = lsctables.instrument_set_from_ifos(opts.instruments)

  opts.injfnames = glob.glob(opts.inj_data_glob)

  return opts, filenames

opts, filenames = parse_command_line()

if opts.veto_segments_name is not None:
  working_filename = dbtables.get_connection_filename(opts.full_data_file, verbose = opts.verbose)
  connection = sqlite3.connect(working_filename)
  dbtables.DBTable_set_connection(connection)
  veto_segments = db_thinca_rings.get_veto_segments(connection, opts.veto_segments_name)
  connection.close()
  dbtables.discard_connection_filename(opts.full_data_file, working_filename, verbose = opts.verbose)
  dbtables.DBTable_set_connection(None)
else:
  veto_segments = segments.segmentlistdict()


# FIXME:  don't you need to request the min(FAR) for the given "on"
# instruments?
FAR, seglists = get_far_threshold_and_segments(opts.full_data_file, opts.live_time_program, verbose = opts.verbose)


# times when only exactly the required instruments are on
seglists -= veto_segments
zero_lag_segments = seglists.intersection(opts.instruments) - seglists.union(set(seglists.keys()) - opts.instruments)

print FAR, abs(zero_lag_segments)


Found, Missed = get_injections(opts.injfnames, FAR, zero_lag_segments, verbose = opts.verbose)


# restrict the sims to a distance range
Found = cut_distance(Found, 1, 2000)
Missed = cut_distance(Missed, 1, 2000)

# get a 2D mass binning
twoDMassBins = get_2d_mass_bins(1, 99, 50)

# get log distance bins
dBin = rate.LogarithmicBins(0.1,2500,100)

gw = rate.gaussian_window2d(3,3,4)

#Get derivative of volume with respect to FAR
dvA = get_volume_derivative(opts.injfnames, twoDMassBins, dBin, FAR, zero_lag_segments, gw)

print >>sys.stderr, "computing volume at FAR " + str(FAR)
vA, vA2 = twoD_SearchVolume(Found, Missed, twoDMassBins, dBin, gw, float(abs(zero_lag_segments)), opts.bootstrap_iterations)

#Trim the array to have sane values outside the total mass area of interest
trim_mass_space(dvA, twoDMassBins, minthresh=0.0, minM=25.0, maxM=100.0)
trim_mass_space(vA, twoDMassBins, vA.min()/10.0, minM=25.0, maxM=100.0)
trim_mass_space(vA2, twoDMassBins, minthresh=0.0, minM=25.0, maxM=100.0)

#output an XML file with the result
xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())
xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(vA, "2DsearchvolumeFirstMoment"))
xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(vA2, "2DsearchvolumeSecondMoment"))
xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(dvA, "2DsearchvolumeDerivative"))
utils.write_filename(xmldoc, "2Dsearchvolume-%s-%s.xml" % (opts.output_name_tag, "".join(sorted(opts.instruments))))
