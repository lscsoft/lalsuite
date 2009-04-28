from pylal import rate
from pylal import SimInspiralUtils
import scipy
from scipy import interpolate
import numpy
#import pylab
from math import *
import sys
import glob
import copy
from glue.ligolw import ligolw
from optparse import OptionParser

try:
        import sqlite3
except ImportError:
        # pre 2.5.x
        from pysqlite2 import dbapi2 as sqlite3
from glue.ligolw import dbtables


class sim:
  def __init__(self, distance, mass1, mass2):
    self.distance = distance
    self.mass1 = mass1
    self.mass2 = mass2

def get_zero_lag_far(zerofname, ifos):
  query = 'SELECT min(false_alarm_rate) FROM coinc_inspiral JOIN coinc_event ON (coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id) WHERE NOT EXISTS(SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset != 0);'
  working_filename = dbtables.get_connection_filename(zerofname, tmp_path = None, verbose = True)
  connection = sqlite3.connect(working_filename)
  cursor = connection.cursor()

  result = cursor.execute(query)
  query = 'SELECT sum(in_end_time - in_start_time) FROM search_summary WHERE (ifos == "' + ifos + '");'
  for r in result:
    far = r[0]

  result = cursor.execute(query)
  for r in result:
    time = r[0]
  cursor.close()
  return far, time

def get_volume_derivative(injfnames, twoDMassBins, dBin, FAR, live_time, ifos, gw):
  FARh = FAR*100000
  FARl = FAR*0.001
  nbins = 5
  FARS = rate.LogarithmicBins(FARl, FARh, nbins)
  #gw = rate.gaussian_window2d(3,3,10)
  vA = []
  vA2 = []
  for far in FARS.centres():
    m, f = get_injections(injfnames, far, ifos)
    print >>sys.stderr, "computing volume at FAR " + str(far)
    vAt, vA2t = twoD_SearchVolume(f, m, twoDMassBins, dBin, gw, live_time, 1)  
    vAt.array = scipy.log10(vAt.array + 0.001)
    vA.append(vAt)
  # the derivitive is calcuated with respect to FAR * t
  FARS = rate.LogarithmicBins(FARl * live_time, FARh * live_time, nbins)
  return derivitave_fit(FARS, FAR * live_time, vA, twoDMassBins)
  
  
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

def get_injections(injfnames, FAR, ifos):
  """
  """
  found = []
  missed = []
  cnt = 0
  print >>sys.stderr, ""
  for f in injfnames:
    print >>sys.stderr, "getting injections below FAR: " + str(FAR) + ":\t%.1f%%\r" % (100.0 * cnt / len(injfnames),),
    cnt= cnt+1
    working_filename = dbtables.get_connection_filename(f, tmp_path = None, verbose = True)
    connection = sqlite3.connect(working_filename)
    cursor = connection.cursor()
    #Create a view of the sim table where the intersection of instrument time restricts the injections (i.e. add a WHERE clause to look at only injections that were made during a given ifo time
    try: cursor.execute('DROP VIEW sims;') 
    except: pass
    try: cursor.execute('DROP VIEW found;')
    except: pass
    sim_query = 'CREATE TEMPORARY VIEW sims AS SELECT * FROM sim_inspiral WHERE EXISTS (SELECT in_start_time, in_end_time FROM search_summary WHERE (ifos=="'+ifos+'" AND in_start_time < sim_inspiral.geocent_end_time AND in_end_time > sim_inspiral.geocent_end_time));'
    cursor.execute(sim_query)
    #cursor.execute('CREATE VIEW sims AS SELECT * FROM sim_inspiral;') 

    #First find injections that were not found even above threshold (these don't participate in coincs)
    query = 'SELECT distance, mass1, mass2 FROM sims WHERE NOT EXISTS(SELECT * FROM coinc_event_map WHERE (sims.simulation_id == coinc_event_map.event_id AND coinc_event_map.table_name="sim_inspiral"));'
    for (distance, mass1, mass2) in cursor.execute(query):
        missed.append( sim(distance,mass1,mass2) )
    
    #Now we have to find the found / missed above the loudest event.  
    query = '''CREATE TEMPORARY VIEW found AS
                 SELECT *
                 FROM coinc_event 
                 JOIN coinc_event_map AS mapa ON (coinc_event.coinc_event_id = mapa.coinc_event_id) 
                 JOIN sims ON (sims.simulation_id == mapa.event_id)
                 WHERE (
                       coinc_event.coinc_def_id == "coinc_definer:coinc_def_id:2" 
                       AND mapa.table_name == "sim_inspiral"
                       );'''

    cursor.execute(query)

    exists_query = '''SELECT * FROM coinc_event_map AS mapa JOIN coinc_event_map AS mapb ON (mapa.coinc_event_id == mapb.coinc_event_id) JOIN coinc_inspiral ON (mapb.table_name == "coinc_event" AND mapb.event_id == coinc_inspiral.coinc_event_id) WHERE mapa.table_name == "sim_inspiral" AND mapa.event_id == found.simulation_id AND coinc_inspiral.false_alarm_rate < ''' + str(FAR)

    # The found injections meet the exists query
    query = 'SELECT distance, mass1, mass2 FROM found WHERE EXISTS(' + exists_query + ');'
    for (distance, mass1, mass2) in cursor.execute(query):
        found.append( sim(distance,mass1,mass2) )
    # The missed injections do not
    query = 'SELECT distance, mass1, mass2 FROM found WHERE NOT EXISTS(' + exists_query + ');'
    for (distance, mass1, mass2) in cursor.execute(query):
        missed.append( sim(distance,mass1,mass2) )
    cursor.close()
  print >>sys.stderr, "\nFound = " + str(len(found)) + " Missed = " + str(len(missed))
  return found, missed



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
  of injection by a random error (log normal)
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
  for k in range(len(sims)):
    if sims[k].distance >= mxd or sims[k].distance <= mnd: sims.pop(k)
 

######################## ACTUAL PROGRAM #######################################
###############################################################################
###############################################################################

parser = OptionParser(
                version = "%prog CVS $Id$",
                usage = "%prog [options] [file ...]",
                description = "%prog computes mass/mass upperlimit"
        )
parser.add_option("-i", "--ifos", default="H1H2L1", metavar = "ifo1ifo2...", help = "Set on instruments to compute upper limit for. Example H1H2L1")
parser.add_option("-t", "--output-name-tag", default = "", metavar = "name", help = "Set the file output name tag, real name is 2Dsearchvolume-<tag>-<ifos>.xml")
parser.add_option("-f", "--full-data-file", default = "FULL_DATACAT_3.sqlite", metavar = "pattern", help = "File in which to find the full data, example FULL_DATACAT_3.sqlite")
parser.add_option("-s", "--inj-data-glob", default = "*INJCAT_3.sqlite", metavar = "pattern", help = "Glob for files to find the inj data, example *INJCAT_3.sqlite")
parser.add_option("-b", "--bootstrap-iterations", default = 1, metavar = "number", help = "Number of iterations to compute mean and variance of volume MUST BE GREATER THAN 1 TO GET USABLE NUMBERS, a good number is 10000")


opts, filenames = parser.parse_args()


ifos = opts.ifos #"H1H2L1"
iters = int(opts.bootstrap_iterations)

injfnames = glob.glob(opts.inj_data_glob)
FAR, livetime = get_zero_lag_far(opts.full_data_file, ifos)
#livetime/=31556926.0 # convert seconds to years

print FAR, livetime

Found, Missed = get_injections(injfnames, FAR, ifos)

#Found = SimInspiralUtils.ReadSimInspiralFromFiles(glob.glob('*FOUND*.xml'),verbose=True)
#Missed = SimInspiralUtils.ReadSimInspiralFromFiles(glob.glob('*MISSED*.xml'),verbose=True)


# replace these with pylal versions ?
fix_masses(Found)
fix_masses(Missed)

# restrict the sims to a distance range
cut_distance(Found, 1, 2000)
cut_distance(Missed, 1, 2000)

# get a 2D mass binning
twoDMassBins = get_2d_mass_bins(1, 99, 50)
# get log distance bins
dBin = rate.LogarithmicBins(0.1,2500,200)


gw = rate.gaussian_window2d(3,3,4)

#Get derivative of volume with respect to FAR
dvA = get_volume_derivative(injfnames, twoDMassBins, dBin, FAR, livetime, ifos, gw)

print >>sys.stderr, "computing volume at FAR " + str(FAR)
vA, vA2 = twoD_SearchVolume(Found, Missed, twoDMassBins, dBin, gw, livetime, iters)

#output an XML file with the result
xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())
xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(vA, "2DsearchvolumeFirstMoment"))
xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(vA2, "2DsearchvolumeSecondMoment"))
xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(dvA, "2DsearchvolumeDerivative"))
xmldoc.write(open("2Dsearchvolume-"+opts.output_name_tag+"-"+opts.ifos+".xml","w"))
