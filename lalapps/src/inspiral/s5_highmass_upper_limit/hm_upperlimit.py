from pylal import rate
from pylal import SimInspiralUtils
import scipy
import numpy
#import pylab
from math import *
import sys
import glob
import copy
from glue.ligolw import ligolw
#try:
#        import sqlite3
#except ImportError:
#        # pre 2.5.x
#        from pysqlite2 import dbapi2 as sqlite3
#from glue.ligolw import dbtables


class sim:
  def __init__(self, distance, mass1, mass2):
    self.distance = distance
    self.mass1 = mass1
    self.mass2 = mass2

def get_injections(injfnames,ifos='H1H2L1'):
  """
  """
  found = []
  missed = []
  for f in injfnames:
    working_filename = dbtables.get_connection_filename(f, tmp_path = None, verbose = True)
    connection = sqlite3.connect(working_filename)
    cursor = connection.cursor()
    #Create a view of the sim table where the intersection of instrument time restricts the injections (i.e. add a WHERE clause to look at only injections that were made during a given ifo time
    try: cursor.execute('DROP VIEW sims;') 
    except: pass
    try: cursor.execute('DROP VIEW found;')
    except: pass
    sim_query = 'CREATE VIEW sims AS SELECT * FROM sim_inspiral WHERE EXISTS (SELECT in_start_time, in_end_time FROM search_summary WHERE (ifos=="'+ifos+'" AND in_start_time < sim_inspiral.geocent_end_time AND in_end_time > sim_inspiral.geocent_end_time));'
    cursor.execute(sim_query)
    #cursor.execute('CREATE VIEW sims AS SELECT * FROM sim_inspiral;') 

    #First find injections that were not found even above threshold (these don't participate in coincs)
    query = 'SELECT distance, mass1, mass2 FROM sims WHERE NOT EXISTS(SELECT * FROM coinc_event_map WHERE (sims.simulation_id == coinc_event_map.event_id AND coinc_event_map.table_name="sim_inspiral"));'
    for (distance, mass1, mass2) in cursor.execute(query):
        missed.append( sim(distance,mass1,mass2) )
    
    #Now we have to find the found / missed above the loudest event.  
    query = '''CREATE VIEW found AS
                 SELECT *
                 FROM coinc_event 
                 JOIN coinc_event_map AS mapa ON (coinc_event.coinc_event_id = mapa.coinc_event_id) 
                 JOIN sims ON (sims.simulation_id == mapa.event_id)
                 WHERE (
                       coinc_event.coinc_def_id == "coinc_definer:coinc_def_id:2" 
                       AND mapa.table_name == "sim_inspiral"
                       );'''

    cursor.execute(query)

    #FIXME Require the coincs to have a FA rate above some value
    exists_query = '''SELECT * FROM found JOIN coinc_event_map AS mapb ON (found.coinc_event_id = mapb.coinc_event_id)
                               JOIN coinc_inspiral ON (mapb.event_id == coinc_inspiral.coinc_event_id)
                               WHERE (mapb.table_name == "coinc_event")'''

    # The found injections meet the exists query
    query = 'SELECT distance, mass1, mass2 FROM found WHERE EXISTS(' + exists_query + ');'
    print query
    for (distance, mass1, mass2) in cursor.execute(query):
        found.append( sim(distance,mass1,mass2) )
    # The missed injections do not
    query = 'SELECT distance, mass1, mass2 FROM found WHERE NOT EXISTS(' + exists_query + ');'
    for (distance, mass1, mass2) in cursor.execute(query):
        missed.append( sim(distance,mass1,mass2) )

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

def twoD_SearchVolume(found, missed, twodbin, dbin, wnfunc, bootnum=1):
  """ 
  Compute the search volume in the mass/mass plane, bootstrap
  and measure the first and second moment (assumes the underlying 
  distribution can be characterized by those two parameters) 
  This is gonna be brutally slow
  """
  #print wnfunc.data
  #print wnfunc.shape
  #index = (wnfunc.shape[0]-1) / 2 * wnfunc.shape[1] + (wnfunc.shape[1]-1) / 2
  #print index
  #print wnfunc
  #print max(wnfunc.data)
  #print len(wnfunc.data)
  wnfunc /= wnfunc[(wnfunc.shape[0]-1) / 2, (wnfunc.shape[1]-1) / 2]
  x = twodbin.shape[0]
  y = twodbin.shape[1]
  z = dbin.n
  rArrays = []
  volArray=rate.BinnedArray(twodbin)
  volArray2=rate.BinnedArray(twodbin)
  #MCErrorArray = rate.BinnedArray(twodbin)
  #one = numpy.ones(twodbin.shape)
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
    sm, sf = scramble_pop(missed, found)
    for l in sf:#found:
      tbin = rArrays[dbin[scramble_dist(l.distance,0.1)]]
      tbin.incnumerator( (l.mass1, l.mass2) )
    for l in sm:#missed:
      tbin = rArrays[dbin[scramble_dist(l.distance,0.1)]]
      tbin.incdenominator( (l.mass1, l.mass2) )
    #print >>sys.stderr, "bootstrapping:\t%.1f%%\r" % (100.0 * n / bootnum),
    # make denom total, regularize and smooth
    
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
      #print tbins.denominator.array
      #MCErrorArray.array += 4.0 * pi * tbins.ratio() * ( 1.0-tbins.ratio() ) / tbins.denominator.array * dbin.centres()[k]**3 * dbin.delta
      print >>sys.stderr, "bootstrapping:\t%.1f%% and Calculating smoothed volume:\t%.1f%%\r" % ((100.0 * n / bootnum), (100.0 * k / z)),
    tmpArray2.array *= tmpArray2.array
    volArray2.array += tmpArray2.array
    
  print >>sys.stderr, "\n\nDone\n" 
  #Mean and variance
  volArray.array /= bootnum
  volArray2.array /= bootnum
  volArray2.array -= volArray.array**2 # Variance
  #MCErrorArray.array /= bootnum
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



#FIXME needs to get injections from DB and separate different types
# read the sim inspiral table

#Found, Missed = get_injections(glob.glob("*INJCAT_3.sqlite"))

#print len(Found), len(Missed)
#sys.exit(0)

Found = SimInspiralUtils.ReadSimInspiralFromFiles(glob.glob('*FOUND*.xml'),verbose=True)
Missed = SimInspiralUtils.ReadSimInspiralFromFiles(glob.glob('*MISSED*.xml'),verbose=True)


# replace these with pylal versions ?
fix_masses(Found)
fix_masses(Missed)

# restrict the sims to a distance range
cut_distance(Found, 1, 2000)
cut_distance(Missed, 1, 2000)

# get a 2D mass binning
twoDMassBins = get_2d_mass_bins(1, 99, 50)
#dBin = rate.LinearBins(0,2000,200)
dBin = rate.LogarithmicBins(0.1,2500,200)


gw = rate.gaussian_window2d(5,5,4)
#FIXME make search volume above loudest event
vA, vA2 = twoD_SearchVolume(Found, Missed, twoDMassBins, dBin, gw, 100)

#output an XML file with the result
xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())
xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(vA, "2DsearchvolumeFirstMoment"))
xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(vA2, "2DsearchvolumeSecondMoment"))
#xmldoc.childNodes[-1].appendChild(rate.binned_array_to_xml(eA, "2DsearchvolumeMCError"))
xmldoc.write(open("2Dsearchvolume.xml","w"))
