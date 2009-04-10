from pylal import rate
from pylal import SimInspiralUtils
import scipy
import numpy
import pylab
from math import *
import sys
import glob
import copy

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
  return dist * float( scipy.exp( relerr * scipy.random.standard_normal(1) ) )


def twoD_SearchVolume(found, missed, twodbin, dbin, wnfunc, bootnum=1):
  """ 
  Compute the search volume in the mass/mass plane, bootstrap
  to marginalize the errors 
  """
  x = twodbin.shape[0]
  y = twodbin.shape[1]
  z = dbin.n
  rArrays = []
  volArray=rate.BinnedArray(twodbin)
  #set up ratio arrays
  for k in range(z):
    rArrays.append(rate.BinnedRatios(twodbin))
  # Bootstrap to account for errors
  for n in range(bootnum):
    sm, sf = scramble_pop(missed, found)
    for l in sf:#found:
      tbin = rArrays[dbin[scramble_dist(l.distance,0.1)]]
      tbin.incnumerator( (l.mass1, l.mass2) )
    for l in sm:#missed:
      tbin = rArrays[dbin[scramble_dist(l.distance,0.1)]]
      tbin.incdenominator( (l.mass1, l.mass2) )
    print >>sys.stderr, "bootstrapping:\t%.1f%%\r" % (100.0 * n / bootnum),
  # Take the mean, marginalize
  tbin.denominator.array /= bootnum
  tbin.numerator.array /= bootnum 
  # make denom total, regularize and smooth
  print >>sys.stderr, "\n"
  for k in range(z): 
    tbins = rArrays[k]
    tbins.denominator.array += tbins.numerator.array
    rate.filter_array(tbins.denominator.array,wnfunc)
    rate.filter_array(tbins.numerator.array,wnfunc)
    tbins.regularize()
    # logarithmic(d)
    volArray.array += 4.0 * pi * tbins.ratio() * dbin.centres()[k]**3 * dbin.delta
    print >>sys.stderr, "Calculating smoothed volume:\t%.1f%%\r" % (100.0 * k / z),
  return volArray
 

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
Found = SimInspiralUtils.ReadSimInspiralFromFiles(glob.glob('*FOUND*.xml'),verbose=True)
Missed = SimInspiralUtils.ReadSimInspiralFromFiles(glob.glob('*MISSED*.xml'),verbose=True)


# replace these with pylal versions ?
fix_masses(Found)
fix_masses(Missed)

# restrict the sims to a distance range
cut_distance(Found, 1, 2000)
cut_distance(Missed, 1, 2000)

# get a 2D mass binning
twoDMassBins = get_2d_mass_bins(1, 99, 100)
#dBin = rate.LinearBins(0,2000,200)
dBin = rate.LogarithmicBins(0.1,2500,200)


gw = rate.gaussian_window2d(11,11)
#FIXME make search volume above loudest event
vA = twoD_SearchVolume(Found, Missed, twoDMassBins, dBin, gw, 1)

#output an XML file with the result
xml = rate.binned_array_to_xml(vA, "2Dsearchvolume")
xml.write(open("2Dsearchvolume.xml","w"))
