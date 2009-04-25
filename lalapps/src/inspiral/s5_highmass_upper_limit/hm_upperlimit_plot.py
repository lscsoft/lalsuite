from pylal import rate
from pylal import SimInspiralUtils
import scipy
import numpy
import pylab
from math import *
import sys
import glob
import copy
from glue.ligolw.utils import ligolw_add
from glue.ligolw import ligolw, utils

def trim_mass_space(minM, maxM, eff, twodbin,minthresh=-2):
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
      #if (c1[i] + c2[j]) > maxM or (c1[i]+c2[j]) < minM: eff[i][j] = 0
      else: numbins+=1
  print "found " + str(numbins) + " bins"


def istrue(arg):
  return True

vA = None
vA2 = None
eA = None
for f in glob.glob('2Dsearchvolume*.xml'):
  xmldoc = utils.load_filename(f, verbose=True, gz = (f or "stdin").endswith(".gz"))
  xmldoc = xmldoc.childNodes[0]
  if vA: vA += rate.binned_array_from_xml(xmldoc.childNodes[0], "2DsearchvolumeFirstMoment")
  else: vA = rate.binned_array_from_xml(xmldoc.childNodes[0], "2DsearchvolumeFirstMoment")
  if vA2: vA2 += rate.binned_array_from_xml(xmldoc.childNodes[1], "2DsearchvolumeSecondMoment")
  else: vA2 = rate.binned_array_from_xml(xmldoc.childNodes[1], "2DsearchvolumeSecondMoment")
  if eA: eA += rate.binned_array_from_xml(xmldoc.childNodes[2], "2DsearchvolumeDerivative")
  else: eA = rate.binned_array_from_xml(xmldoc.childNodes[2], "2DsearchvolumeDerivative")


twoDMassBins = rate.bins_from_xml(xmldoc.childNodes[0])

#bin edges Number of bins + 1 for pcolor
X = numpy.array( list(twoDMassBins.lower()[0]) + [twoDMassBins.upper()[0][-1]] )
Y = numpy.array( list(twoDMassBins.lower()[1]) + [twoDMassBins.upper()[1][-1]] )

print len(X), len(Y), vA.array.shape

#ul = pylab.log10(2.3 / (vA.array + 1)) + 3 # just a temporary scaling for the sake of making a plot that won't get everyone all worked up until the real calculation is done

log_vol = pylab.log10(vA.array)

#log_vol = vA.array

vol_error = vA2.array**0.5 / (vA.array + 0.0001)
#vol_error = eA.array**0.5

#der = eA.array

trim_mass_space(25, 100, log_vol, twoDMassBins, 5)
trim_mass_space(25, 100, vol_error, twoDMassBins, 0)
#trim_mass_space(25, 100, der, twoDMassBins, 0)


pylab.figure(1)
pylab.gray()
pylab.pcolor(X,Y, log_vol)
#pylab.hold(1)
#cn = pylab.contour(twoDMassBins.centres()[0], twoDMassBins.centres()[1], ul)
#pylab.clabel(cn)
#pylab.hold(0)
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Log10[< Volume >] in Mpc^3",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.grid()
#pylab.show()

pylab.figure(2)
pylab.pcolor(X,Y, vol_error)
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Fractional Error on Volume [std/mean]",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.grid()

pylab.figure(3)
pylab.pcolor(X,Y, der)
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Volume derivative",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.grid()


pylab.show()


