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
from glue.ligolw import ligolw

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

vA = None
for f in glob.glob('2Dsearchvolume*.xml'):
  xml = ligolw_add.ligolw_add(ligolw.Document(), [f], verbose = True, non_lsc_tables_ok=True)
  if vA: vA += rate.binned_array_from_xml(xml, "2Dsearchvolume")
  else: vA = rate.binned_array_from_xml(xml, "2Dsearchvolume")

twoDMassBins = rate.bins_from_xml(xml)

#bin edges Number of bins + 1 for pcolor
X = numpy.array( list(twoDMassBins.lower()[0]) + [twoDMassBins.upper()[0][-1]] )
Y = numpy.array( list(twoDMassBins.lower()[1]) + [twoDMassBins.upper()[1][-1]] )


ul = pylab.log10(2.3 / (vA.array + 1)) + 3 # just a temporary scaling for the sake of making a plot that won't get everyone all worked up until the real calculation is done

trim_mass_space(25, 100, ul, twoDMassBins)

pylab.gray()
pylab.pcolor(X,Y, ul)
#pylab.hold(1)
#cn = pylab.contour(twoDMassBins.centres()[0], twoDMassBins.centres()[1], ul)
#pylab.clabel(cn)
#pylab.hold(0)
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Log10(Rate)",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.grid()
pylab.show()

