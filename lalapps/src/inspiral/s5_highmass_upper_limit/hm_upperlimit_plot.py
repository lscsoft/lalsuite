#!/usr/bin/python
import sys
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

# FIXME, I apparently don't know how to read XML as cleverly as I thought, how do I find the table
# without specifying the childnode number?
def get_combined_array(tablename, childnode):
  # FIXME assumes that all the xml files have the same binned array tables
  # Figure out the shape of the arrays in the file, make an array with one more
  # dimension, the number of files from sys.argv[1:]
  xmldoc = utils.load_filename(sys.argv[1], verbose=True, gz = (f or "stdin").endswith(".gz"))
  xmldoc = xmldoc.childNodes[0]
  A  = rate.binned_array_from_xml(xmldoc.childNodes[0], tablename) 
  bins = rate.bins_from_xml(xmldoc.childNodes[0])
  out = numpy.zeros((len(sys.argv[1:]),)+A.shape,dtype="float")
  # Read the data
  for i, f in enumerate(sys.argv[1:]):
    xmldoc = utils.load_filename(f, verbose=True, gz = (f or "stdin").endswith(".gz"))
    xmldoc = xmldoc.childNodes[0]
    out[i] = rate.binned_array_from_xml(xmldoc.childNodes[childnode], tablename).array
  A.array = numpy.zeros(A.shape)
  return bins, out, 


def istrue(arg):
  return True

# bins is the same for each call and ulA is an empty binnedArray that has the right shape
# that can hold the upperlimit when we get around to computing it later, so it is okay
# that bins, and ulA are overwritten in each call. vA, vA2 and dvA are the important ones
bins, vA, ulA = get_combined_array("2DsearchvolumeFirstMoment", 0)
bins, vA2, ulA = get_combined_array("2DsearchvolumeSecondMoment", 1)
bins, dvA, ulA = get_combined_array("2DsearchvolumeDerivative", 2)

#bin edges Number of bins + 1 for pcolor
X = numpy.array( list(bins.lower()[0]) + [bins.upper()[0][-1]] )
Y = numpy.array( list(bins.lower()[1]) + [bins.upper()[1][-1]] )

print len(X), len(Y), vA[0].shape


log_vol = pylab.log10(vA[0])


vol_error = vA2[0]**0.5 / (vA[0] + 0.0001)

der = dvA[0] #pylab.log10(eA.array)

pylab.figure(1)
pylab.gray()
pylab.pcolor(X,Y, log_vol)
#pylab.hold(1)
#cn = pylab.contour(bins.centres()[0], bins.centres()[1], ul)
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
pylab.pcolor(X,Y, der )
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Volume derivative",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.grid()


pylab.show()


