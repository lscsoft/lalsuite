#!/usr/bin/python

import sys
import os
import shutil
import re
import exceptions
import glob
import string

from optparse import *
from types import *
from pylab import *
from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils


##############################################################################
# Function to read in the source file and return a sorted list of its
# components according to the luminosity
##############################################################################
def read_source_file( source_file ):
  """
  read in the galaxy list, from an inspsrcs.dat type file
  @param source_file: input file name
  """
  f = open( source_file , "r")
  lines = f.readlines()
  f.close()

  distance = []
  luminosity = []

  f = open( 'inspsrcs.tmp', 'w')
  for line in lines:
    f.write(line)
    if line[0] != '#':
      b  = string.split(line)
      distance.append( float(b[3])/1000.0 )
      luminosity.append( float(b[4]) )
  f.close()

  return asarray(distance),asarray(luminosity)

##############################################################################
# Define a Hist Function
##############################################################################
def histng(xdata,ydata,xedges):

  """
  histogram xdata with edges specified xedges and yedges.  
  Can rescale the entries by lum_weight
  @param xdata:  array of data for parameter x
  @param ydata:  array of data for parameter y
  @param xedges: bin boundaries for parameter x
  """
  ng_x = zeros(len(xedges),'d')
  xstep = xedges[1] - xedges[0]
  
  for i in range(len(xdata)):
    l = int((xdata[i] - xedges[0])/xstep)

    if (l>=0 and l<len(xedges)):
      ng_x[l] += ydata[i]
 
  return ng_x

##############################################################################
usage = """usage: %prog [options] file1 (file2 file3)

Generate a fake population extrapolated from the Galaxy catalog put
together by Chad and Ravi

"""
parser = OptionParser( usage )
parser.add_option("-d","--min-distance",action="store",type="float",\
    default=25.0, metavar=" MINDIST",help="the closest galaxy to generate" ) 
parser.add_option("-D","--max-distance",action="store",type="float",\
    default=100.0, metavar=" MAXDIST",help="the furthest galaxy to generate" ) 
parser.add_option("-e","--seed",action="store",type="int",\
    default=123, metavar=" SEED",help="seed for random number generator" ) 
parser.add_option("-n","--nbins",action="store",type="int",\
    default=5, metavar=" NBINS",help="number of bins to use for extrapolation" ) 
parser.add_option("-N","--nmax",action="store",type="int",\
    default=5, metavar=" NMAX",help="bin at largest distance to use in extrapolation" ) 
parser.add_option("-s","--show-plot",action="store_true",default=False,\
    help="display the figures on the terminal" )
parser.add_option("-f","--fig-name",action="store",type="string",\
    default=None, metavar=" FNAME",\
    help="generate png figures with name FNAME-fig.png" ) 
parser.add_option("-S","--source-file",action="store",type="string",\
    default=None, metavar=" SOURCEFILE",help="galaxy source file" ) 
parser.add_option("-o","--outfile",action="store",type="string",\
    default=None, metavar=" OUTFILE",\
    help="name of the file to write simulated galaxy list to" ) 
parser.add_option("-V","--verbose",action="store_true",default=False,\
    help="print some more verbose output" )

(opts,args) = parser.parse_args()

# generate an injection file
command = "lalapps_bbhinj --seed " + str(opts.seed) 
command = command + "  --min-distance " + str(1.0e3 * opts.min_distance) 
command = command + "  --max-distance " + str(1.0e3 * opts.max_distance) 
command = command + "  --d-distr 2"
os.system(command)

# read the injection file and create a list of galaxies from it
bbhinjFile = "HL-INJECTIONS_123-729273613-5094000.xml"
doc = utils.load_filename(bbhinjFile)
sims = None
try: 
  simInspiralTable = \
      table.get_table(doc, lsctables.SimInspiralTable.tableName)
  sims = simInspiralTable
except: simInspiralTable = None
sims.sort(lambda a, b: cmp(a.distance, b.distance))

# make distance, latitute, longitude arrays
distance=sims.getColumnByName('distance').asarray()
latitude=sims.getColumnByName('latitude').asarray()
longitude=sims.getColumnByName('longitude').asarray()
ngalaxies=float(len(longitude))

figure()
plot( longitude, sin(latitude), 'bx')
xlabel(r'Longitude')
ylabel(r'Sin(Latitude)')
axis([0, 2*pi, -1, 1])
savefig('long-v-lat.png')

# compute the extrapolation constant
if opts.source_file:
  galdistance,galluminosity = read_source_file( opts.source_file )

  ledges=arange(0.0,25.0,1.0)
  hedges=arange(1.0,26.0,1.0)
  dv = (4.0/3.0) * pi * ( (hedges)**3 - (ledges)**3 )
  lumsum=histng(galdistance,galluminosity,ledges)

  suml = 0.0
  for i in arange(opts.nbins):
    n=opts.nmax-i
    l=(lumsum[n]/(dv[n]))
    suml = suml+l
    if opts.verbose:
      print "%f (%f) L_10/Mpc^3" % (l,suml/(float(i+1)))

  luminosity_extrapolation_factor = suml/float(opts.nbins)

  lpergalaxy = (luminosity_extrapolation_factor * \
      (4.0*pi/3.0)*(100.0**3-25.0**3) / ngalaxies )
else:
  lpergalaxy = 20.0


if opts.outfile:
  fp = open('inspsrcs.tmp', 'a')

  rahrmin = 24.0*longitude/6.28
  rahr = floor(rahrmin)
  ramin = 60.0*(rahrmin-rahr)

  decdegmin = 180.0*latitude/3.14
  decdeg = floor(abs(decdegmin))
  decmin = 60.0*(abs(decdegmin)-decdeg)
  
  for i,e in enumerate(decdegmin):
    if ( e < 0.0 ):
      decdeg[i] = -decdeg[i]
    fp.write("PB%05d \t %+02d:%02.2f \t %+02d:%02.2f \t %0.1f \t %0.3f \
        \t 1.0 \t 0.3 \t 0.3 \n" % (i+1,rahr[i],ramin[i],decdeg[i],decmin[i], \
          distance[i]*1000.0,lpergalaxy))

  fp.close()
  shutil.move('inspsrcs.tmp', opts.outfile)

galdistance,galluminosity = read_source_file( opts.outfile )


figure()
loglog(galdistance,cumsum(galluminosity),'ro',galdistance,luminosity_extrapolation_factor*(4.0*pi/3.0)*(galdistance)**3,'b-')
xlabel(r'Distance (Mpc)')
ylabel(r'Blue Luminosity (L$_{10}$)')
axvline(25.0, linewidth=2, color='g',linestyle='--')
axis([1.0,100.0,1.0,1.0e6])
savefig('number-v-distance.png')
show()
