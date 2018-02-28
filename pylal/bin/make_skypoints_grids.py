#!/usr/bin/python

from pylal import git_version

__author__ = "Larry Price <larry.price@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


import sys
import cPickle
from optparse import *
from pylal import skylocutils

usage = """
usage: %prog [options]

Create a pickle of grids for use by run_skypoints.py
See ligovirgocvs/cbc/protected/projects/s6/sky_localization for an example.

"""

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser(usage=usage)
  
  parser.add_option("-c","--coarse-resolution",action="store",type="float",\
      default=4.0, metavar=" CRES",help="side of a square pixel (in degrees) on the coarse grid" )
  parser.add_option("-f","--fine-resolution",action="store",type="float",\
      default=0.4, metavar=" FRES",help="side of a square pixel (in degrees) on the fine grid" )
  (options,args) = parser.parse_args()

  return options, sys.argv[1:]

opts, args = parse_command_line()

#make the grids
coarse_grid = skylocutils.gridsky(opts.coarse_resolution,True)
fine_grid = skylocutils.gridsky(opts.fine_resolution,True)

#create a dictionary mapping the points in the grids
shifted_grids, leftovers = skylocutils.map_grids(coarse_grid,fine_grid,opts.coarse_resolution)

#make sure it worked
if leftovers:
  print >>sys.stderr, "Could not map all the points in the fine grid" 
  print >>sys.stderr, str(len(leftovers)) + " leftover points"
  sys.exit(1)

#FIXME: this should be moved upstream
#relabel the latitudes

def _shift_lats(pt):
  """
  shit latitudes in pt by -pi/2
  """
  return (pt[0]-skylocutils.pi/2,pt[1])

grids = {}
for cpt in shifted_grids.keys():
  grids[(cpt[0]-skylocutils.pi/2,cpt[1])] = [(fpt[0]-skylocutils.pi/2,fpt[1]) for fpt in shifted_grids[cpt]]

griddata = {}
griddata['grids'] = grids
griddata['coarse_res'] = opts.coarse_resolution
griddata['fine_res'] = opts.fine_resolution

f = open('skygrids.pkl','w')
cPickle.dump(griddata,f,protocol=2)
f.close


