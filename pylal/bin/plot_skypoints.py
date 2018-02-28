#!/usr/bin/python

from glue import git_version

__author__ = "Larry Price <larry.price@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

import sys
from optparse import *
import glob
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import numpy as np
from glue.ligolw import utils, table
from pylal import skylocutils

usage = """
usage: %prog [options] 

Make a skymap from the output of run_skypoints.py.
Requires that the xml and corresponding .txt files are in the same directory.
NOTE: Requires the python basemap module!

"""

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage, version="%prog "+__version__ )

  # options related to input and output
  parser.add_option("-g","--glob",action="store",type="string",\
      default=None, metavar=" GLOB",help="GLOB of skypoints xml files" )

  (options,args) = parser.parse_args()

  return options, sys.argv[1:]
  
opts, args = parse_command_line()

#deal with the glob 
files = []
if opts.glob is not None:
  for gl in opts.glob.split(" "):
    files.extend(glob.glob(gl))
  if len(files) < 1:
    print >>sys.stderr, "The glob for " + opts.glob + " returned no files" 
    sys.exit(1)
else:
  print >>sys.stderr, "Need to specify a glob"
  sys.exit(1)

for file in files:
  xmldoc =utils.load_filename(file)
  try:
    sltab = table.get_table(xmldoc,skylocutils.SkyLocInjTable.tableName)
    if sltab[0].grid:
      print "injection"    
      inj = True
    else:
      print "no injection"
      sltab = table.get_table(xmldoc,skylocutils.SkyLocTable.tableName)
      inj = False  
  except:
    print "no injection"
    sltab = table.get_table(xmldoc,skylocutils.SkyLocTable.tableName)
    inj = False
  
  print "Plotting "+file
  for row in sltab:
    for fname in [row.grid,row.galaxy_grid]:
      if fname:
        if inj:
          injlat = row.dec*180/np.pi
          injlon = row.ra*180/np.pi
        snr = row.comb_snr
        latitude = []
        longitude = []
        color = []
        f = open(fname,'r')
        #set limits for the colors
        #these are chosen by trial and error and likely needs more tweaking
        cutoff = .00010
        minc = cutoff
        maxc = 0.003
        for line in f:
          if '#' in line:
            continue
          if '=' in line:
            continue
          elif line.strip():
            try:
              lon, lat, c = line.strip().split()
            except ValueError:
              continue
            #remove everything with a probability below the cutoff
            #this makes the skymaps cleaner
            if float(c) >= cutoff:
              latitude.append(float(lat)*180/np.pi)
              longitude.append(float(lon)*180/np.pi)
              color.append(float(c))

        f.close()

        #make sure we've got something to plot!
        if not latitude or not longitude or not color:
          print >>sys.stderr, "Unable to plot " + fname
          continue

        #sort them from least to most likely to avoid 
        #strange plotting artifacts
        sortme = zip(color,latitude,longitude)
        sortme.sort()
        bestlat = sortme[-1][1]
        bestlon = sortme[-1][2]
        color, latitude, longitude = zip(*sortme)

        pyplot.clf()
        m = Basemap(projection='moll', lon_0=180.0, lat_0 = 0.0)
        plx,ply = m(np.asarray(longitude),np.asarray(latitude))
        pyplot.scatter(plx,ply,s=5,c=np.asarray(color),faceted=False,cmap=cm.jet,vmin=minc,vmax=maxc)
        m.drawmapboundary()
        m.drawparallels(np.arange(-90.,120.,45.),labels=[1,0,0,0],labelstyle='+/-') # draw parallels
        m.drawmeridians(np.arange(0.,420.,90.),labels=[0,0,0,1],labelstyle='+/-') # draw meridians
        pyplot.title("Skymap for "+str(fname) + "\n"+"SNR = "+str(snr)) # add a title
        pyplot.colorbar()
        if inj:
          injx,injy = m(np.asarray(injlon),np.asarray(injlat))
          pyplot.scatter([injx],[injy],s=60,c='#8B008B',marker=(6,1,0),faceted=False)
        bestx,besty = m(np.asarray(bestlon),np.asarray(bestlat))
        pyplot.scatter([bestx],[besty],s=60,c='#FF1493',marker=(6,1,0),faceted=False)
        pyplot.savefig(fname.replace('txt','png'))
