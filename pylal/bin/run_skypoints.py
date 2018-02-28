#!/usr/bin/python

from glue import git_version

__author__ = "Larry Price <larry.price@ligo.org> and Patrick Brady <patrick.brady@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

import os
import sys
import glob
import cPickle
from optparse import *
from math import sqrt
from numpy import zeros, ceil
from pylal import skylocutils
from glue.ligolw import ligolw, lsctables

##############################################################################
#
#          options
#
##############################################################################

usage = """
usage: %prog [options] 

Estimate the sky position from a coincident trigger.

"""


def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage, version="%prog "+__version__ )

  # options related to input and output
  parser.add_option("-g","--glob",action="store",type="string",\
      default=None, metavar=" GLOB",help="GLOB of files to read" )
  parser.add_option("-G","--grids",action="store",type="string",\
      default=None, metavar=" GRID",help="pickled sky grids (generated with make_skypoints_grids.py)")
  parser.add_option("-R","--ranks",action="store",type="string",\
      default=None, metavar=" RANKS",help="pickled ranking object (generated with make_skypoints_rankings.py)")
  parser.add_option("-u","--galaxy-priors-dir",action="store",type="string",\
      default=None, metavar=" PRIDIR", help="path to a directory containg pickles for using galaxy catalog priors (generated with make_skypoints_galaxy_priors.py)")
  parser.add_option("-o","--output-prefix",action="store",type="string",default='',\
                    help="appends ouput-prefix to output file names")
  parser.add_option("-z","--input-type",action="store",default="coinctable",\
                    help="specify the type of input in the glob.  valid options are coinctable (DEFAULT) and coire")
  parser.add_option("-c","--coarse-cut",action="store",type="float", default=0.01, metavar=" COARSE_CUT", \
    help="probability (based on timing only) required for moving from the coarse to the fine grid. default is 0.01 = 1% probability." )

  (options,args) = parser.parse_args()

  return options, sys.argv[1:]



##############################################################################
#
#          i/o setup
#
##############################################################################

opts, args = parse_command_line()

#dump the args into a single string to write with the grids
argstring = " ".join(args)


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
#put the files into the coinc data structure
coincs = skylocutils.Coincidences(files,opts.input_type)

#setup the output xml tables
xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())
skyloctable = lsctables.New(skylocutils.SkyLocTable)
skylocinjtable = lsctables.New(skylocutils.SkyLocInjTable)
xmldoc.childNodes[0].appendChild(skyloctable)
xmldoc.childNodes[0].appendChild(skylocinjtable)
#FIXME: insert a process params table!


#make it work with the XSL stylesheet
ligolw.Header += u"""\n\n"""\
                 +u"""<?xml-stylesheet type="text/xsl" href="ligolw.xsl"?>"""\
                 +u"""\n\n"""

#setup the output filenames
base_name = 'SKYPOINTS' + opts.output_prefix
post_fname = base_name + '_posterior_GPSTIME.txt'
prob_fname = base_name + '_probability_GPSTIME.txt'
gal_fname = base_name + '_posterior_galaxy_prior_GPSTIME.txt'
outfile = base_name + '_GPSTIME.xml'

##############################################################################
#
#          convenience functions
#
##############################################################################

def get_unique_filename(name):
  """
  use this to avoid name collisions
  """
  counter = 1
  base_name, ext = os.path.splitext(name)
  while os.path.isfile(name):
    name = base_name + '_' + str(counter) + ext
    counter += 1

  return name

##############################################################################
#
#          main program
#
##############################################################################

#open up the pickled grids
gridfile = open(opts.grids,'r')
griddata = cPickle.load(gridfile)
grid = griddata['grids']
#fbins = griddata['skybins']
coarse_res = griddata['coarse_res']
fine_res = griddata['fine_res']
gridfile.close()

#open up the pickled pdfs
rankfile = open(opts.ranks,'r')
rankings = cPickle.load(rankfile)
rankfile.close()
dtr = rankings['dt']
dDr = rankings['dD']
Pdt = rankings['Pdt']
ref_freq = rankings['ref_freq']
snr_threshold = rankings['snr_threshold']

#add reference frequency and snr threshold arguments to metadata
argstring += ' --reference-frequency='+str(ref_freq)+\
             ' --snr-threshold='+str(snr_threshold) + '\n'

#the area of each pixel on the fine grid in square degrees
#this gets recorded for each point and makes computing areas simple
fine_area = fine_res*fine_res
coarse_area = coarse_res*coarse_res

for coinc in coincs:
  if len(coinc.ifo_list) < 3:
    continue
  sp = skylocutils.SkyPoints()
  
  #compute combined snr if snr dependent thresholds are specified
  rhosquared = 0.0
  for ifo in coinc.ifo_list:
    rhosquared += coinc.snr[ifo]*coinc.snr[ifo]
  combsnr = sqrt(rhosquared)
  if snr_threshold:
    dtsnrfac = combsnr/10.0
  else:
  #otherwise just multiply by unity
    dtsnrfac = 1.0
  
  #open up the necessary pickle with info on the galaxy prior
  if opts.galaxy_priors_dir:
    mineffD = ceil(min(coinc.eff_distances.values()))
    if mineffD > 50.: 
      mineffD = 50
    f = open(opts.galaxy_priors_dir+'/galaxy_prior_'+str(int(mineffD))+'Mpc.pkl','r')
    gal_prior = cPickle.load(f)
    f.close()

  print >>sys.stdout, 'Processing trigger at '+str(coinc.time)
  #main loop over the coarse grid
  for coarse_pt in grid.keys():

    #use timing alone to determine if we should move to the fine grid
    dtrss_coarse = dtsnrfac*skylocutils.get_delta_t_rss(coarse_pt,coinc,ref_freq)
    Pdt_coarse = 1 - Pdt.get_rank(dtrss_coarse)
    if Pdt_coarse >= opts.coarse_cut:
      #loop over points on the fine grid 
      for fine_pt in grid[coarse_pt]:
        dtrss_fine = dtsnrfac*skylocutils.get_delta_t_rss(fine_pt,coinc,ref_freq)
        dtrank = dtr.get_rank(dtrss_fine)
        dDrss_fine = skylocutils.get_delta_D_rss(fine_pt,coinc)
        dDrank = dDr.get_rank(dDrss_fine)
        L = dtrank*dDrank
        pval = 0.0
        if opts.galaxy_priors_dir:
          try:
            pval = gal_prior[fine_pt]
          except KeyError:
            pass
        sp.append([fine_pt,L,Pdt.get_rank(dtrank),dtrank,pval*L])

  fnames = {}
  if opts.input_type == 'coinctable':
    fnames['posterior'] = 'skymap_no_galaxies.txt'
    if opts.galaxy_priors_dir:
      fnames['galaxy'] = 'skymap.txt'
    else:
      fnames['galaxy'] = None
  else:
    fnames['posterior'] = get_unique_filename(post_fname.replace('GPSTIME',str(coinc.time.seconds)))
    if opts.galaxy_priors_dir:
      fnames['galaxy'] = get_unique_filename(gal_fname.replace('GPSTIME',str(coinc.time.seconds)))
    else:
      fnames['galaxy'] = None

  if sp:
    print >>sys.stdout, 'Populating sky localization table...'
    #populate the output tables
    #list of points has been sorted so the best one is at the top
    #FIXME: replace None with a link to the skymap file name!!!
    skylocutils.populate_SkyLocTable(skyloctable,coinc,sp,fine_area,fnames['posterior'],None,fnames['galaxy'])
  else:
    print >>sys.stdout, 'Unable to localize.'
    sys.exit(1)
  if coinc.is_injection:
    print >>sys.stdout, 'Recording injection data...'
    #NB: using the *recovered* snr for the snr dependent threshold
    inj_pt = (coinc.latitude_inj,coinc.longitude_inj)
    dtrss_inj = dtsnrfac*skylocutils.get_delta_t_rss(inj_pt,coinc,ref_freq)
    dtrank_inj = dtr.get_rank(dtrss_inj)
    dDrss_inj = skylocutils.get_delta_D_rss(inj_pt,coinc)
    dDrank_inj = dDr.get_rank(dDrss_inj)
    rank_inj = dtrank_inj*dDrank_inj
    area = {}
    if  opts.galaxy_priors_dir:
      try:
        galpt = fbins[skylocutils.sbin(fbins,inj_pt,fine_res)]
        pval_inj = gal_prior[galpt]
      except KeyError:
        pval_inj = 0.0
      gal_area = fine_area*len([pt for pt in sp if pt[4] >= rank_inj*pval_inj])
      area['gal'] = gal_area
    else:
      area['gal'] = None
    dt_area = fine_area*len([pt for pt in sp if pt[3] >= dtrank_inj])
    rank_area = fine_area*len([pt for pt in sp if pt[1] >= rank_inj])
    area['dt'] = dt_area
    area['rank'] = rank_area
    skylocutils.populate_SkyLocInjTable(skylocinjtable,coinc,rank_inj,area,\
                                        dtrss_inj,dDrss_inj,fnames['posterior'],fnames['galaxy'])

  #check for name collisions and then write the grid
  #use seconds of the smallest gpstime to label the event
  print >>sys.stdout, 'Writing skymap...'
  post_dat = {}
  post_dat['normfac'] = sp.normalize(1)
  post_dat['snr'] = combsnr
  post_dat['FAR'] = coinc.FAR
  post_dat['gps'] = coinc.time
  if opts.galaxy_priors_dir:
    post_dat['gnormfac'] = sp.normalize(4)
  sp.write(fnames,post_dat,argstring)
  
  print >>sys.stdout, 'Finished processing trigger.'

#name the xml file according to the range of gps times
if opts.input_type == 'coinctable':
  output = 'skypoints.xml'
else:
  if len(coincs) > 1:
    tmin = min([min(c.gps.values()) for c in coincs]).seconds
    tmax = max([max(c.gps.values()) for c in coincs]).seconds
    ofname=outfile.replace('GPSTIME',str(tmin)+'-'+str(tmax))
  #or the single time if only one coinc is present
  else:
    tmin = min([c for c in coincs[0].gps.values()])
    ofname=outfile.replace('GPSTIME',str(tmin))
  output = get_unique_filename(ofname)
#write the xml file and we're done
f = open(output,'w')
xmldoc.write(f)
f.close()

