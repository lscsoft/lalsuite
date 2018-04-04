#!/usr/bin/python
#
# Copyright (C) 2009  Chad Hanna, Kari Hodge, Anand Sengupta
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
import sys,os
from optparse import *
import glob
from time import clock,time
import re
from numpy import *
import math
import copy

from glue import iterutils
from glue import lal
from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from pylal import CoincInspiralUtils
from pylal import SnglInspiralUtils
from pylal import SearchSummaryUtils
from pylal import git_version
from pylal.tools import XLALCalculateEThincaParameter
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

def new_coincs_from_coincs(coincTable, coinc_stat):
  """
  We not only want to analyze each triple coincident trigger, but also the 3 double coincident
  triggers that you can make from the triple. Similarly, for the quadruples, we also want the 4 
  triples and 6 doubles that you can make from the quadruple. However, in order to be able to store
  these new (sub)coincidences, we have to assign each of them a unique event ID.
  """
  newCoincTable = CoincInspiralUtils.coincInspiralTable()
  id_generator = SnglInspiralUtils.SnglInspiralID_old()
  for row in coincTable:
    break_up_coinc(row, newCoincTable, coinc_stat, id_generator)
  return newCoincTable

def break_up_coinc(coincRow, coincTable, coinc_stat, id_generator):
  """
  This is the function that splits up the triples and quadruples into their respective subcoincidences
  """
  ifostr,ifolist=coincRow.get_ifos()
  new_ifos=[ifos for n in range(2, len(ifolist) + 1) for ifos in iterutils.choices(ifolist, n)]
  n_subcombs=len(new_ifos)
  for comb in new_ifos:
    n_ifos_in_comb=len(comb)
    t=getattr(coincRow,comb[0])
    event_id=t.event_id
    new_event_id=id_generator.new(t)
    myrow=CoincInspiralUtils.coincInspiralTable.row(new_event_id)
    for ifo in comb:
      myrow.add_trig(getattr(coincRow,ifo),coinc_stat)
      myt=getattr(myrow,comb[0])
    coincTable.append(myrow)

def ifos_to_list(ifos):
  ifofound = []
  ifolist = ['H1','H2','L1','V1']
  for i in range(len(ifos)/2):
    ifo = ifos[(2*i):(2*i+2)]
    ifofound.append(ifo)
  ifofound.sort()
  return ifofound

def coinc_to_info_param(coinc, param):
  if param.strip()[-1] == ')': return str(eval('coinc.'+param))
  else: return repr(int(getattr(coinc,param)))

def coinc_to_sngl_param_list(coinc, param):
  c_ifos,ifolist=coinc.get_ifos()
  ifolist.sort() # guarantee the order (you must keep ifo combinations in alphabetical order) 
  out = []
  for ifo in ifolist:
    t = getattr(coinc,ifo)
    if param.strip()[-1] == ')': out.append(str(eval('t.'+param)))
    else: out.append(repr(getattr(t,param)))
  return out

def coinc_to_metric_param(coinc,param):
  c_ifos,ifolist=coinc.get_ifos()
  ifolist.sort() # guarantee the order (you must keep ifo combinations in alphabetical order) 
  out=[]
  n=len(ifolist)
  for i in range(n):
    for j in range(i+1,n):
      t1=getattr(coinc,ifolist[i])
      t2=getattr(coinc,ifolist[j])
      out.append(repr(XLALCalculateEThincaParameter(t1,t2)))
  return out
           

def coinc_to_delta_param(coinc, param):
  par = coinc_to_sngl_param_list(coinc, param)
  return list_to_delta(par)

def coinc_to_rel_delta_param(coinc, param):
  par = coinc_to_sngl_param_list(coinc, param)
  return list_to_rel_delta(par)
 
def coinc_to_class_param(coinc):
  is_slide = int(bool(coinc._get_slide_num()))
  is_slide = 1 - is_slide
  is_slide = int(bool(is_slide))
  return str(is_slide) 
   
def list_to_delta(par):
  out = []
  for i in range(len(par)):
    for j in range(i+1,len(par)):
      out.append(repr(abs(float(par[i]) - float(par[j]))))
  return out

def list_to_rel_delta(par):
  out = []
  for i in range(len(par)):
    for j in range(i+1,len(par)):
      out.append(repr( (abs(float(par[i]) - float(par[j]))) / (abs(float(par[i]) + float(par[j])))/2.0 ))
  return out

# FIXME use the proper function
def mod_list(par):
  out = []
  for f in par: out.append(repr(float(f) % 1) )
  return out

def time_list_from_coinc(coinc):
  out = []
  for t in coinc.get_gps_times().values():
    out.append(float(t))
  return out
    
def parse_coinc(coincs, table, params):
  for coinc in coincs:
    tr = []
    ifos,ifolist=coinc.get_ifos()
    for p in params["single"]:
      tr.extend( coinc_to_sngl_param_list(coinc, p) )
    for p in params["metricInfo"]:
        tr.extend(coinc_to_metric_param(coinc,p))
    for p in params["coincDelta"]:
      # FIXME delta t corrected for time slides is a special case
      if p == 'time':
        tr.extend( mod_list(list_to_delta(time_list_from_coinc(coinc))) )
      else: tr.extend( coinc_to_delta_param(coinc, p) )
    for p in params["coincRelativeDelta"]:
      tr.extend( coinc_to_rel_delta_param(coinc, p) )
    for p in params["coincInfo"]:
      if p == 'class':
        tr.append(coinc_to_class_param(coinc))
      else:
       tr.append(coinc_to_info_param(coinc, p))
    table[ifos].append(tr)

def delta_ifo_str(par):
  out = []
  for i in range(len(par)):
    for j in range(i+1,len(par)):
      out.append(par[i] + par[j])
  return out

def make_table_header(ifos, params):
  # WARNING the order in which this is done must be kept in sync
  # with the parse coinc function
  ifolist = ifos_to_list(ifos)
  deltaifolist = delta_ifo_str(ifolist)
  tr = []
  for p in params["single"]:
    for ifo in ifolist:
      tr.append(ifo + p)
  for p in params["metricInfo"]:
    for ifo in deltaifolist:
      tr.append(ifo + p)
  for p in params["coincDelta"]:
    for ifo in deltaifolist: 
      tr.append(ifo + p)
  for p in params["coincRelativeDelta"]:
    for ifo in deltaifolist:
      tr.append(ifo + p)
  for p in params["coincInfo"]:
    if p=='class': continue
    tr.append(p)
  return tr

def table_dict_to_file(table,params,basename,types):
  totalrows = 0
  for t in types.keys():
    totalrows += len(types[t])
  rownum = 0
  cnt = 0
 
  for k in table.keys():
    f = {}
    for t in types.keys():
      f[t] = open(k+t+basename,'w')
      header = make_table_header(k, params)
      f[t].write(str(len(header))+"\n")
      if header: f[t].write(" ".join(header)+"\n")
    for r in table[k]:
      rownum = cnt % totalrows
      cnt += 1
      for t in types.keys():
        if rownum in types[t]:
          f[t].write(" ".join(r)+"\n")
    for t in types.keys():  
      f[t].close()

def merge_table_dict(t1,t2):
  for k in t2.keys():
    t1[k].extend(t2[k])
      
def get_coincs_from_cache(cachefile, pattern, match, verb, coinc_stat):
  cache = cachefile.sieve(description=pattern, exact_match=match)
  found, missed = cache.checkfilesexist()
  files = found.pfnlist()
  if not len(files):
    print >>sys.stderr, "cache contains no files with " + pattern + " description"
    return None
  # extract the coinc table
  coinc_table = SnglInspiralUtils.ReadSnglInspiralFromFiles(files, mangle_event_id=True, verbose=verb, non_lsc_tables_ok=False)
  # extract the list of coinc triggers
  return CoincInspiralUtils.coincInspiralTable(coinc_table,coinc_stat)

def get_slide_coincs_from_cache(cachefile, pattern, match, verb, coinc_stat):
  full_coinc_table = []
  cache = cachefile.sieve(description=pattern, exact_match=match)
  found, missed = cache.checkfilesexist()
  files = found.pfnlist()
  if not len(files):
    print >>sys.stderr, "cache contains no files with " + pattern + " description"
    return None
  # split the time slide files into 105 groups to aid with I/O
  num_files=len(files)

  #Changed by Tristan Miller as a memory fix
  #groups_of_files = split_seq(files,105)
  groups_of_files = split_seq(files,50)
  for filegroup in groups_of_files:
    if filegroup:  
      # extract the coinc table
      coinc_table = SnglInspiralUtils.ReadSnglInspiralFromFiles(filegroup, mangle_event_id=False, verbose=verb, non_lsc_tables_ok=False)
      segDict = SearchSummaryUtils.GetSegListFromSearchSummaries(filegroup)
      rings = segments.segmentlist(iterutils.flatten(segDict.values()))
      rings.sort()
      for k,ring in enumerate(rings):
        rings[k] = segments.segment(rings[k][0], rings[k][1] + 10**(-9))
      shift_vector = {"H1": 0, "H2": 0, "L1": 5, "V1": 5}
      if coinc_table:
        SnglInspiralUtils.slideTriggersOnRingWithVector(coinc_table, shift_vector, rings)
        full_coinc_table.extend(CoincInspiralUtils.coincInspiralTable(coinc_table,coinc_stat))
  return full_coinc_table

def split_seq(seq,p):
  newseq = []
  n = len(seq) / p # min items per subsequence
  r = len(seq) % p # remaindered items 
  b,e = 0l,n + min(1,r) # first split
  for i in range(p):
    newseq.append(seq[b:e])
    r = max(0, r-1)
    b,e = e,e + n + min(1,r) # min(1,r) is always 1 or 0
  return newseq
   
## MAIN PROGRAM ###
###################
usage="""%prog [options] 
This code takes the result of ihope and puts it into files readable by SprBaggerDecisionTreeApp, a C++ program in the software package StatPatternRecognition.

SprBaggerDecisionTreeApp creates a random forest of bagged decision trees. Each tree is a series of binary splits that classify an event (represented by a vector of parameters: SNR for each detector in the coincident ifo combination, chi^2 for each, r^2 veto duration for each, effectiveSNR^2 for each, pairwise differences in end time, pairwise relative differences in chirp mass) as either signal (0) or background (1). The forest trains on time slides (1) and injections (0). Zero lag triggers are then run through the trained forest, and the trees are averaged together, creating a ranking statsitic that is close to 0 for signal and 1 for background. This code creates the training and validation sets, which are ASCII tables (.pat extension) where the first row is the number of parameters, the second row is a list of the parameters, and subsequent rows contain the values of the parameters for a specific event, with the final entry being a 0 for injections and a 1 for time slides. Similarly for the testing set, whose events are all given the class 0 to start with. 
"""

__author__ = "Kari Hodge <khodge@ligo.caltech.edu>"
__prog__ = "pylal_ihope_to_randomforest_input"
__title__ = "Create input files for random forest from the ihope cache"


# parse the options and their arguments that you provide when you run the program
parser=OptionParser(usage=usage, version=git_version.verbose_msg)

parser.add_option("", "--cache-file", default="*ihope.cache",metavar="CACHEFILE", help="cache pointing to files of interest")
parser.add_option("","--playground-zerolag-pattern", default="COIRE_SECOND*_PLAYGROUND_CAT_2_VETO", action="store",type="string", metavar="PLAYZLTRIGPTTRN", help="sieve pattern for trig-files" )
parser.add_option("","--full-data-zerolag-pattern", default="COIRE_SECOND*CAT_2_VETO", action="store",type="string", metavar="FULLZLTRIGPTTRN", help="sieve pattern for trig-files" )
parser.add_option("","--slide-pattern", default="COIRE_SLIDE_SECOND_*_FULL_DATA_CAT_2_VETO", action="store",type="string", metavar="SLIDEPTTRN", help="sieve pattern for background" )
parser.add_option("","--found-pattern", default="COIRE_INJECTIONS_*_FOUND_SECOND_*_*INJ_CAT_2_VETO", metavar="FOUNDPTTRN", help="sieve pattern for found injection files")
parser.add_option("", "--ifo-times", default=None, action="store", type="string", metavar="IFOTIMES", help="sieve a cache file according to a particular ifo type")
parser.add_option("","--veto-category", default="CAT_2", action="store", type="string", metavar="VETOCAT", help="in the format CAT_x where x=1,2,3,4,5,6") 
parser.add_option("", "--exact-match",default=False, action="store_true", help="the pattern should match exactly if this option is used" )
parser.add_option("--statistic", action="store",type="string", default="effective_snr",help="choice of statistic used in making plots, valid arguments are: snr, snr_over_chi, s3_snr_chi_stat, effective_snr")
parser.add_option("--verbose",action="store_true",default=True,help="print extra info to the screen")
(opts,args)=parser.parse_args()

# get tables ready for putting things in
#instruments = ("H1", "H2", "L1", "V1")
#ZeroTable = dict(("".join(key), []) for n in range(2, len(instruments + 1)) for key in iterutils.choices(instruments, n))
PlaygroundZeroTable = {'H1H2':[],'H1L1':[],'H1V1':[],'H2L1':[],'H2V1':[],'L1V1':[],'H1H2L1':[],'H1L1V1':[],'H2L1V1':[],'H1H2V1':[],'H1H2L1V1':[]}
FullDataZeroTable = {'H1H2':[],'H1L1':[],'H1V1':[],'H2L1':[],'H2V1':[],'L1V1':[],'H1H2L1':[],'H1L1V1':[],'H2L1V1':[],'H1H2V1':[],'H1H2L1V1':[]}
SlideTable = {'H1H2':[],'H1L1':[],'H1V1':[],'H2L1':[],'H2V1':[],'L1V1':[],'H1H2L1':[],'H1L1V1':[],'H2L1V1':[],'H1H2V1':[],'H1H2L1V1':[]}
InjTable = {'H1H2':[],'H1L1':[],'H1V1':[],'H2L1':[],'H2V1':[],'L1V1':[],'H1H2L1':[],'H1L1V1':[],'H2L1V1':[],'H1H2V1':[],'H1H2L1V1':[]}
KnownTable = {'H1H2':[],'H1L1':[],'H1V1':[],'H2L1':[],'H2V1':[],'L1V1':[],'H1H2L1':[],'H1L1V1':[],'H2L1V1':[],'H1H2V1':[],'H1H2L1V1':[]}

# too keep I/O functioning smoothly, save the relevant lines from the cache into a new, smaller cache
cache_file = glob.glob(opts.cache_file)[0]
os.system('cat ' + cache_file + ' | grep COIRE | grep SECOND | grep ' + opts.veto_category + ' > mvsc_cache.cache') 
# now we can open the cache file
cachefile = lal.Cache.fromfile(open('mvsc_cache.cache')).sieve(ifos=opts.ifo_times, exact_match=opts.exact_match)
# initialize CoincInspiralUtils, which is going to pull coincidences out of the xml files that you provide as arguments to the options
coinc_stat=CoincInspiralUtils.coincStatistic(opts.statistic)
# Read in the files
SlideCoincs = get_slide_coincs_from_cache(cachefile, opts.slide_pattern, opts.exact_match, opts.verbose, coinc_stat)
if SlideCoincs: SlideCoincs = new_coincs_from_coincs(SlideCoincs, coinc_stat)
PlaygroundZeroCoincs = get_coincs_from_cache(cachefile, opts.playground_zerolag_pattern, opts.exact_match, opts.verbose, coinc_stat)
if PlaygroundZeroCoincs: PlaygroundZeroCoincs = new_coincs_from_coincs(PlaygroundZeroCoincs, coinc_stat)
FullDataZeroCoincs = get_coincs_from_cache(cachefile, opts.full_data_zerolag_pattern, opts.exact_match, opts.verbose, coinc_stat)
if FullDataZeroCoincs: FullDataZeroCoincs = new_coincs_from_coincs(FullDataZeroCoincs, coinc_stat)
InjCoincs = get_coincs_from_cache(cachefile, opts.found_pattern, opts.exact_match, opts.verbose, coinc_stat)
if InjCoincs: InjCoincs = new_coincs_from_coincs(InjCoincs,coinc_stat)

params = {"single":['snr','chisq','rsqveto_duration','get_effective_snr()','eff_distance','get_end()'], "metricInfo":['ethinca'],"coincRelativeDelta":['mchirp','eta'], "coincDelta":['time'], "coincInfo":['class']}

if SlideCoincs: parse_coinc(SlideCoincs,SlideTable,params)
if PlaygroundZeroCoincs: parse_coinc(PlaygroundZeroCoincs,PlaygroundZeroTable,params)
if FullDataZeroCoincs: parse_coinc(FullDataZeroCoincs,FullDataZeroTable,params)
if InjCoincs: parse_coinc(InjCoincs,InjTable,params)

# need to put the triggers we know to represent signal (injections) and background (timeslides) into the same file for testing/validation
merge_table_dict(KnownTable,InjTable)
merge_table_dict(KnownTable,SlideTable)

table_dict_to_file(KnownTable,params,'Known.pat',{'setA':[0],'setB':[1],'setC':[2],'setD':[3]})
table_dict_to_file(PlaygroundZeroTable,params,'ZeroLag_playground.pat',{'set':[0]})
table_dict_to_file(FullDataZeroTable,params,'ZeroLag_fulldata.pat',{'set':[0]})
