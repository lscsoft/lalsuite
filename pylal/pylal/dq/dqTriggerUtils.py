#!/usr/bin/env python

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import sys,os,re,math,datetime,glob,copy
from socket import getfqdn

from glue.ligolw import ligolw,table,lsctables,utils
from glue.ligolw.utils import process as ligolw_process
from glue import segments

from glue.lal import Cache as LALCache
from glue.lal import CacheEntry as LALCacheEntry

from pylal import date
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

from glue import git_version

from scipy import special
import numpy

__author__  = "Duncan Macleod <duncan.macleod@astro.cf.ac.uk>"
__version__ = "git id %s" % git_version.id
__date__    = git_version.date

"""
This module provides a bank of useful functions for manipulating triggers and trigger files for data quality investigations.
"""

# global regular expressions
trigsep = re.compile('[\t\s,]+')
_cchar_regex = re.compile('[-#%<!()_\[\]-{}:;\'\"\ ]')

# =============================================================================
# Define ETG options
# =============================================================================

_trig_regex  = re.compile('(burst|inspiral|ring)', re.I)
_burst_regex = re.compile('(burst|omega|kleine|kw|cwb|hacr)', re.I)
_cbc_regex   = re.compile('(ihope|inspiral|cbc)', re.I)
_ring_regex  = re.compile('(ring)', re.I)

def SnglTriggerTable(etg, columns=None):

  """
    Handy function to return the correct type of ligolw table for the given ETG.
  """

  if _burst_regex.search(etg):
    t = lsctables.New(lsctables.SnglBurstTable, columns=columns)
  elif _cbc_regex.search(etg):
    t = lsctables.New(lsctables.SnglInspiralTable, columns=columns)
  elif _ring_regex.search(etg):
    t = lsctables.New(lsctables.SnglRingdownTable, columns=columns)
  else:
    raise AttributeError("etg=%s not recognised by SnglTriggerTable." % etg)

  # set columns
  if columns:
    if isinstance(columns[0], str): columns = map(str.lower, columns)
    if isinstance(columns[0], unicode): columns = map(unicode.lower, columns)
    for c in t.columnnames:
      if c.lower() not in columns:
        idx = t.columnnames.index(c)
        t.columnnames.pop(idx)
        t.columntypes.pop(idx)

  return t

def SnglTrigger(etg):
  """
    Handy function to return the correct type of ligolw table row object for
    the given ETG
  """

  if _burst_regex.search(etg):
    return lsctables.SnglBurst()
  elif _cbc_regex.search(etg):
    return lsctables.SnglInspiral()
  elif _ring_regex.search(etg):
    return lsctables.SnglRingdown()
  else:
    raise AttributeError("etg=%s not recognised by SnglTrigger." % etg)

# =============================================================================
# Define get_time choice
# =============================================================================

def def_get_time(tableName, ifo=None):

  """
    Define the get_time() function for given table
  """

  get_time = None

  if ifo:  ifo = ifo[0]

  # if given an injection table:
  if re.match('sim', tableName):

    if re.search('inspiral',tableName):
      get_time = lambda row: row.get_end(site=ifo)  
    else:
      get_time = lambda row: row.get_time_geocent()


  # if given a sngl trigger table
  elif re.match('(sngl|multi|coinc)', tableName):

    if re.search('inspiral',tableName):
      get_time = lambda row: row.get_end()
    elif re.search('ringdown',tableName):
      get_time = lambda row: row.get_start()
    else:
      get_time = lambda row: row.get_peak()

  return get_time

# =============================================================================
# Convert list from text file into Sngl{Burst,Inspiral} object
# =============================================================================

def trigger(data, etg, ifo=None, channel=None, loadcolumns=None):

  """
    Reads the list object data and returns a Sngl{Burst,Inspiral} object
    with the correct columns seeded, based on the given string object etg.

    Arguments:

      data : [ string | list ]
        string or list object containing the data to be parsed.

      etg : [ "ihope" | "kw" | "omega" | "omegadq" ]
        Defines how the data list is understood and parsed, based on the
        standard column layouts (comma-, space-, or tab-delimited) for each ETG.
        "ihope" is read as per the .csv files from ihope_daily. "omegadq"
        is read as per the .txt or .clst produced by the detchar scripts used to
        plot omega online events.

  """

  # set up trig object
  trig = SnglTrigger(etg)

  # set up load columns
  if loadcolumns==None:
    loadcolumns = trig.__slots__

  # if given string, split on space, tab or comma
  if isinstance(data,str):
    data = trigsep.split(data.rstrip())

  # =====
  # ihope
  # =====
  if re.match('ihope',etg):
    # comma separated values are:
    # end_time,end_time_ns,ifo,snr,mass1,mass2,mtotal,eta,event_duration,
    # template_duration,eff_distance,chisq,chisq_dof,bank_chisq,bank_chisq_dof
    # cont_chisq,cont_chisq_dof

    trig.end_time              = int(data[0])
    trig.end_time_ns           = int(data[1])
    trig.ifo                   = str(data[2])
    trig.snr                   = float(data[3])
    trig.mass1                 = float(data[4])
    trig.mass2                 = float(data[5])
    trig.mtotal                = float(data[6])
    trig.eta                   = float(data[7])
    trig.event_duration        = float(data[8])
    trig.template_duration     = float(data[9])
    trig.eff_distance          = float(data[10])
    trig.chisq                 = float(data[11])
    trig.chisq_dof             = float(data[12])
    trig.bank_chisq            = float(data[13])
    try:
      trig.bank_chisq_dof      = float(data[14])
    except:
      trig.bank_chisq_dof = None
    try:
      trig.cont_chisq          = float(data[15])
      trig.cont_chisq_dof      = float(data[16])
    except ValueError:
      trig.cont_chisq = None
      trig.cont_chisq_dof = None

  # ============
  # kleine-welle
  # ============
  elif etg=='kw':
    # space separated values are:
    # peak_time.peak_time_ns start_time.start_time_ns stop_time.stop_time_ns 
    # central_freq,energy,amplitude,n_pix,significance,N

    # peak time
    peak  = LIGOTimeGPS(data[2])
    trig.peak_time             = peak.seconds
    trig.peak_time_ns          = peak.nanoseconds
    # start time
    start = LIGOTimeGPS(data[0])
    ms_start = start
    trig.start_time            = start.seconds
    trig.start_time_ns         = start.nanoseconds
    trig.ms_start_time         = ms_start.seconds
    trig.ms_start_time_ns      = ms_start.nanoseconds
    # end time
    stop = LIGOTimeGPS(data[1])
    ms_stop = stop
    trig.stop_time             = stop.seconds
    trig.stop_time_ns          = stop.nanoseconds
    trig.ms_stop_time          = ms_stop.seconds
    trig.ms_stop_time_ns       = ms_stop.nanoseconds
    # duration
    trig.duration              = stop-start
    trig.ms_duration           = ms_stop-ms_start
    # others
    trig.central_freq          = float(data[3])
    trig.peak_frequency        = float(data[3])
    energy                = float(data[4])
    trig.amplitude             = float(data[5])
    n_pix                 = float(data[6])
    #significance          = float(data[7])
    #N                     = float(data[8])
    trig.snr                   = math.sqrt(trig.amplitude-n_pix)

  # =====
  # omega
  # =====
  if etg=='omega' or etg=='wpipe':
    # space separated values are:
    # peak_time.peak_time_ns peak_frequency duration bandwidth amplitude
    # cluster_size cluster_norm_energy cluster_number

    # peak time
    peak = LIGOTimeGPS(data[0])
    trig.peak_time           = peak.seconds
    trig.peak_time_ns        = peak.nanoseconds
    # central frequency
    trig.central_freq        = float(data[1])
    trig.peak_frequency      = trig.central_freq
    # duration
    trig.duration            = LIGOTimeGPS(float(data[2]))
    halfduration             = LIGOTimeGPS(float(data[2])/2)
    # start time
    start = peak-halfduration
    ms_start = start
    trig.start_time          = start.seconds
    trig.start_time_ns       = start.nanoseconds
    trig.ms_start_time       = ms_start.seconds
    trig.ms_start_time_ns    = ms_start.nanoseconds
    # end time
    stop = peak+halfduration
    ms_stop = stop
    trig.stop_time           = stop.seconds
    trig.stop_time_ns        = stop.nanoseconds
    trig.ms_stop_time        = ms_stop.seconds
    trig.ms_stop_time_ns     = ms_stop.nanoseconds

    trig.ms_duration         = ms_stop-ms_start
    # bandwidth and flow,fhigh
    trig.bandwidth           = float(data[3])
    trig.flow                = trig.peak_frequency - 0.5*trig.bandwidth
    trig.fhigh               = trig.peak_frequency + 0.5*trig.bandwidth
    # energy
    trig.amplitude           = float(data[4])

    # cluster parameters
    if len(data)>5:
      trig.param_one_name      = 'cluster_size'
      trig.param_one_value     = float(data[5])
      trig.param_two_name      = 'cluster_norm_energy'
      trig.param_two_value     = float(data[6])
      trig.param_three_name    = 'cluster_number'
      trig.param_three_value   = float(data[7])

    # SNR
    trig.snr                 = math.sqrt(2*trig.amplitude)

  # =====
  # omega
  # =====
  if etg=='omegaspectrum':
    # space separated values are:
    # peak_time.peak_time_ns peak_frequency duration bandwidth amplitude
    # cluster_size cluster_norm_energy cluster_number

    # peak time
    peak = LIGOTimeGPS(data[0])
    trig.peak_time           = peak.seconds
    trig.peak_time_ns        = peak.nanoseconds
    # central frequency
    trig.central_freq        = float(data[1])
    trig.peak_frequency      = trig.central_freq

    trig.amplitude           = float(data[2])

    # SNR
    trig.snr                 = math.sqrt(trig.amplitude)

    # use chisq to stor PSD variance, 
    if len(data)>3:
      trig.chisq               = math.sqrt(math.sqrt(float(data[3])))/trig.snr

  # ==============
  # omegadq
  # ==============
  # follow Cadonati's clustering output
  if etg=='omegadq':
    # start time
    start = LIGOTimeGPS(data[0])
    trig.start_time              = start.seconds
    trig.start_time_ns           = start.nanoseconds
    # end time
    stop  = LIGOTimeGPS(data[1])
    trig.stop_time               = stop.seconds
    trig.stop_time_ns            = stop.nanoseconds
    # peak time
    peak  = LIGOTimeGPS(data[2])
    trig.peak_time               = peak.seconds
    trig.peak_time_ns            = peak.nanoseconds
    # duration
    trig.duration                = stop-start
    # bandwidth, and flow,fhigh,central_freq
    trig.flow                    = float(data[3])
    trig.fhigh                   = float(data[4])
    trig.bandwidth               = trig.fhigh - trig.flow
    trig.central_freq            = trig.flow  + 0.5*trig.bandwidth

    # MS params
    ms_start = LIGOTimeGPS(data[6])
    trig.ms_start_time           = ms_start.seconds
    trig.ms_start_time_ns        = ms_start.nanoseconds
    ms_stop  = LIGOTimeGPS(data[7])
    trig.ms_stop_time            = ms_stop.seconds
    trig.ms_stop_time_ns         = ms_stop.nanoseconds
    trig.ms_duration             = ms_stop-ms_start
    trig.ms_flow                 = float(data[8])
    trig.ms_fhigh                = float(data[9])
    trig.ms_bandwidth            = trig.ms_fhigh - trig.ms_flow
    trig.peak_frequency          = trig.ms_flow  + 0.5*trig.ms_bandwidth

    # SNR
    trig.amplitude               = float(data[11])
    trig.snr                     = math.sqrt(2*float(data[11]))
    trig.ms_snr                  = math.sqrt(2*float(data[12]))

  # ====
  # HACR
  # ====

  # based on output of hacr web interface
  if etg=='hacr':
    # peak time
    trig.peak_time               = int(data[0])
    trig.peak_time_ns            = int(float(data[1])*math.pow(10,9))
    trig.param_one_name          = 'peak_time_offset'
    trig.param_one_value         = float(data[1])
    # duration
    trig.duration                = float(data[4])
    # start time
    start = trig.get_peak()-trig.duration
    trig.start_time              = start.seconds
    trig.start_time_ns           = start.nanoseconds
    # end time
    stop = trig.get_peak()+trig.duration
    trig.stop_time               = stop.seconds
    trig.stop_time_ns            = stop.nanoseconds
    # bandwidth, and flow,fhigh,central_freq
    trig.central_freq            = float(data[2])
    trig.peak_frequency          = trig.central_freq
    trig.bandwidth               = float(data[3])
    trig.peak_frequency_error    = trig.bandwidth/trig.central_freq
    trig.flow                    = trig.central_freq - 0.5*trig.bandwidth
    trig.fhigh                   = trig.central_freq + 0.5*trig.bandwidth
    #snr
    trig.snr                     = float(data[6])
    trig.ms_snr                  = float(data[6])

    # ms extras
    ms_start = start
    trig.ms_start_time           = ms_start.seconds
    trig.ms_start_time_ns        = ms_start.nanoseconds
    ms_stop = stop
    trig.ms_stop_time            = ms_stop.seconds
    trig.ms_stop_time_ns         = ms_stop.nanoseconds
    trig.ms_duration             = ms_stop-ms_start
    trig.ms_fhigh                = trig.fhigh
    trig.ms_flow                 = trig.flow
    trig.ms_bandwidth            = trig.ms_fhigh-trig.ms_flow

    # others
    trig.param_two_name = 'numPixels'
    trig.param_two_value = int(data[5])
    trig.param_three_name = 'totPower'
    trig.param_three_value = float(data[7])

    #trig.process_id              = int(float(data[10]))

  # sundries
  if ifo:
    trig.ifo = ifo
  if channel:
    trig.channel = channel

  return trig

# =============================================================================
# Write triggers to file in etg standard form
# =============================================================================

def totrigxml(file, table, program=None, params=[]):

  """
    Write the given lsctables compatible table object to the file object file
    in xml format.
    
    Arguments:

      file: file object
        file object describing output file

      table: lsctable
        glue.ligolw.lsctables compatible table

    Keyword Arguments:

      program: string
        name of generating executable, defaults to self
      params: list
        list of (param,type,value) tuples to append to process_params:table

  """

  xmldoc = ligolw.Document()
  xmlligolw = ligolw.LIGO_LW()
  xmldoc.appendChild(ligolw.LIGO_LW())
  xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessTable))
  xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessParamsTable))

  # append process to table
  if not program:
    program='pylal.dq.dqDataUtils.totrigxml'

  process = ligolw_process.append_process(xmldoc, program=program,\
                                  version=__version__,\
                                  cvs_repository = 'lscsoft',\
                                  cvs_entry_time = __date__)

  ligolw_process.append_process_params(xmldoc, process, params)

  # append trig table to file
  xmldoc.childNodes[-1].appendChild(table)

  # write triggers to file object file
  ligolw_process.set_process_end_time(process)
  with utils.SignalsTrap():
    utils.write_fileobj(xmldoc, file, gz=file.name.endswith('gz'))

# =============================================================================
# Write triggers to text file in etg standard form
# =============================================================================

def totrigfile(file,table,etg,header=True,columns=None):

  """
    Write the lsctables.table to the given file object file, in the standard
    format for the given etg.

    If columns is given as a list, the standard format will be overwritten
    and only the given columns will be processed.
  """

  etg = etg.lower()

  # set columns
  if not columns:
    if re.match('ihope',etg):
      # comma separated values are:
      columns = ['end_time','end_time_ns','ifo','snr','mass1','mass2','mtotal',\
                 'eta','event_duration','template_duration','eff_distance',\
                 'chisq','chisq_dof','bank_chisq','bank_chisq_dof',\
                 'cont_chisq','cont_chisq_dof']

    elif re.match('omegadq',etg) or re.match('omega_dq',etg):
      columns = ['start_time','stop_time','peak_time','flow','fhigh',\
                 'cluster_length','ms_start_time','ms_stop_time','ms_flow',\
                 'ms_fhigh','cluster_size','amplitude','amplitude']

    elif re.match('omegaspectrum',etg):
      columns = ['peak_time','peak_frequency','amplitude','chisq']

    elif re.match('omega',etg) or re.match('wpipe',etg):
      if len(table) and hasattr(table[0], 'param_one_value'):
        columns = ['peak_time','peak_frequency','duration',\
                   'bandwidth','amplitude',\
                   'param_one_value','param_two_value','param_three_value']
      else:
        columns = ['peak_time','peak_frequency','duration',\
                   'bandwidth','amplitude']

    elif re.match('kw',etg.lower()):
      columns = ['peak_time','start_time','stop_time','peak_frequency',\
                 'energy','amplitude','num_pixels','significance','N']

    elif re.match('hacr',etg.lower()):
      columns = ['peak_time','param_one_value','peak_frequency','bandwidth',\
                'duration','param_two_value','snr','param_three_value']

  # set delimiter
  if etg=='ihope':
    d = ','
  else:
    d=' '

  # print header
  if header:
    cols = []
    for c in columns:
      if re.search('param_[a-z]+_value', c) and len(table)>0\
      and hasattr(table[0], c.replace('value','name')):
          cols.append(getattr(table[0], c.replace('value','name')))
      else:
        cols.append(c)
    file.write('%s\n' % d.join(['#']+cols))

  # work out columns
  if len(table)>0:
    t = table[0]
    columnnames = [t for t in table[0].__slots__ if hasattr(table[0], t)]

  # print triggers
  for row in table:
    line = []
    for col in columns:
       if col not in columnnames:
         line.append('-1')
         continue
       entry = ''
       # if ihope, print column
       if re.match('(ihope|hacr)',etg.lower()):
         # HACR default is to have peak_time_ns in seconds, not ns
         if re.match('hacr',etg.lower()) and col=='peak_time_ns':
           entry = str(getattr(row, col))
         else:
           entry = str(getattr(row, col))
       # if not ihope, check for time and print full GPS
       else:
         if col=='peak_time':
           entry = str(row.get_peak())
         elif col=='start_time':
           entry = str(row.get_start())
         elif col=='ms_start_time':
           entry = str(row.get_ms_start())
         elif col=='stop_time':
           entry = str(row.get_stop())
         elif col=='ms_stop_time':
           entry = str(row.get_ms_stop())
         else:
           entry = str(getattr(row, col))

       line.append(entry)

    file.write('%s\n' % d.join(line))

# =============================================================================
# Function to load triggers from xml
# =============================================================================

def fromtrigxml(file,tablename='sngl_inspiral:table',start=None,end=None,\
                columns=None):

  """
    Reads a trigger table from the given table from the xml
    file object file

    Arguments:

      file : file object
   
    Keyword arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      tablename : string
        name of requested trigger table in xml file, defaults to
        'sngl_inspiral:table'
  """

  # set times
  if not start:
    start=0
  if not end:
    end=9999999999

  span = segments.segment(start,end)

  # set columns
  if columns!=None:
    if re.search('sngl_burst', tablename):
      lsctables.SnglBurstTable.loadcolumns = columns
    elif re.search('sngl_inspiral', tablename):
      lsctables.SnglInspiralTable.loadcolumns = columns
    if re.search('multi_burst', tablename):
      lsctables.MultiBurstTable.loadcolumns = columns
    elif re.search('multi_inspiral', tablename):
      lsctables.MultiInspiralTable.loadcolumns = columns

  # set tablename
  if not tablename.endswith(':table'):
    tablename = ':'.join([tablename,'table'])

  # crack open xml file
  xmldoc,digest = utils.load_fileobj(file,gz=file.name.endswith('gz'))
  alltriggers = table.get_table(xmldoc,tablename)

  triggers = table.new_from_template(alltriggers)
  append = triggers.append

  get_time = def_get_time(triggers.tableName)
  # parse triggers in time
  for row in alltriggers:
    if float(get_time(row)) in span:
      append(row)

  # sort table in time
  triggers.sort(key=lambda trig: float(get_time(trig)))

  # reset columns
  if columns:
    type(triggers).loadcolumns = None

  return triggers

# =============================================================================
# Load triggers from text file
# =============================================================================

def fromtrigfile(file,etg,start=None,end=None,ifo=None,channel=None,\
                 tabletype=None, columns=None, virgo=False):

  """
    Reads the file object file containing standard columns for the given etg and
    returns either a corresponding lsctable.

    Arguments:

      file : file object
      etg : [ "ihope" | "kw" | "omega" | "omegadq" ]
        string defining how to read text file.
        "ihope" is read as per the .csv files from ihope_daily. "omegadq"
        is read as per the .txt or .clst produced by the detchar scripts used to
        plot omega online events.

    Keyword arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      tabletype : type
        Specific ligolw table type for output. By default tables will be
        SnglInspiralTable or SnglBurstTable type depending on ETG
  """

  if re.search('omegaspectrum', etg, re.I):
    return fromomegaspectrumfile(file, start=start, end=end, ifo=ifo,\
                           channel=channel, columns=columns)
  elif re.search('omegadq', etg, re.I):
    return fromomegadqfile(file, start=start, end=end, ifo=ifo,\
                           channel=channel,columns=columns)
  elif re.search('omega', etg, re.I):
    return fromomegafile(file, start=start, end=end, ifo=ifo, channel=channel,\
                         columns=columns, virgo=virgo)
  elif re.search('kw', etg, re.I):
    return fromkwfile(file, start=start, end=end, ifo=ifo, channel=channel,\
                      columns=columns)
  elif re.search('hacr', etg, re.I):
    return fromhacrfile(file, start=start, end=end, ifo=ifo, channel=channel,\
                        columns=columns)
  elif re.search('ihope', etg, re.I):
    return fromihopefile(file, start=start, end=end, ifo=ifo, channel=channel,\
                         columns=columns)

# =============================================================================
# Load triggers from a cache
# =============================================================================

def fromLALCache(cache, etg, start=None, end=None, columns=None,\
                 virgo=False, verbose=False, snr=False):

  """
    Extract triggers froa given ETG from all files in a glue.lal.Cache object.
    Returns a glue.ligolw.Table relevant to the given trigger generator etg.
  """

  # set up counter
  if verbose:
    sys.stdout.write("Extracting %s triggers from %d files...     "\
                     % (etg, len(cache)))
    sys.stdout.flush()
    delete = '\b\b\b'
    num = len(cache)/100

  trigs = SnglTriggerTable(etg, columns=columns)

  # load files
  for i,e in enumerate(cache):
    if re.search('(xml|xml.gz)\Z', e.path):
      trigsTmp = fromtrigxml(open(e.path), tablename=trigs.tableName,\
                               start=start, end=end, columns=columns)
    else:
      trigsTmp = fromtrigfile(open(e.path), etg=etg, start=start, end=end,\
                                columns=columns, virgo=virgo)
    # keep only triggers above SNR threshold if requested
    if snr:
      trigs.extend(t for t in trigsTmp if t.snr > snr)
    else:
      trigs.extend(trigsTmp)

    # print verbose message
    if verbose and len(cache)>1:
      progress = int((i+1)/num)
      sys.stdout.write('%s%.2d%%' % (delete, progress))
      sys.stdout.flush()

  if verbose: sys.stdout.write("\n")  
  return trigs

# =============================================================================
# Generate a daily ihope cache 
# =============================================================================

def daily_ihope_cache(start,end,ifo,cluster=None,filetype='xml',cat=0):

  """
    Generates glue.lal.Cache containing CacheEntires for all daily ihope
    INSPIRAL files for given ifo and clustering between start and end time.

    Arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      ifo : [ "H1" | "L1" | "V1" ]
        IFO

    Keyword arguments:
      cluster : [ "unclustered" | "100ms" | "30ms" | "16s" ]
        clustering time in human format
      filetype : [ "xml" | "csv" ]
        file format of desired cache contents
      cat : [ 0 | 1 | 2 | 3 | 4 ]
        int veto category of trigger files requested
  """

  # daily path
  ihope_daily_path = os.path.expanduser('~cbc/ihope_daily')

  # set clustering tag
  if cluster==None or cluster.upper()=='UNCLUSTERED':
    cluster_tag='UNCLUSTERED'
  elif cluster.upper()=='100MS':
    cluster_tag='100MILLISEC_CLUSTERED'
  elif cluster.upper()=='30MS':
    cluster_tag='30MILLISEC_CLUSTERED'
  elif cluster.upper()=='16S':
    cluster_tag='16SEC_CLUSTERED'

  # work out days
  days = gps_day_list(start,end)
  span = segments.segment(start,end)
  cache = LALCache()
  # loop over days gathering files
  for day in days:
    utc = datetime.datetime(*date.XLALGPSToUTC(day)[:6])
    day_path = os.path.join(ihope_daily_path,utc.strftime("%Y%m"),
                                             utc.strftime("%Y%m%d"))

    if filetype=='xml':
      filenames = glob.glob(os.path.join(day_path,
                                      ifo+'-INSPIRAL_'+cluster_tag+'*.xml.gz'))

      for filename in filenames:
        e = LALCacheEntry.from_T050017(filename)
        if span.intersects(e.segment):  cache.append(e)

    elif filetype=='csv':
      csvfile = os.path.join(day_path,ifo+'-'+str(cat)+'-INSPIRAL_'+\
                                      cluster_tag+'.csv')
      if os.path.isfile(csvfile):
        e = LALCacheEntry.from_T050017(csvfile)
        if span.intersects(e.segment):  cache.append(e)

  cache.sort(key=lambda e: e.path)

  return cache

# =============================================================================
# Function to generate an omega online cache
# =============================================================================

def omega_online_cache(start, end, channel, mask='DOWNSELECT',\
                       check_files_exist=False, **kwargs):

  """
    Returns a glue.lal.Cache contatining CacheEntires for all omega online
    trigger files between the given start and end time for the given ifo.
    For S6 triggers are only available for each IFO on it's own site cluster.

    Arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      channel name or ifo : [ "H1" | "L1" | "V1" | "G1:SEI_TCC_STS2x" | ... ]
        IFO
  """

  # format channel
  if re.match('\w\d:', channel):
    ifo, channel = channel.split(':', 1)
  else:
    ifo = channel
    channel = None

  cache = LALCache()

  # verify host
  host = getfqdn()
  ifo_host = { 'G1':'atlas', 'H1':'ligo-wa', 'H2':'ligo-wa', 'L1':'ligo-la'}
  if (not kwargs.has_key('directory') and not re.search(ifo_host[ifo],host)):
    print >>sys.stderr, "Error: Omega online files are not available for "+\
                        "IFO=%s on this host." % ifo
    return cache

  span = segments.segment(start,end)
  if ifo == 'G1':
    if channel:
      kwargs.setdefault('directory', '/home/omega/online/%s_%s/segments' % (ifo, channel))
    else:
      kwargs.setdefault('directory', '/home/omega/online/%s/segments' % ifo)
    kwargs.setdefault('epoch', 0)
  else:
    kwargs.setdefault('directory',\
                      '/home/omega/online/%s/archive/S6/segments' % ifo)
    kwargs.setdefault('epoch', 931211808)
  kwargs.setdefault('duration', 64)
  kwargs.setdefault('overlap', 8)

  # optimise
  append       = cache.append
  splitext     = os.path.splitext
  isfile   = os.path.isfile
  intersects   = span.intersects
  segment      = segments.segment
  from_T050017 = LALCacheEntry.from_T050017
  basedir      = kwargs['directory']
  basetime     = kwargs['epoch']
  triglength   = kwargs['duration']
  overlap      = kwargs['overlap']

  # get times
  start_time = int(start-math.fmod(start-basetime,triglength-overlap))
  t = start_time

  # loop over time segments constructing file paths and appending to the cache
  while t<end:
    if ifo == 'G1':
      trigfile = '%s/%d/%d-%d/%s-OMEGA_TRIGGERS_%s-%.10d-%d.txt'\
          % (basedir, t/100000, t, t+triglength, ifo, mask, t, triglength)
    else:
      trigfile = '%s/%.10d-%10.d/%s-OMEGA_TRIGGERS_%s-%.10d-%d.txt'\
          % (basedir, t, t+triglength, ifo, mask, t, triglength)
    if intersects(segment(t, t+triglength))\
    and (not check_files_exist or isfile(trigfile)):
      append(from_T050017(trigfile))
    t+=triglength-overlap

  cache.sort(key=lambda e: e.path)

  return cache

# =============================================================================
# Function to generate an omega spectrum online cache
# =============================================================================

def omega_spectrum_online_cache(start, end, ifo, **kwargs):

  """
    Returns a glue.lal.Cache contatining CacheEntires for all omega online
    trigger files between the given start and end time for the given ifo.
    For S6 triggers are only available for each IFO on it's own site cluster.

    Arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      ifo : [ "H1" | "L1" | "V1" ]
        IFO
  """

  return omega_online_cache(start, end, ifo, mask='SPECTRUM', **kwargs)

# =============================================================================
# DetChar 'omegadq' cache
# =============================================================================

def omega_dq_cache(start,end,ifo):

  """
    Arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      ifo : [ "H1" | "L1" | "V1" ]
        IFO
  """

  # verify host
  host = getfqdn()
  ifo_host = {'H1':'ligo-wa','H2':'ligo-wa','L1':'ligo-la'}
  if not re.search(ifo_host[ifo],host):
    print >>sys.stderr, "Error: OmegaClustered files are not available for "+\
                        "IFO="+ifo+" on this host."
    return LALCache()

  cache = LALCache()
  basedir = os.path.expanduser('~detchar/public_html/S6/glitch/Wdata')

  # work out days
  days = gps_day_list(start,end)
  span = segments.segment(start,end)
  triglength=86400
  cache = LALCache()

  for day in days:
    dirstart = day
    dirend   = day+triglength
    dirpath  = os.path.join(basedir,'%s_%s' % (dirstart,dirend))
    trigfile = os.path.join(dirpath,'clusters.txt')
    if os.path.isfile(trigfile):

      e = LALCacheEntry(ifo,'OMEGADQ',segments.segment(dirstart,dirend),\
                        os.path.realpath(trigfile))
      if span.intersects(e.segment):  cache.append(e)

  cache.sort(key=lambda e: e.path)

  return cache

# ==============================================================================
# Class for KWCacheEntry  
# ==============================================================================

class KWCacheEntry(LALCacheEntry):

  _regex = re.compile(r"\A\s*(?P<obs>\S+)\s+(?P<dsc>\S+)\s*\Z")

  def from_KWfilename(cls, url, coltype = LIGOTimeGPS):
    """
    Parse a URL in the style of KW filenames into a FrameCacheEntry.
    The KW file name format is, essentially,

    /path/to/start_end/observatory_description.txt


    """

    try:
      head,tail = os.path.split(url)
      observatory,description = re.split('_',os.path.splitext(tail)[0],\
                                         maxsplit=1)
      observatory = observatory
      start,end = [coltype(t) for t in os.path.basename(head).split('_')]
      duration = end-start

      segment = segments.segment(start,end)

    except:
      raise ValueError, "could not convert %s to KWCacheEntry" % repr(url)

    return cls(observatory, description, segment, url)

  from_KWfilename = classmethod(from_KWfilename)

  def from_T050017(cls,url,coltype = LIGOTimeGPS):
    """
    Redirects to from_KWfilename for KWCacheEntry objects due to KW not
    following T50017-00 conventions.
    """
    return KWCacheEntry.from_KWfilename(url,coltype=coltype)

  from_T050017 = classmethod(from_T050017)

# =============================================================================
# Function to generate a KW DARM_ERR cache
# =============================================================================

def kw_cache(start, end, channel='H1:LSC-DARM_ERR', frequency=None):

  """
    Returns a list of KW trigger files between the given start and end
    time for the given ifo. For S6 triggers are only available for each IFO on
    it's own site cluster.

    Arguments:

      start : [ float | int | LIGOTimeGPS ]
        GPS start time of requested period
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of requested period
      ifo : [ "H1" | "L1" | "V1" ]
        IFO
  """

  # format channel
  if re.match('\w\d:', channel):
    ifo, channel = channel.split(':', 1)
    #_cchar_regex.sub('_', channel)
  else:
    raise ValueError("Please give channel in the form \"IFO:CHANNEL-NAME\"")
  
  # verify host
  host = getfqdn()
  ifo_host = {'H0':'(ligo-wa|ligo\.)', 'H1':'(ligo-wa|ligo\.)',\
              'H2':'ligo-wa', 'L0':'(ligo-la|ligo\.)', 'L1':'(ligo-la|ligo\.)',\
              'V1':'ligo\.'}
  if not re.search(ifo_host[ifo],host):
    print >>sys.stderr, "Warning: KW files are not available for "+\
                        "IFO="+ifo+" on this host."
    return LALCache()

  cache = LALCache()
  if ifo == 'V1':
    basedir = os.path.expanduser('~mabizoua/public_html/KW')
  elif re.search('ligo\.', host):
    basedir = os.path.expanduser('~lindy/public_html/triggers/s6-merged')
  else:
    basedir = os.path.expanduser('~lindy/public_html/triggers/s6')

  # times are numbere from a given start, which for S6 is:
  if ifo=='V1':
    base = LIGOTimeGPS(938736015)
  else:
    base = LIGOTimeGPS(938736000)
  triglength = 86400

  # construct span
  span = segments.segment(start, end)

  # get fbin
  if not frequency:
    f = '*_*'
  elif frequency.lower() == 'low':
    f = '32_2048'
  elif frequency.lower() == 'high':
    f = '1024_4096'
  elif not isinstance(frequency ,str) and len(frequency)==2:
    f = '%d_%d' % tuple(frequency)
  else:
    f = '*_*'

  # get times
  start_time = int(start-math.fmod(start-base,triglength))
  t = start_time

  # loop over time segments constructing file paths and appending to the cache
  while t<end:
    dirstart = str(t)
    dirend   = str(t+triglength)
    dirpath  = os.path.join(basedir,dirstart+'_'+dirend)
    trigfile = '%s/%s_%s_%s.trg' % (dirpath, ifo, channel, f)
    if '*' in trigfile:
      trigfiles = glob.glob(trigfile)
    else:
      trigfiles = [trigfile]
    for trigfile in trigfiles:
      if os.path.isfile(trigfile):
        e = KWCacheEntry.from_KWfilename(trigfile)
        if span.intersects(e.segment):
          cache.append(e)

    t+=triglength

  cache.sort(key=lambda e: e.path) 

  return cache

# ==============================================================================
# Function to construct list of days from start and end times
# ==============================================================================

def gps_day_list(start,end):

  """
    This script will construct a list of LIGOTimeGPS days encompassing the start and end GPS times.
  """

  start=LIGOTimeGPS(start)

  start_d = date.XLALUTCToGPS(datetime.datetime(*date.XLALGPSToUTC(start)[:6])\
                                  .replace(hour=0,minute=0,second=0)\
                                  .timetuple())
  days = []
  day = start_d
  while day<=end:
    days.append(day)
    day+=86400

  return days

# =============================================================================
# Function to cluster triggers
# =============================================================================

def cluster(triggers,params=[('time',1)],rank='snr'):

  """
    Cluster the lsctable triggers in the each of the pairs (column,width),
    using the rank column.

    Arguments:

      triggers: glue.ligowl.Table
        Table containing trigger columns for clustering

      params : list
        List object containing (column,width) pairs. Clustering is nested in
        order

      rank: string
        Column by which to rank clusters
      test
  """

  outtrigs = table.new_from_template(triggers)

  j = 0

  cols = [p[0] for p in params]
  coldata = dict((p, get_column(triggers, p)) for p in cols+[rank])
  # need bandwidth and duration for all burst triggers
  if _burst_regex.search(triggers.tableName):
    coldata['stop_time'] = get_column(triggers, 'stop_time') +\
                           get_column(triggers, 'stop_time_ns')*1e-9
    coldata['start_time'] = get_column(triggers, 'start_time') +\
                            get_column(triggers, 'start_time_ns')*1e-9
    coldata['flow'] = get_column(triggers, 'flow')
    coldata['fhigh'] = get_column(triggers, 'fhigh')
  # need time for time clustering
  elif 'time' in cols:
    coldata['stop_time'] = coldata['time']
    coldata['start_time'] = coldata['time']

  for key in coldata.keys():
    coldata[key] = coldata[key].astype(float)

  # for each parameter break the clusters generated using the previous
  # parameter into smaller clusters by sorting triggers and clustering
  # when all parameters have been used, pick the loudest in each cluster

  clusters = [range(len(triggers))]

  while j < len(params):

    col,width = params[j]

    newclusters = []

    for k,subcluster in enumerate(clusters):

      # sort triggers incluster parameter
      subcluster.sort(key=lambda i: coldata[col][i])
      subsubcluster = []

      for i in subcluster:

        # get value of param
        if col=='time':
	  valueStop = coldata['stop_time'][i]
          valueStart = coldata['start_time'][i]
        elif col=='peak_frequency':
          valueStop = coldata['fhigh'][i]
          valueStart = coldata['flow'][i]
        else:
          valueStop = coldata[col][i]
          valueStart = valueStop

        # if subcluster is empty, simply add the first trigger
        if not subsubcluster:
          subsubcluster = [i]
          prevStop = valueStop
          prevStart = valueStart
          continue

        # if current trig is inside width, append to cluster
        if (valueStart-prevStop)<width:
          subsubcluster.append(i)

        # if not the subcluster is complete, append it to list and start again
        else:
          newclusters.append(subsubcluster)
          subsubcluster=[i]

        prevStart = valueStart
        prevStop = valueStop

      # append final subsubcluster
      newclusters.append(subsubcluster)

    clusters = copy.deepcopy(newclusters)
    j += 1

  # process clusters
  for cluster in clusters:
    if len(cluster)==1:
      outtrigs.append(copy.deepcopy(triggers[cluster[0]]))
    elif len(cluster) > 1:
      carray = numpy.asarray(cluster)
      cluster.sort(key=lambda i: coldata[rank][i], reverse=True)
      t = copy.deepcopy(triggers[cluster[0]])    
      # reset burst params for a clustered event
      if _burst_regex.search(triggers.tableName):
        # record most significant trigger
        t.ms_start_time = t.start_time
        t.ms_start_time_ns = t.start_time_ns
        t.ms_stop_time = t.stop_time
        t.ms_stop_time_ns = t.stop_time_ns
        t.ms_duration = t.duration
        t.ms_bandwidth = t.bandwidth
        t.ms_flow = t.flow
        t.ms_fhigh = t.fhigh
        t.ms_snr = t.snr
        # record cluster
        start = LIGOTimeGPS(min(coldata['start_time'][carray]))
        t.start_time = start.seconds
        t.start_time_ns = start.nanoseconds
        stop = LIGOTimeGPS(max(coldata['stop_time'][carray]))
        t.stop_time = stop.seconds
        t.stop_time_ns = stop.nanoseconds
        t.duration = float(t.get_stop()-t.get_start())
        t.flow = min(coldata['flow'][carray])
        t.fhigh = max(coldata['fhigh'][carray])
        t.bandwidth = t.fhigh-t.flow
        t.tfvolume = t.bandwidth * t.duration
      outtrigs.append(t)

  # resort trigs in first parameter
  outtrigs.sort(key=lambda t: get(t, cols[0]))

  return outtrigs

# =============================================================================
# Compute trigger auto-correlation
# =============================================================================

def autocorr(triggers,column='time',timeStep=0.02,timeRange=60):

  """
    Compute autocorrelation of lsctable triggers in the each of the pairs
    (column,width), using the rank column.

    Arguments:

      triggers: glue.ligowl.Table
        Table containing trigger columns for clustering

      column:
        On which trigger column to auto-correlate, almost always time

      timeStep:
        Step (bin width) for the autocorrelation

      timeRange:
        Longest time to consider for autocorrelation
        
  """

  # time sort triggers before proceeding
  get_time = def_get_time(triggers.tableName)
  triggers.sort(key=lambda trig: get_time(trig))

  previousTimes = []
  histEdges = numpy.arange(timeStep,timeRange,timeStep);
  delayHist = numpy.zeros(int(math.ceil(timeRange/timeStep)))

  for trig in triggers:
    curTime = trig.peak_time + 1e-9*trig.peak_time_ns
    # remove previous times which are beyond the considered timeRange
    while len(previousTimes) > 0 and curTime - previousTimes[0] > timeRange:
      previousTimes.pop()
    for t in previousTimes:
      pos = int(math.floor((curTime - t)/timeStep))
      if pos < len(delayHist):
        delayHist[pos] += 1
    previousTimes.append(curTime)
  
  delayHistFFT = numpy.abs(numpy.fft.fft(delayHist))
  freqBins = numpy.fft.fftfreq(len(delayHist), d=timeStep)

  return delayHistFFT, freqBins, delayHist, histEdges

# =============================================================================
# Get coincidences between two tables
# =============================================================================

def get_coincs(table1, table2, dt=1, returnsegs=False, timeshift=0):

  """
    Returns the table of those entries in table1 whose time is within +-dt of
    an entry in table2.
  """

  t1 = get_column(table1, 'time')
  t2 = get_column(table2, 'time')

  coincsegs  = segments.segmentlist(segments.segment(t-dt+timeshift, t+dt+timeshift) for t in t2)\
                   .coalesce()
  coincsegs.sort()
  coinctrigs = table.new_from_template(table1)
  coinctrigs.extend(t for i,t in enumerate(table1) if t1[i] in coincsegs)

  if returnsegs:
    return coinctrigs,coincsegs
  else:
    return coinctrigs

# =============================================================================
# Get number of coincidences between two tables
# =============================================================================

def get_number_coincs(table1, table2, dt=1, timeshift=0, tabletype='trigger'):

  """
    Returns the numbers of entries in table1 whose time is within +-dt of
    an entry in table2.
  """

  # t1 = get_column(table1, 'time')
  # t2 = get_column(table2, 'time')

  # coincsegs  = segments.segmentlist(segments.segment(t-dt+timeshift, t+dt+timeshift) for t in t2)\
  #                  .coalesce()
  # coincsegs.sort()
  # coinctrigs = table.new_from_template(table1)
  # coinctrigs.extend(t for i,t in enumerate(table1) if t1[i] in coincsegs)

  # return len(coinctrigs)

  if tabletype == 'trigger':
    time1 = get_column(table1, 'time')
    time2 = get_column(table2, 'time')
  elif tabletype == 'time':
    time1 = table1
    time2 = table2
  else:
    raise ValueError("Unrecognized table type for coincidence number: %s" % tabletype)

  time1.sort()
  time2.sort()

  i2 = 0
  ncoinc = 0
  for i1,t1 in enumerate(time1):
    while  i2 < len(time2) and t1 > time2[i2]+timeshift+dt:
      i2 = i2 + 1
    if i2 >= len(time2):
      break
    if t1 >= time2[i2]+timeshift-dt:
      ncoinc = ncoinc + 1

  return ncoinc

# ==============================================================================
# Calculate poisson significance of time coincidences
# ==============================================================================

def coinc_significance_times(gwtrigtime, auxtrigtime, window=1, livetime=None):


  # get livetime
  if not livetime:
    start    = min(gwtrigtime)
    end      = max(gwtrigtime)
    livetime = end-start

  # calculate probability of a GW trigger falling within the window
  gwprob = len(gwtrigtime) * 2.0 * float(window) / float(livetime)

  # calculate mean of Poisson distribution
  mu = gwprob * len(auxtrigtime)

  # get coincidences
  ncoinc = get_number_coincs(gwtrigtime, auxtrigtime, dt=window, tabletype='time')

  g = special.gammainc(ncoinc, mu)

  # if no coincidences, set significance to zero
  if ncoinc<1:
    significance = 0
  # if significance would blow up, use other formula (ref. hveto_significance.m)
  elif g == 0:
    significance = -ncoinc * math.log10(mu) + \
                   mu * math.log10(math.exp(1)) +\
                   special.gammaln(ncoinc + 1) / math.log(10)
  # otherwise use the standard formula
  else:
    significance = -math.log(g, 10)
  
  return significance


# ==============================================================================
# Calculate poisson significance of trigger coincidences
# ==============================================================================

def coinc_significance(gwtriggers, auxtriggers, window=1, livetime=None,\
                       returnsegs=False):

  get_time = def_get_time(gwtriggers.tableName)
  aux_get_time = def_get_time(auxtriggers.tableName)

  # get livetime
  if not livetime:
    start    = min([get_time(t) for t in gwtriggers])
    end      = max([get_time(t) for t in gwtriggers])
    livetime = end-start

  # calculate probability of a GW trigger falling within the window
  gwprob = len(gwtriggers) * 2.0 * float(window) / float(livetime)

  # calculate mean of Poisson distribution
  mu = gwprob * len(auxtriggers)

  # get coincidences
  if returnsegs:
    coinctriggers,coincsegs = get_coincs(gwtriggers, auxtriggers, dt=window,\
                                           returnsegs=True)
    ncoinc = len(coinctriggers)
  else:
    ncoinc = get_number_coincs(gwtriggers, auxtriggers, dt=window)

  g = special.gammainc(ncoinc, mu)

  # if no coincidences, set significance to zero
  if ncoinc<1:
    significance = 0
  # if significance would blow up, use other formula (ref. hveto_significance.m)
  elif g == 0:
    significance = -ncoinc * math.log10(mu) + \
                   mu * math.log10(math.exp(1)) +\
                   special.gammaln(ncoinc + 1) / math.log(10)
  # otherwise use the standard formula
  else:
    significance = -math.log(g, 10)

  if returnsegs:
    return significance,coincsegs
  else:
    return significance

# =============================================================================
# Extract column from generic table
# =============================================================================

def get_column(lsctable, column):

  """
    Extract column from the given glue.ligolw.table lsctable as numpy array.
    Tries to use a 'get_col() function if available, otherwise uses
    getColumnByName(), treating 'time' as a special case for the known
    Burst/Inspiral/Ringdown tables.
  """
 
  # format column
  column = str(column).lower()
  obj_type = str(type(lsctable))

  # if there's a 'get_' function, use it
  if hasattr(lsctable, 'get_%s' % column):
    return numpy.asarray(getattr(lsctable, 'get_%s' % column)())

  # treat 'time' as a special case
  elif column == 'time' and _trig_regex.search(obj_type):
    if _burst_regex.search(obj_type):
      tcol = 'peak_time'
    if _cbc_regex.search(obj_type):
      tcol = 'end_time'
    if _ring_regex.search(obj_type):
      tcol = 'start_time'
    return numpy.asarray(lsctable.getColumnByName(tcol)) + \
           numpy.asarray(lsctable.getColumnByName('%s_ns' % tcol))*10**-9

  return numpy.asarray(lsctable.getColumnByName(column))

def get(self, parameter):

  """
    Extract parameter from given ligolw table row object. 
  """

  # format
  parameter = parameter.lower()
  obj_type = str(type(self))

  # if there's a 'get_' function, use it
  if hasattr(self, 'get_%s' % parameter):
    return getattr(self, 'get_%s' % parameter)()

  # treat 'time' as a special case
  elif parameter == 'time' and _trig_regex.search(obj_type):
    if _burst_regex.search(obj_type):
      tcol = 'peak_time'
    elif _cbc_regex.search(obj_type):
      tcol = 'end_time'
    elif _ring_regex.search(obj_type):
      tcol = 'start_time'
    return getattr(self, tcol)+getattr(self, '%s_ns' % tcol)*10**-9

  else:
   return getattr(self, parameter)

# =============================================================================
# Veto triggers from a ligolw table
# =============================================================================

def veto(self, seglist, inverse=False):

  """
    Returns a ligolw table of those triggers outwith the given seglist.
    If inverse=True is given, the opposite is returned, i.e. those triggers
    within the given seglist.
  """

  get_time = def_get_time(self.tableName)

  keep = table.new_from_template(self)
  if inverse:
    keep.extend(t for t in self if float(get_time(t)) in seglist)
  else:
    keep.extend(t for t in self if float(get_time(t)) not in seglist)

  return keep

def vetoed(self, seglist):

  """
    Returns the opposite of veto, i.e. those triggers that lie within the given
    seglist.
  """

  return veto(self, seglist, inverse=True)

# =============================================================================
# read triggers from file
# =============================================================================

def fromomegafile(fname, start=None, end=None, ifo=None, channel=None,\
                  columns=None, virgo=False):

  """
    Load triggers from an Omega format text file into a SnglBurstTable object.
    Use start and end to restrict the returned triggers, and give ifo and
    channel to fill those columns in the table.

    If columns is given as a list, only those columns in the table will be
    filled. This is advisable to speed up future operations on this table.

    Arguments :

      fname : file or str
        file object or filename path to read with numpy.loadtext

    Keyword arguments :

      start : float
        minimum peak time for returned triggers
      end : float
        maximum peak time for returned triggers
      ifo : str
        name of IFO to fill in table
      channel : str
        name of channel to fill in table
      columns : iterable
        list of columnnames to populate in table
  """

  # set columns
  if columns==None: columns = lsctables.SnglBurst.__slots__
  if start or end:
    if not start:
      start = 0
    if not end:
      end   = numpy.inf
    span = segments.segment(start, end)
    if 'peak_time' not in columns: columns.append('peak_time')
    if 'peak_time_ns' not in columns: columns.append('peak_time_ns')
    check_time = True
  else:
    check_time = False

  if 'snr' in columns and not 'amplitude' in columns:
    columns.append('amplitude')

  # generate table
  out = SnglTriggerTable('omega', columns=columns)

  # force filename not file object
  if hasattr(fname, 'readline'):
    fh = fname
  else:
    fh = open(fname, 'r')

  dat = loadtxt(fh)

  if not hasattr(fname, 'readline'):
    fh.close()

  if numpy.shape(dat) == (0,):
    return out

  if virgo:
    start, stop, peak, freq, bandwidth, cln, cle, snr = dat
    duration = stop-start
    amplitude = snr**2/2
    omega_clusters = False
  if len(dat)==11:
    peak, freq, duration, bandwidth, amplitude, cls, cle, cln, av_freq, av_bandwidth, err_freq = dat
    omega_clusters = True
  else:
    if len(dat)==8:
      peak, freq, duration, bandwidth, amplitude, cls, cle, cln = dat
      omega_clusters = True
    elif len(dat)==5:
      peak, freq, duration, bandwidth, amplitude = dat
      omega_clusters = False
    else:
      raise ValueError("Wrong number of columns in omega format file. "\
                         "Cannot read.")
    av_freq = freq
    av_bandwidth = bandwidth
    err_freq = av_bandwidth/av_freq

  numtrigs = len(peak)
  attr_map = dict()

  if 'start_time' in columns or 'start_time_ns' in columns:
    start = map(LIGOTimeGPS, peak - duration/2)
    attr_map['start_time'], attr_map['start_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in start])
  if 'stop_time' in columns or 'stop_time_ns' in columns:
    stop = map(LIGOTimeGPS, peak + duration/2)
    attr_map['stop_time'], attr_map['stop_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in stop])
  if 'peak_time' in columns or 'peak_time_ns' in columns:
    peak = map(LIGOTimeGPS, peak)
    attr_map['peak_time'], attr_map['peak_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in peak])

  if 'ms_start_time' in columns or 'ms_start_time_ns' in columns:
    ms_start = map(LIGOTimeGPS, peak-duration/2)
    attr_map['ms_start_time'], attr_map['ms_start_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in ms_start])
  if 'ms_stop_time' in columns or 'ms_stop_time_ns' in columns:
    ms_stop = map(LIGOTimeGPS, peak+duration/2)
    attr_map['ms_stop_time'], attr_map['ms_stop_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in ms_stop])

  if 'central_freq' in columns:   attr_map['central_freq']   = freq
  if 'peak_frequency' in columns: attr_map['peak_frequency'] = av_freq
  if 'peak_frequency_error' in columns: attr_map['peak_frequency_error'] = err_freq
  if 'bandwidth' in columns:      attr_map['bandwidth']      = av_bandwidth
  if 'ms_bandwidth' in columns:   attr_map['ms_bandwidth']   = bandwidth
  if 'flow' in columns:           attr_map['flow']           = freq-bandwidth/2
  if 'fhigh' in columns:          attr_map['fhigh']          = freq+bandwidth/2
  if 'ms_flow' in columns:        attr_map['ms_flow']        = freq-bandwidth/2
  if 'ms_fhigh' in columns:       attr_map['ms_fhigh']       = freq+bandwidth/2

  if 'duration' in columns:       attr_map['duration']       = duration
  if 'ms_duration' in columns:    attr_map['ms_duration']    = duration
  if 'snr' in columns:            attr_map['snr']         = (2*amplitude)**(1/2)
  if 'ms_snr' in columns:         attr_map['ms_snr']      = (2*amplitude)**(1/2)

  if 'cluster_size' in columns or 'param_one_value' in columns:
    attr_map['param_one_name'] = ['cluster_size'] * numtrigs
    if omega_clusters:
      attr_map['param_one_value'] = cls
    else:
      attr_map['param_one_value'] = [numpy.NaN] * numtrigs
  if 'cluster_norm_energy' in columns or 'param_two_value' in columns:
    attr_map['param_two_name'] = ['cluster_norm_energy'] * numtrigs
    if omega_clusters:
      attr_map['param_two_value'] = cls
    else:
      attr_map['param_two_value'] = [numpy.NaN] * numtrigs
  if 'cluster_size' in columns or 'param_three_value' in columns:
    attr_map['param_three_name'] = ['cluster_number'] * numtrigs
    if omega_clusters:
      attr_map['param_three_value'] = cls
    else:
      attr_map['param_three_value'] = [numpy.NaN] * numtrigs

  cols   = attr_map.keys()
  append = out.append
  for i in range(numtrigs):
    t = lsctables.SnglBurst()
    for c in cols: setattr(t, c, attr_map[c][i])
    if not check_time or (check_time and float(t.get_peak()) in span):
      append(t)
  
  return out

def fromkwfile(fname, start=None, end=None, ifo=None, channel=None,\
               columns=None):

  """
    Load triggers from a KW format text file into a SnglBurstTable object.
    Use start and end to restrict the returned triggers, and give ifo and
    channel to fill those columns in the table.

    If columns is given as a list, only those columns in the table will be
    filled. This is advisable to speed up future operations on this table.

    Arguments :

      fname : file or str
        file object or filename path to read with numpy.loadtext

    Keyword arguments :

      start : float
        minimum peak time for returned triggers
      end : float
        maximum peak time for returned triggers
      ifo : str
        name of IFO to fill in table
      channel : str
        name of channel to fill in table
      columns : iterable
        list of columnnames to populate in table
  """

  # set columns
  if columns==None: columns = lsctables.SnglBurst.__slots__
  if start or end:
    if not start:
      start = 0
    if not end:
      end   = numpy.inf
    span = segments.segment(start, end)
    if 'peak_time' not in columns: columns.append('peak_time')
    if 'peak_time_ns' not in columns: columns.append('peak_time_ns')
    check_time = True
  else:
    check_time = False

  # generate table
  out = SnglTriggerTable('kw', columns=columns)

  # force filename not file object
  if hasattr(fname, 'readline'):
    fh = fname
  else:
    fh = open(fname, 'r')

  # load data from file
  dat = loadtxt(fh, usecols=[0,1,2,3,4,5,6,7])

  # close file if we opened it
  if not hasattr(fname, 'readline'):
    fh.close()

  if numpy.shape(dat) == (0,):
    return out

  if len(dat)==8:
    st, stop, peak, freq, energy, amplitude, n_pix, sig = dat
  else:
    raise ValueError("Wrong number of columns in KW format file. "\
                     "Cannot read.")

  numtrigs = len(peak)

  attr_map = dict()

  if 'duration' in columns:       attr_map['duration']       = stop-st
  if 'ms_duration' in columns:    attr_map['ms_duration']    = stop-st

  if 'start_time' in columns or 'start_time_ns' in columns:
    start = map(LIGOTimeGPS, st)
    attr_map['start_time'], attr_map['start_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in start])
  if 'stop_time' in columns or 'stop_time_ns' in columns:
    stop = map(LIGOTimeGPS, stop)
    attr_map['stop_time'], attr_map['stop_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in stop])
  if 'peak_time' in columns or 'peak_time_ns' in columns:
    peak = map(LIGOTimeGPS, peak)
    attr_map['peak_time'], attr_map['peak_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in peak])

  if 'ms_start_time' in columns or 'ms_start_time_ns' in columns:
    ms_start = map(LIGOTimeGPS, st)
    attr_map['ms_start_time'], attr_map['ms_start_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in ms_start])
  if 'ms_stop_time' in columns or 'ms_stop_time_ns' in columns:
    ms_stop = map(LIGOTimeGPS, stop)
    attr_map['ms_stop_time'], attr_map['ms_stop_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in ms_stop])

  if 'central_freq' in columns:   attr_map['central_freq']   = freq
  if 'peak_frequency' in columns: attr_map['peak_frequency'] = freq
  if 'bandwidth' in columns:      attr_map['bandwidth'] = numpy.zeros(len(freq))
  if 'ms_bandwidth' in columns: attr_map['ms_bandwidth'] = attr_map['bandwidth']
  if 'flow' in columns:           attr_map['flow']           = freq
  if 'fhigh' in columns:          attr_map['fhigh']          = freq
  if 'ms_flow' in columns:        attr_map['ms_flow']        = freq
  if 'ms_fhigh' in columns:       attr_map['ms_fhigh']       = freq

  if 'duration' in columns:       attr_map['duration']       = stop-st
  if 'ms_duration' in columns:    attr_map['ms_duration']    = stop-st
  if 'snr' in columns:         attr_map['snr']     = (amplitude-n_pix)**(1/2)
  if 'ms_snr' in columns:      attr_map['ms_snr']  = (amplitude-n_pix)**(1/2)
  if 'confidence' in columns:      attr_map['confidence']  = sig
  
  if 'n_pix' in columns or 'param_one_value' in columns:
    attr_map['param_one_name'] = ['n_pix'] * numtrigs
    attr_map['param_one_value'] = n_pix
  """
  if 'significance' in columns or 'param_two_value' in columns:
    attr_map['param_two_name'] = ['significance'] * numtrigs
    attr_map['param_two_value'] = sig
  """

  cols   = attr_map.keys()
  append = out.append
  for i in range(numtrigs):
    t = lsctables.SnglBurst()
    for c in cols: setattr(t, c, attr_map[c][i])
    if ifo!=None:
      t.ifo = ifo
    if channel!=None:
      t.channel = channel
    if not check_time or (check_time and t.get_peak() in span):
      append(t)
  
  return out

def fromomegaspectrumfile(fname, start=None, end=None, ifo=None, channel=None,\
                          columns=None):

  """
    Load triggers from an OmegaSpectrum format text file into a SnglBurstTable
    object.
    Use start and end to restrict the returned triggers, and give ifo and
    channel to fill those columns in the table.

    If columns is given as a list, only those columns in the table will be
    filled. This is advisable to speed up future operations on this table.

    Arguments :

      fname : file or str
        file object or filename path to read with numpy.loadtext

    Keyword arguments :

      start : float
        minimum peak time for returned triggers
      end : float
        maximum peak time for returned triggers
      ifo : str
        name of IFO to fill in table
      channel : str
        name of channel to fill in table
      columns : iterable
        list of columnnames to populate in table
  """

  # set columns
  if columns==None: columns = lsctables.SnglBurst.__slots__
  if start or end:
    if not start:
      start = 0
    if not end:
      end   = numpy.inf
    span = segments.segment(start, end)
    if 'peak_time' not in columns: columns.append('peak_time')
    if 'peak_time_ns' not in columns: columns.append('peak_time_ns')
    check_time = True
  else:
    check_time = False

  # generate table
  out = SnglTriggerTable('omegaspectrum', columns=columns)

  # force filename not file object
  if hasattr(fname, 'readline'):
    fh = fname
  else:
    fh = open(fname, 'r')

  # load data from file
  dat = loadtxt(fh)

  # close file if we opened it
  if not hasattr(fname, 'readline'):
    fh.close()

  if numpy.shape(dat) == (0,):
    return out

  if len(dat)==3:
    peak, freq, amplitude = dat
  elif len(dat)==4:
    peak, freq, amplitude, chisq = dat
  else:
    raise ValueError("Wrong number of columns in omega spectrum format file. "\
                     "Cannot read.")
  numtrigs = len(peak)

  attr_map = dict()

  if 'peak_time' in columns or 'peak_time_ns' in columns:
    peak = map(LIGOTimeGPS, peak)
    attr_map['peak_time'], attr_map['peak_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in peak])

  if 'central_freq' in columns:   attr_map['central_freq']   = freq
  if 'peak_frequency' in columns: attr_map['peak_frequency'] = freq
  if 'amplitude' in columns:      attr_map['amplitude'] = amplitude
  if 'snr' in columns:            attr_map['snr'] = amplitude**(1/2)
  if 'chisq' in columns:          attr_map['chisq'] = chisq**(0.25)/attr_map['snr'] 

  cols   = attr_map.keys()
  append = out.append
  for i in range(numtrigs):
    t = lsctables.SnglBurst()
    for c in cols: setattr(t, c, attr_map[c][i])
    if ifo!=None:
      t.ifo = ifo
    if channel!=None:
      t.channel = channel
    if not check_time or (check_time and t.get_peak() in span):
      append(t)
  
  return out

def fromomegadqfile(fname, start=None, end=None, ifo=None, channel=None,\
                    columns=None):

  """
    Load triggers from an OmegaDQ format text file into a SnglBurstTable object.
    Use start and end to restrict the returned triggers, and give ifo and
    channel to fill those columns in the table.

    If columns is given as a list, only those columns in the table will be
    filled. This is advisable to speed up future operations on this table.

    Arguments :

      fname : file or str
        file object or filename path to read with numpy.loadtext

    Keyword arguments :

      start : float
        minimum peak time for returned triggers
      end : float
        maximum peak time for returned triggers
      ifo : str
        name of IFO to fill in table
      channel : str
        name of channel to fill in table
      columns : iterable
        list of columnnames to populate in table
  """

  # set columns
  if columns==None: columns = lsctables.SnglBurst.__slots__
  if start or end:
    if not start:
      start = 0
    if not end:
      end   = numpy.inf
    span = segments.segment(start, end)
    if 'peak_time' not in columns: columns.append('peak_time')
    if 'peak_time_ns' not in columns: columns.append('peak_time_ns')
    check_time = True
  else:
    check_time = False

  # generate table
  out = SnglTriggerTable('omegadq', columns=columns)

  # force filename not file object
  if hasattr(fname, 'readline'):
    fh = fname
  else:
    fh = open(fname, 'r')

  # load data from file
  dat = loadtxt(fh)

  # close file if we opened it
  if not hasattr(fname, 'readline'):
    fh.close()

  if numpy.shape(dat) == (0,):
    return out

  if len(dat)==13:
    st, stop, peak, flow, fhigh, nev, ms_start, ms_stop, ms_flow, ms_fhigh,\
    cls, cle, ms_cle = dat
  else:
    raise ValueError("Wrong number of columns in OmegaDQ format file. "\
                     "Cannot read.")
  numtrigs = len(peak)

  attr_map = dict()

  if 'duration' in columns:       attr_map['duration']       = stop-st
  if 'ms_duration' in columns:    attr_map['ms_duration']    = ms_stop-ms_start

  if 'start_time' in columns or 'start_time_ns' in columns:
    start = map(LIGOTimeGPS, st)
    attr_map['start_time'], attr_map['start_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in start])
  if 'stop_time' in columns or 'stop_time_ns' in columns:
    stop = map(LIGOTimeGPS, stop)
    attr_map['stop_time'], attr_map['stop_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in stop])
  if 'peak_time' in columns or 'peak_time_ns' in columns:
    peak = map(LIGOTimeGPS, peak)
    attr_map['peak_time'], attr_map['peak_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in peak])

  if 'ms_start_time' in columns or 'ms_start_time_ns' in columns:
    ms_start = map(LIGOTimeGPS, ms_start)
    attr_map['ms_start_time'], attr_map['ms_start_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in ms_start])
  if 'ms_stop_time' in columns or 'ms_stop_time_ns' in columns:
    ms_stop = map(LIGOTimeGPS, ms_stop)
    attr_map['ms_stop_time'], attr_map['ms_stop_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in ms_stop])

  if 'flow' in columns:           attr_map['flow']           = flow
  if 'fhigh' in columns:          attr_map['fhigh']          = fhigh
  if 'bandwidth' in columns:      attr_map['bandwidth']      = fhigh-flow
  if 'ms_flow' in columns:        attr_map['ms_flow']        = flow
  if 'ms_fhigh' in columns:       attr_map['ms_fhigh']       = fhigh
  if 'ms_bandwidth' in columns:   attr_map['ms_bandwidth']   = ms_fhigh-ms_flow

  if 'central_freq' in columns:   attr_map['central_freq']   = (flow+fhigh)/2
  if 'peak_frequency' in columns: attr_map['peak_frequency'] = (flow+fhigh)/2

  if 'snr' in columns:            attr_map['snr']            = cle**(1/2)
  if 'ms_snr' in columns:         attr_map['ms_snr']         = ms_cle**(1/2)

  if 'cluster_size' in columns or 'param_one_value' in columns:
    attr_map['param_one_name'] = ['cluster_size'] * numtrigs
    attr_map['param_one_value'] = cls
  if 'cluster_number' in columns or 'param_two_value' in columns:
    attr_map['param_two_name'] = ['cluster_number'] * numtrigs
    attr_map['param_two_value'] = nev

  cols   = attr_map.keys()
  append = out.append
  for i in range(numtrigs):
    t = lsctables.SnglBurst()
    for c in cols: setattr(t, c, attr_map[c][i])
    if ifo!=None:
      t.ifo = ifo
    if channel!=None:
      t.channel = channel
    if not check_time or (check_time and t.get_peak() in span):
      append(t)
  
  return out

def fromhacrfile(fname, start=None, end=None, ifo=None, channel=None,\
                 columns=None):

  """
    Load triggers from a HACR format text file into a SnglBurstTable object.
    Use start and end to restrict the returned triggers, and give ifo and
    channel to fill those columns in the table.

    If columns is given as a list, only those columns in the table will be
    filled. This is advisable to speed up future operations on this table.

    Arguments :

      fname : file or str
        file object or filename path to read with numpy.loadtext

    Keyword arguments :

      start : float
        minimum peak time for returned triggers
      end : float
        maximum peak time for returned triggers
      ifo : str
        name of IFO to fill in table
      channel : str
        name of channel to fill in table
      columns : iterable
        list of columnnames to populate in table
  """

  # set columns
  if columns==None: columns = lsctables.SnglBurst.__slots__
  if start or end:
    if not start:
      start = 0
    if not end:
      end   = numpy.inf
    span = segments.segment(start, end)
    if 'peak_time' not in columns: columns.append('peak_time')
    if 'peak_time_ns' not in columns: columns.append('peak_time_ns')
    check_time = True
  else:
    check_time = False

  # generate table
  out = SnglTriggerTable('hacr', columns=columns)

  # force filename not file object
  if hasattr(fname, 'readline'):
    fh = fname
  else:
    fh = open(fname, 'r')
      
  # load data from file
  dat = loadtxt(fh)

  # close file if we opened it
  if not hasattr(fname, 'readline'):
    fh.close()

  if numpy.shape(dat) == (0,):
    return out
  elif numpy.shape(dat) == (8,):
    peak_time, peak_time_offset, freq, bandwidth, duration, n_pix, snr,\
    totPower = map(lambda n: numpy.asarray([n]), dat)
  elif len(dat)==8:
    peak_time, peak_time_offset, freq, bandwidth, duration, n_pix, snr,\
    totPower = dat
  else:
    raise ValueError("Wrong number of columns in HACR format file. "\
                     "Cannot read.")

  numtrigs = len(peak_time)

  if columns==None: columns = lsctables.SnglBurst.__slots__   
  if start or end:
    if not start:
      start = 0
    if not end:
      end   = numpy.inf
    span = segments.segment(start, end)
    if 'peak_time' not in columns: columns.append('peak_time')
    if 'peak_time_ns' not in columns: columns.append('peak_time_ns')
    check_time = True
  else:
    check_time = False

  attr_map = dict()

  peak = peak_time+peak_time_offset
  if 'start_time' in columns or 'start_time_ns' in columns:
    start = map(LIGOTimeGPS, peak-duration/2)
    attr_map['start_time'], attr_map['start_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in start])
  if 'stop_time' in columns or 'stop_time_ns' in columns:
    stop = map(LIGOTimeGPS, peak+duration/2)
    attr_map['stop_time'], attr_map['stop_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in stop])
  if 'peak_time' in columns or 'peak_time_ns' in columns:
    peak = map(LIGOTimeGPS, peak)
    attr_map['peak_time'], attr_map['peak_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in peak])

  if 'ms_start_time' in columns or 'ms_start_time_ns' in columns:
    ms_start = map(LIGOTimeGPS, peak-duration/2)
    attr_map['ms_start_time'], attr_map['ms_start_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in ms_start])
  if 'ms_stop_time' in columns or 'ms_stop_time_ns' in columns:
    ms_stop = map(LIGOTimeGPS, peak+duration/2)
    attr_map['ms_stop_time'], attr_map['ms_stop_time_ns'] =\
        zip(*[(s.seconds, s.nanoseconds) for s in ms_stop])

  if 'duration' in columns:       attr_map['duration']       = duration
  if 'ms_duration' in columns:    attr_map['ms_duration']    = duration
  if 'central_freq' in columns:   attr_map['central_freq']   = freq
  if 'peak_frequency' in columns: attr_map['peak_frequency'] = freq
  if 'peak_frequency_error' in columns: attr_map['peak_frequency_error'] = bandwidth/freq

  if 'flow' in columns:           attr_map['flow']           = freq-bandwidth/2
  if 'fhigh' in columns:          attr_map['fhigh']          = freq+bandwidth/2
  if 'bandwidth' in columns:      attr_map['bandwidth']      = bandwidth
  if 'ms_flow' in columns:        attr_map['ms_flow']        = freq-bandwidth/2
  if 'ms_fhigh' in columns:       attr_map['ms_fhigh']       = freq+bandwidth/2
  if 'ms_bandwidth' in columns:   attr_map['ms_bandwidth']   = bandwidth

  if 'snr' in columns:            attr_map['snr']            = snr
  if 'ms_snr' in columns:         attr_map['ms_snr']         = snr

  if 'peak_time_offset' in columns or 'param_one_value' in columns\
  or 'peak_time_ns' in columns:
    attr_map['param_one_name'] = ['peak_time_offset'] * numtrigs
    attr_map['param_one_value'] = peak_time_offset
  if 'numPixels' in columns or 'param_two_value' in columns:
    attr_map['param_two_name'] = ['numPixels'] * numtrigs
    attr_map['param_two_value'] = n_pix
  if 'totPower' in columns or 'param_three_value' in columns:
    attr_map['param_three_name'] = ['totPower'] * numtrigs
    attr_map['param_three_value'] = totPower

  cols   = attr_map.keys()
  append = out.append
  for i in range(numtrigs):
    t = lsctables.SnglBurst()
    for c in cols: setattr(t, c, attr_map[c][i])
    if ifo!=None:
      t.ifo = ifo
    if channel!=None:
      t.channel = channel
    if not check_time or (check_time and t.get_peak() in span):
      append(t)
  
  return out

def fromihopefile(fname, start=None, end=None, ifo=None, channel=None,\
                  columns=None):

  """
    Load triggers from an iHope format CSV file into a SnglInspiralTable object.
    Use start and end to restrict the returned triggers, and give ifo and
    channel to fill those columns in the table.

    If columns is given as a list, only those columns in the table will be
    filled. This is advisable to speed up future operations on this table.

    Arguments :

      fname : file or str
        file object or filename path to read with numpy.loadtext

    Keyword arguments :

      start : float
        minimum peak time for returned triggers
      end : float
        maximum peak time for returned triggers
      ifo : str
        name of IFO to fill in table
      channel : str
        name of channel to fill in table
      columns : iterable
        list of columnnames to populate in table
  """

  # get columns
  def_cols = dict(enumerate(['end_time','end_time_ns','ifo','snr','mass1',\
                             'mass2', 'mtotal','eta','event_duration',\
                             'template_duration','eff_distance','chisq',\
                             'chisq_dof','bank_chisq','bank_chisq_dof',\
                             'cont_chisq','cont_chisq_dof']))

  if columns==None: columns = lsctables.SnglInspiral.__slots__   
  if start or end:
    if not start:
      start = 0
    if not end:
      end   = numpy.inf
    span = segments.segment(start, end)
    if 'end_time' not in columns: columns.append('end_time')
    if 'end_time_ns' not in columns: columns.append('end_time_ns')
    check_time = True
  else:
    check_time = False

  usecols = [t for t in def_cols if def_cols[t] in columns]
  
  # force filename not file object
  if hasattr(fname, 'readline'):
    fh = fname
  else:
    fh = open(fname, 'r')

  # load data from file
  dat = loadtxt(fh, usecols)

  # close file if we opened it
  if not hasattr(fname, 'readlin'):
    fh.close()

  if usecols:
    numtrigs = len(dat[0])
  else:
    numtrigs = 0

  # generate table
  out = SnglTriggerTable('ihope', columns=columns)

  cols   = numpy.arange(len(dat))
  append = out.append
  for i in range(numtrigs):
    t = lsctables.SnglInspiral()
    for c in cols:
      setattr(t, def_cols[usecols[c]],\
              int(data[c][i]) if re.search('time', def_cols[usecols[c]])\
              else data[c][i])
    if ifo!=None:
      t.ifo = ifo
    if channel!=None:
      t.channel = channel
    if not check_time or (check_time and t.get_peak() in span):
      append(t)
  
  return out

# ==============================================================================
# Time shift trigger table
# ==============================================================================

def time_shift(lsctable, dt=1):

  """
    Time shift lsctable by time dt.
  """

  out = table.new_from_template(lsctable)
  get_time = def_get_time(lsctable.tableName)

  for t in lsctable:
    t2 = copy.deepcopy(t)
    if re.search('inspiral', lsctable.tableName, re.I):
      t2.set_end(t2.get_end()+dt)
    elif re.search('burst', lsctable.tableName, re.I):
      t2.set_peak(t2.get_peak()+dt)
    elif re.search('ringdown', lsctable.tableName, re.I):
      t2.set_start(t2.get_start()+dt)
    out.append(t2)

  return out

# =============================================================================
# Read file
# =============================================================================

def loadtxt(fh, usecols=None):

  """
    Stripped down version of numpy.loadtxt to work with empty files.
  """

  _comment = re.compile('[#%]')
  _delim   = re.compile('[\t\,\s]+')
  output = []
  nVals = 0
  for i,line in enumerate(fh):
    if _comment.match(line): continue
    vals = _delim.split(line.rstrip())
    if nVals>0 and len(vals) != nVals: 
      print "Warning, line %d of file %s was skipped, uncorrect column number" % (i, fh)
      continue
    nVals = len(vals)
    if usecols is not None:
      output.append(tuple(map(float, [vals[j] for j in usecols])))
    else:
      output.append(tuple(map(float, vals)))
  return numpy.squeeze(numpy.array(output, float)).T

