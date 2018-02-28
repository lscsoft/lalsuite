#!/usr/bin/env python

# Copyright (C) 2011 Duncan Macleod
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import re,os,sys,numpy,subprocess,datetime,shlex,urlparse,glob,fnmatch,\
       httplib,cjson,copy
from socket import getfqdn

from glue import segments,git_version
from glue.lal import Cache as LALCache
from glue.lal import CacheEntry as LALCacheEntry

from pylal import Fr,date
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.dq.dqDataUtils import make_external_call

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

"""
This module provides frame manipulation routines for use in data quality investigations, and cache manipulation routines.
"""

# ==============================================================================
# Grab data from frames
# ==============================================================================

def fromframefile(filename, channel, start=None, end=None):

  """
    Extract 1D data for given channel from the GWF format frame filename.
    Returns (x, data) pair of numpy arrays (x is array of times or frequencies.

    Arguments:

      filename : string
        path to GWF frame file
      channel : string
        channel name to extract

    Keyword arguments:

      start : float
        GPS start time (s) or minimum frequency (Hz) to return
      end : float
        GPS end time (s) or maximum frequency (Hz) to return
  """

  # try to extract data from frame
  y, fstart, offset, dt = Fr.frgetvect1d(filename, str(channel))[:4]
  x = fstart+dt*numpy.arange(len(y))+offset

  # apply constraint on x-axis
  if start or end:
    if not start: start=-numpy.infty
    if not end:   end=numpy.infty
    condition = (x>=start) & (x<end)
    y = y[condition]
    x = x[condition]

  return x,y

def toframefile(filename, channel, data, start, dx, **frargs):

  """
    Write numpy array data to GWF frame file using the given arguments.

    Arguments:

      filename : string
        name of file to write
      channel : string
        name of channel to write
      data : numpy.array
        array of data to write
      start : float
        GPS start time (s) or minimum frequency (Hz)
      dx : float
        GPS time step (s) or frequency step (Hz)

    Unnamed arguments are held in frargs. For usage, see documentation for
    pylal.Fr.frputvect.
  """

  datadict = frargs
  datadict['name']  = channel
  datadict['data']  = data
  datadict['start'] = start
  datadict['dx']    = dx

  Fr.frputvect(filename, [datadict], verbose=False)

def fromLALCache(cache, channel, start=None, end=None, verbose=False):

  """
    Extract data for given channel from glue.lal.Cache object cache. Returns 
    (time, data) pair of numpy arrays.

    Arguments:

      cache : glue.lal.Cache
        Cache list of frame files to read
      channel : string
        channel name to extract

    Keyword arguments:

      start : float
        GPS start time (s) or minimum frequency (Hz) to return
      end : float
        GPS end time (s) or maximum frequency (Hz) to return
  """

  # initialise data
  time = numpy.ndarray((1,0))
  data = numpy.ndarray((1,0))

  # set up counter
  if verbose:
    sys.stdout.write("Extracting data from %d frames...     " % len(cache))
    sys.stdout.flush()
    delete = '\b\b\b'
    num = len(cache)/100

  # loop over frames in cache
  for i,frame in enumerate(cache):
    # check for read access
    if os.access(frame.path, os.R_OK):
      # get data
      frtime, frdata = fromframefile(frame.path, channel, start=start,\
                                     end=end)
      # resize array and extend
      op = len(time[0])
      np = len(frtime)
      time.resize((1,op+np))
      time[0][op:] = frtime
      data.resize((1,op+np))
      data[0][op:] = frdata

      # print verbose message
      if verbose and len(cache)>1:
        progress = int((i+1)/num)
        sys.stdout.write('%s%.2d%%' % (delete, progress))
        sys.stdout.flush()

    else:
      raise RuntimeError("Cannot read frame\n%s" % frame.path)

  if verbose: sys.stdout.write("\n")

  return time[0], data[0]

def grab_data(start, end, channel, type, nds=False, dmt=False, verbose=False):

  """
    This function will return the frame data for the given channel of the given
    type in the given [start,end] time range and will construct a gps time
    vector to go with it. The nds option is not yet supported, and the dmt
    option will return data for dmt channels in frames not available on the
    datafind server.

    >>> grab_data(997315215, 997315225, 'G1:DER_DATA_QUALITY', 'G1_RDS_C01_L3')
    (array([  9.97315215e+08,   9.97315216e+08,   9.97315217e+08,
             9.97315218e+08,   9.97315219e+08,   9.97315220e+08,
             9.97315221e+08,   9.97315222e+08,   9.97315223e+08,
             9.97315224e+08]),
     array([ 256.,  256.,  256.,  256.,  256.,  256.,  256.,  256.,  256.,
             256.]))

    Arguments:

      start : float
        GPS start time (s).
      end : float
        GPS end time (s).
      channel : string
        channel name to extract, e.g. 'G1:DER_DATA_H'.
      type : string
        frame data type to use, e.g. 'G1_RDS_C01_L3'.

    Keyword arguments:

      nds : [ True | False ]
        use NDS connection to data server (UNSUPPORTED).
      dmt : [ True | False ]
        frame type is DMT product (DMT data not found by datafind server).
  """

  time = []
  data = []

  # generate framecache
  ifo = channel[0]
  if not dmt:
    cache = get_cache(start, end, ifo, type, verbose=verbose)
  else:
    cache = dmt_cache(start, end, ifo, type)

  time, data = fromLALCache(cache, channel, start=start, end=end,\
                            verbose=verbose)

  return time,data

# ==============================================================================
# Generate data cache
# ==============================================================================

def get_cache(start, end, ifo, ftype, framecache=False, server=None,\
              verbose=False):

  """
    Queries the LSC datafind server and returns a glue.lal.Cache object
    containing the frame file paths in the given GPS (start, end) interval
    for the given ifo and type (can be lists).

    framecache=True returns a pylal.dq.dqFrameUTils.FrameCache object in stead.

    Arguments:

      start : float
        GPS start time (s).
      end : float
        GPS end time (s).
      ifo : [ string | list ]
        ifo (or list of) to find, e.g. 'G1'.
      ftype : [ string | list ]
        frame data type (or list of) to find, e.g. 'G1_RDS_C01_L3'.

  """

  # set lists
  if isinstance(ftype, str):
    types = [ftype]
  else:
    types = ftype
  if isinstance(ifo, str):
    ifos = [ifo]
  else:
    ifos = ifo

  # construct span
  span = segments.segment(start,end)

  # set options
  cache = LALCache()
  entry_class = LALCacheEntry

  # try querying the ligo_data_find server
  if not server:
    server = _find_datafind_server()

  if verbose: sys.stdout.write("Opening connection to %s...\n" % server)

  if re.search(':', server):
    port = int(server.split(':')[-1])
  else:
    port = None

  cert, key = _get_grid_proxy()

  # if we have a credential then use it when setting up the connection
  if cert and key and port!=80:
    h = httplib.HTTPSConnection(server, key_file=key, cert_file=cert)
  else:
    h = httplib.HTTPConnection(server)

  if verbose: sys.stdout.write("Querying server for frames...\n")

  # loop over ifos and types
  for ifo in ifos:
    for t in types:
      # construct the URL for a simple data find query
      url = "/LDR/services/data/v1/gwf/%s/%s/%s,%s/file.json" % (ifo[0], t,\
                                                                 str(start),\
                                                                 str(end))
      # query the server
      h.request("GET", url)
      response = h.getresponse()
      _verify_response(response)
      # unravel the response
      urlList = cjson.decode(response.read())
      for url in urlList:
        cache.append(entry_class.from_T050017(url))

  # close the server connection
  h.close()
  if verbose: sys.stdout.write("Connection to %s closed.\n" % server)

  # convert to FrameCache if needed
  if framecache:
    cache = LALCachetoFrameCache(cache)

  return cache

# =============================================================================
# Return latest frame for given type
# =============================================================================

def get_latest_frame(ifo, ftype):

  """
    Returns the latest frame available in the LSC datafind server for the given
    ifo and frame type.

    Arguments:

      ifo : string
        observatory to find
      ftype : string
        frame data type to find
  """

  url = '/LDR/services/data/v1/gwf/%s/%s/latest/file.json' % (ifo[0], ftype)
  frame = query_datafind_server(url)

  if isinstance(frame, list) and len(frame)==1:
    return frame[0]
  else:
    return None

# =============================================================================
# Find ifos
# =============================================================================

def find_ifos():

  """
    Query the LSC datafind server and return a list of sites for which data
    is available. Does not differentiate between H1 and H2.

    Example:
   
    >>> find_ifos()
    ['G', 'H', 'L', 'V']
  """
  query_url = "/LDR/services/data/v1/gwf.json"
  reply = query_datafind_server(query_url)
  if reply:
    return [i for i in reply if len(i)==1]
  else:
    return []

# =============================================================================
# Find types
# =============================================================================

def find_types(ifo=[], ftype=[], search='standard'):

  """
    This function will return a valid list of LIGO frame types given the list of
    type strings. The search option defines the breadth of the search, to speed
    up the search, the following search options are supported:
    'standard','short','full'. 

    The 'R', 'T', and 'M' (raw, raw second trends, and raw minute trends) are 
    treated as special cases, so as not to return all types containing those 
    letters. 

    Example:

    >>>find_types(ftype='H1_RDS')
    ['H1_RDS_C01_LX',
     'H1_RDS_C02_LX',
     'H1_RDS_C03_L1',
     'H1_RDS_C03_L2',
     'H1_RDS_C03_L2_ET',
     'H1_RDS_C03_L2_ET2',
     'H1_RDS_C03_L2_ET30',
     'H1_RDS_C04_LX',
     'H1_RDS_R_L1',
     'H1_RDS_R_L3',
     'H1_RDS_R_L4']

    >>>find_types(ftype=['H1_RDS','R'],search='short')
    ['H1_RDS_R_L1', 'H1_RDS_R_L3', 'H1_RDS_R_L4', 'R']

    Keyword arguments:

      ifo : [ string | list ]
        ifo (or list of) to find, e.g. 'G1'.
      ftype : [ string | list ]
        frame data type (or list of) to find, e.g. 'G1_RDS_C01_L3'.
      search : string
        descriptor of how deep to search for frame type.
  """

  # make sure types is a list
  if isinstance(ftype, str):
    types = [ftype]
  else:
    types = ftype
  if isinstance(ifo, str):
    ifos = [ifo]
  else:
    ifos = ifo

  if not types:
    types = None

  # treat 'R','M' and 'T' as special cases,
  special_types = ['M','R','T']
  foundtypes = []

  # there are thousands of GRBXXXXXX frame types, so ignore them
  ignore = []
  if search!='full': 
    ignore.extend(['GRB'])
  if search=='short':
    # all of these strings are part of frame types that can be ignored for a
    # short search
    short_ignore = ['CAL','BRST','Mon','SG','IMR','DuoTone','Concat',\
                    'BH','WNB','Lock','_M','_S5','Multi','Noise','_C0']
    ignore.extend(short_ignore)
  ignore = '(%s)' % '|'.join(ignore)

  # set up server connection
  server    = _find_datafind_server()
  cert, key = _get_grid_proxy()
  if re.search(':', server):
    port = int(server.split(':')[-1])
  else:
    port = None

  if cert and key and port!=80:
    h = httplib.HTTPSConnection(server, key_file=key, cert_file=cert)
  else:
    h = httplib.HTTPConnection(server)

  # query for individual sites in datafind server
  if not ifos:
    query_url = "/LDR/services/data/v1/gwf.json"
    h.request("GET", query_url)
    response = h.getresponse()
    _verify_response(response)
    ifos = [i for i in cjson.decode(response.read()) if len(i)==1]

  # query for types for each ifo in turn
  datafind_types = []
  for ifo in ifos:
    # find types from datafind server
    query_url = "/LDR/services/data/v1/gwf/%s.json" % ifo[0]
    h.request("GET", query_url)
    response = h.getresponse()
    _verify_response(response)
    datafind_types.extend(cjson.decode(response.read()))
  
  # close connection
  h.close()

  # find special types first, otherwise they'll corrupt the general output
  r = 0
  for i,t in enumerate(datafind_types):
    if (types==None and t in special_types) or\
       (types!=None and t in types and t in special_types):
      foundtypes.append(t)
      datafind_types.pop(i-r)
      if types is not None:
        types.pop(types.index(t))
      r+=1;
      continue

  # find everything else
  for t in datafind_types:
    if re.search(ignore, t):
      continue
    # if no types have been specified, return all
    if types==None:
      foundtypes.append(t)
    # else check for a match
    else:
      if len(types)>=1 and re.search('(%s)' % '|'.join(types), t):
        foundtypes.append(t)

  foundtypes.sort()
  types = []
  for t in foundtypes:
    if t not in types:  types.append(t)

  return types

# =============================================================================
# Functions to find channels
# =============================================================================

def get_channels(framefile):

  """
    Extract the channels held in a frame file.
  """

  # get type
  type = os.path.basename(framefile).split('-')[1]

  # get channels contained in frame, grepping for input channel string
  frchannels,err = make_external_call('FrChannels %s' % framefile)

  channels = ChannelList()
  for line in frchannels.splitlines():
    name, samp = line.split()
    channels.append(Channel(name, sampling=samp, type=type))

  return channels

def find_channels(name=[], ftype=[], ifo=[], not_name=[], not_ftype=[],\
                  exact_match=False, time=None, unique=False):

  """
    Returns a ChannelList containing Channels for all data channels matching
    the given attributes from frames matching the given GPS time (defaults to
    now).

    Keyword arguments:

      name : [ string | list ]
        channel name (or list of) to match in search. Can be part of name,
        e.g. ['STRAIN', 'DARM_ERR']
      ftype : [ string | list ]
        frame data type (or list of) to find, e.g. 'G1_RDS_C01_L3'.
      ifo : [ string | list ]
        ifo (or list of) to find, e.g. 'G1'.
      not_name : [ string | list ]
        channel name (or list of) to negatively match in search. Can be part of
        name, e.g. 'ETMY_EXC_DAQ'
      not_ftype : [ string | list ]
        frame data type (or list of) to remove from search.
      exact_match : [ True | False]
        require complete match with given name list, not just partial match.
      time : float
        GPS time to which to restrict search. Data transfer latency means that
        the very latest data is not always available, best to give a 'recent'
        time.
      unique : [ True | False ]
        return unique list of channels from different types (since the same
        channel can exist in multiple types).
  """

  # check list status
  if isinstance(name, str):      names = [name]
  else:                          names = name
  if isinstance(ftype, str):     types = [ftype]
  else:                          types = ftype
  if isinstance(ifo, str):      ifos  = [ifo]
  else:                          ifos  = ifo
  if isinstance(not_ftype, str): not_types = [not_ftype]
  else:                          not_types = not_ftype
  if isinstance(not_name, str): not_names = [not_name]
  else:                          not_names = not_name

  # find ifos
  if not ifos:
    ifos = find_ifos()

  # find types
  if not types:
    types = find_types(ifo=ifos, ftype=types)

  # remove types we don't want
  if not_types:
    if exact_match:
      notmatch = re.compile("\A(%s)\Z" % '|'.join(not_types))
    else:
      notmatch = re.compile("(%s)" % '|'.join(not_types))
    types = [t for t in types if not re.search(notmatch, t)]

  channels = ChannelList()

  # loop over each ifo
  for ifo in ifos:
    for type in types:
      frchannels  = ChannelList()

      # find first frame file for type
      if time:
        frame = get_cache(time, time, ifo, type)
        if len(frame):
          frame = frame[-1].path
      else:
        frame = get_latest_frame(ifo, type)

      # if frames not found, move on
      if frame:
        frame = urlparse.urlparse(frame)[2]
      else:
        continue

      # access frame
      if os.access(frame, os.R_OK):
        # get channels contained in frame, sieveing when necessary
        allchannels = get_channels(frame)
        # test ifo
        allchannels = allchannels.sieve(ifo=ifo) 

        # test names
        if names:
          for ch in names:
            if re.match('%s\d:' % ifo[0], ch):
              ch = ch.split(':',1)[1]
            frchannels.extend(allchannels.sieve(name=ch,\
                                                exact_match=exact_match))
        else:
          frchannels = allchannels

        # test excluded names
        for ch in not_names:
          frchannels = frchannels.sieve(not_name=ch)

      channels.extend(frchannels)

  if unique:
    channels = channels.unique()

  channels.sort(key=lambda c: str(c))

  return channels

# ==============================================================================
# Class to generate channel structure
# ==============================================================================

class Channel:
  """
    The Channel class defines objects to represent LIGO data channels. Each
    Channel has a 'name' attribute and can be assigned 'type' and 'sampling'
    attributes if relevant.

    Example:

    >>>GWChannel = Channel('H1:LSC-DARM_ERR, 'H1_RDS_R_L3', 4096)
    >>>GWChannel.name, GWChannel.type, GWChannel.sampling
    ('H1:LSC-DARM_ERR', 'H1_RDS_R_L3', 4096)
  """

  def __init__(self,name,type=None,sampling=None):
    """Initialise the Channel object.
    """

    attributes = ['ifo','site','name',\
                  'system','subsystem','signal',\
                  'type',\
                  'sampling']

    __slots__ = attributes

    # extract name attributes
    if ':' in name:
      self.ifo,self.name         = re.split(':',name,maxsplit=1)
      self.site                  = self.ifo[0]
    else:
      self.name = name
      self.ifo  = ""
      self.site = ""

    tags = re.split('[-_]',self.name,maxsplit=3)

    self.system = tags[0]

    if len(tags)>1:
      self.subsystem = tags[1]
    else:
      self.subsystem = ""

    if len(tags)>2:
      self.signal = tags[2]
    else:
      self.signal = ""

    if type:
      self.type = str(type)
    else:
      self.type = ""

    if sampling:
      self.sampling = float(sampling)
    else:
      self.sampling = 0

  def __getattribute__(self,name):

    return self.__dict__[name]

  def __str__(self):

    return '%s:%s' % (self.ifo, self.name)

class ChannelList(list):

  """
    Wrapper for a list of Channel objects, with helper functions.
  """

  def find(self, item):
    """
    Return the smallest i such that i is the index of an element that wholly
    contains item.  Raises ValueError if no such element exists.
    """
    for i,chan in enumerate(self):
      if item.name == chan.name:
        return i
    raise ValueError(item)

  def sieve(self, ifo=None, name=None, type=None, sampling=None,\
            sampling_range=None, not_name=None, exact_match=False):

    """
    Return a ChannelList object with those Channels that match the given
    attributes If exact_match is True, then non-None ifo, name, and type
    patterns must match exactly.

    If sampling is given, it will always test an exact match.

    If sampling_range is given, it will test for Channels with sampling in
    [min,max) of range.

    Bash-style wildcards (*?) are allowed for ifo, name and type.
    """

    if not exact_match:
      if ifo is not None:      ifo      = "*%s*" % ifo
      if name is not None:     name     = "*%s*" % name
      if not_name is not None: not_name = "*%s*" % not_name

    c = self
    if ifo is not None:
      ifo_regexp = re.compile(fnmatch.translate(ifo))
      c = [entry for entry in c if ifo_regexp.match(entry.ifo) is not None]

    if name is not None:
      name_regexp = re.compile(fnmatch.translate(name))
      c = [entry for entry in c if name_regexp.match(entry.name) is not None]

    if not_name is not None:
      name_regexp = re.compile(fnmatch.translate(not_name))
      c = [entry for entry in c if name_regexp.match(entry.name) is None]


    if sampling is not None:
      c = [entry for entry in c if entry.sampling==sampling]

    if sampling_range is not None:
      c = [entry for entry in c if\
           sampling_range[0] <= entry.sampling < sampling_range[1]]

    return self.__class__(c)

  def __isub__(self, other):
    """
    Remove elements from self that are in other.
    """
    end = len(self) - 1
    for i, elem in enumerate(self[::-1]):
      if elem in other:
        del self[end - i]
    return self

  def __sub__(self, other):
    """
    Return a ChannelList containing the entries of self that are not in other.
    """
    return self.__class__([elem for elem in self if elem not in other])

  def __ior__(self, other):
    """
    Append entries from other onto self without introducing (new) duplicates.
    """
    self.extend(other - self)
    return self

  def __or__(self, other):
    """
    Return a ChannelList containing all entries of self and other.
    """
    return self.__class__(self[:]).__ior__(other)

  def __iand__(self, other):
    """
    Remove elements in self that are not in other.
    """
    end = len(self) - 1
    for i, elem in enumerate(self[::-1]):
      if elem not in other:
        del self[end - i]
    return self

  def __and__(self, other):
    """
    Return a ChannelList containing the entries of self that are also in other.
    """
    return self.__class__([elem for elem in self if elem in other])

  def unique(self):
    """
    Return a ChannelList which has every element of self, but without
    duplication of channel name. Preserve order. Does not hash, so a bit slow.
    """
    new = self.__class__([])
    for elem in self:
      if str(elem) not in [str(e) for e in new]:
        new.append(elem)
    return new

# ==============================================================================
# Function to generate a framecache of /dmt types
# ==============================================================================
def dmt_cache(start,end,ifo,type,framecache=False):
  """
  This function will return a list of frame files in the given start and stop 
  time interval for the give IFO using the given DMT frame type. This is
  required if ligo_data_find will not return the dmt frames.

  Example:

  >>>dmt_cache(960000000,960010000,'H1','LockLoss_H1')
  ['/archive/frames/dmt/LHO/LockLoss_H1/H-M-960/H-LockLoss_H1_M-960001200-3600.gwf',
   '/archive/frames/dmt/LHO/LockLoss_H1/H-M-960/H-LockLoss_H1_M-960004800-3600.gwf']
  """

  # find dmt frames path
  host = getfqdn()
  if re.search('ligo-',host):
    dmt_dir = '/dmt'
  elif re.search('ligo.',host):
    site = {'H':'LHO','L':'LLO','V':'V1'}
    dmt_dir = os.path.join('/archive','frames','dmt',site[ifo[0]])


  span    = segments.segment(start,end)

  # set options
  if framecache:
    cache = FrameCache()
    entry_class = FrameCacheEntry
    ldfopt = '--frame-cache'
  else:
    cache = LALCache()
    entry_class = LALCacheEntry
    ldfopt = '--lal-cache'

  basedir = os.path.join(dmt_dir,type)

  # frames are 3600 seconds long, so round
  tmp = int(str(start)[0:3]+'000000')
  cache_start = tmp+3600*int((start-tmp)/3600)
  tmp = int(str(end)[0:3]+'000000')
  cache_end = tmp+3600*int((end-tmp)/3600)

  # directories are listed with three time digits
  start_three = int(str(cache_start)[0:3])
  end_three = int(str(cache_end)[0:3])
  first_three = numpy.arange(start_three,end_three+1,1)

  #loop over directories
  for t in first_three:
    querydir = os.path.join(basedir,'%s-%s-%s' % (ifo[0:1],'M',str(t)),'*')
    filenames = glob.glob(querydir)
    for filename in filenames:
      try:
        e = entry_class.from_T050017(filename)
        if span.intersects(e.segment):  cache.append(e)
      except ValueError:
        sys.stderr.write("Could not convert %s to %s\n"\
                         % (filename,'.'.join([entry_class.__module__,\
                                               entry_class.__name__])))

  return cache

# ==============================================================================
# Class for wCacheEntry
# ==============================================================================

class FrameCacheEntry(LALCacheEntry):
  """
    An object representing one line in a frame cache file.

    Each line in a frame cache identifies multiple files, and the line consists
    of six columns of white-space delimited text.

    The first column, "observatory", generally stores the name of an 
    observatory site or one or more instruments (preferably delimited by ",",
    but often there is no delimiter between instrument names).

    The second column, "description", stores a short string tag that is
    usually all capitals with "_" separating components, in the style of the
    description part of the LIGO-Virgo frame filename format.

    The third and fourth columns store the start time and stop time in GPS
    seconds of the interval spanned by the file identified by the cache line.

    The fifth column stored the duration of each frame identified in the cache
    line. 

    The sixth (last) column stores the file's URL.

    The values for these columns are stored in the .observatory,
    .description, .segment and .url attributes, respectively.  The
    .segment attribute stores a glue.segments.segment object describing
    the interval spanned by the file.  Any of these attributes except
    the URL is allowed to be None.
  """

  # How to parse a line in a frame cache file. Six white-space
  # delimited columns.
  _regex = re.compile(r"\A\s*(?P<obs>\S+)\s+(?P<dsc>\S+)\s+(?P<strt>\S+)\s+(?P<end>\S+)\s+(?P<dur>\S+)\s+(?P<url>\S+)\s*\Z")
  _url_regex = re.compile(r"\A((.*/)*(?P<obs>[^/]+)-(?P<dsc>[^/]+)-(?P<strt>[^/]+)-(?P<dur>[^/\.]+)\.[^/]+)\Z")

  def __init__(self, *args, **kwargs):
    """
    Intialize a FrameCacheEntry object. The arguments can take two
    forms:  a single string argument, which is interpreted and
    parsed as a line from a frame cache file, or four arguments
    used to explicitly initialize the observatory, description,
    segment and URL in that order.  When parsing a single line
    of text from a frame cache, an optional key-word argument
    "coltype" can be provided to set the type the start, end and
    durations are parsed as.  The default is glue.lal.LIGOTimeGPS.

    """
    if len(args) == 1:
      # parse line of text as an entry in a cache file
      match = self._regex.search(args[0])
      coltype = kwargs.pop("coltype", LIGOTimeGPS)

      try:
        match = match.groupdict()
      except AttributeError:
        raise ValueError("could not convert %s to FrameCacheEntry"\
                         % repr(args[0]))
      self.observatory = match["obs"]
      self.description = match["dsc"]
      start            = match["strt"]
      end              = match["end"]
      self.duration    = coltype(match["dur"])

      if start == "-" and end == "-":
        # no segment information
        self.segment = None
      else:
        self.segment = segments.segment(coltype(start),coltype(end))
      self.url = match["url"]

      if kwargs:
        raise TypeError("unrecognized keyword arguments: %s"\
                        % ", ".join(kwargs))
    elif len(args) == 5:
      # parse arguments as observatory, description,
      # segment, duration, url
      if kwargs:
        raise TypeError("invalid arguments: %s" % ", ".join(kwargs))
      self.observatory, self.description, self.segment, self.duration, self.url\
          = args
    else:
      raise TypeError("invalid arguments: %s" % args)

    # "-" indicates an empty column
    if self.observatory == "-":
      self.observatory = None
    if self.description == "-":
      self.description = None

  def __str__(self):
    """
    Convert the FrameCacheEntry to a string in the format of a line
    in a frame cache. Used to write the FrameCacheEntry to a file.

    """
    if self.segment is not None:
      start,end = [str(t) for t in self.segment]
    else:
      start    = "-"
      end      = "-"
      duration = "-"

    return "%s %s %s %s %s %s" % (self.observatory or "-", self.description or "-", start, end, self.duration, self.url)

  def __cmp__(self, other):
    """
    Compare two FrameCacheEntry objects by observatory, then
    description, then segment, then duration, then URL.
    """
    if type(other) != FrameCacheEntry:
      raise TypeError("can only compare FrameCacheEntry to FrameCacheEntry)")
    return cmp((self.observatory, self.description, self.segment,\
                self.duration, self.url), \
               (other.observatory, other.description, other.segment,\
                other.duration, other.url))

  def get_files(self):
    """
    Return Find all files described by this FrameCacheEntry.
    """

    filenames = glob.glob(os.path.join(self.path,\
                                         '%s-%s*-%s.*' % (self.observatory,\
                                                           self.description,\
                                                           self.duration)))
    cache = [e.path for e in\
                 LALCache([LALCacheEntry.from_T050017(f) for f in filenames])\
             if e.observatory==self.observatory and\
                e.description==self.description and\
                self.segment.intersects(e.segment) and\
                abs(e.segment)==self.duration]

    return cache
    
  def from_T050017(cls, url, coltype = LIGOTimeGPS):      

    """      
    Parse a URL in the style of T050017-00 into a FrameCacheEntry.      
    The T050017-00 file name format is, essentially,  
  
    observatory-description-start-dur.ext  
  
    """      
    match = cls._url_regex.search(url)      
    if not match:      
            raise ValueError("could not convert %s to CacheEntry" % repr(url))      
    observatory = match.group("obs")      
    description = match.group("dsc")      
    start = match.group("strt")      
    duration = match.group("dur")      
    if start == "-" and duration == "-":      
            # no segment information      
            segment = None      
    else:      
            segment = segments.segment(coltype(start), coltype(start) + coltype(duration))      
    return cls(observatory, description, segment, duration,\
               os.path.split(url)[0])    


  from_T050017 = classmethod(from_T050017)


# ==============================================================================
# Class for FrameCacheEntry 
# ==============================================================================

class FrameCache(LALCache):
  """
    An object representing a frame cache file. Currently it is possible to
    add anything to a FrameCache. This method should check that the thing you
    are adding is a FrameCacheEntry and throw and error if it is not.
  """
  entry_class = FrameCacheEntry
  
  def sieve(self, ifos=None, description=None, segment=None, duration=None,
    exact_match=False):
    """
    Return a FrameCache object with those FrameCacheEntries that contain the
    given patterns (or overlap, in the case of segment).  If
    exact_match is True, then non-None ifos, description, and
    segment patterns must match exactly.
    
    Bash-style wildcards (*?) are allowed for ifos and description.
    """
    if exact_match:
      segment_func = lambda e: e.segment == segment
    else:
      if ifos is not None: ifos = "*" + ifos + "*"
      if description is not None: description = "*" + description + "*"
      segment_func = lambda e: segment.intersects(e.segment)
    
    c = self
    
    if ifos is not None:
      ifos_regexp = re.compile(fnmatch.translate(ifos))
      c = [entry for entry in c if ifos_regexp.match(entry.observatory) is not None]
    
    if description is not None:
      descr_regexp = re.compile(fnmatch.translate(description))
      c = [entry for entry in c if descr_regexp.match(entry.description) is not None]
    
    if segment is not None:
      c = [entry for entry in c if segment_func(entry)]
    
    if duration is not None:
      c = [entry for entry in c if entry.duration==duration]

    return self.__class__(c)

  def get_files(self):
    c = []
    for e in self:
      c.extend(e.get_files())
    return c

def FrameCachetoLALCache(fcache):

  lcache = LALCache()

  files = fcache.get_files()

  for f in files:
    lcache.append(LALCacheEntry.from_T050017(f))
  
  return lcache

def LALCachetoFrameCache(lcache):

  lcache.sort(key=lambda e: (e.path,e.segment[0]))
  fcache = FrameCache()

  for e in lcache:
    matched = False
    dir = os.path.split(e.path)[0]

    # if path found in FrameCache try to coalesce with other entries
    dirs = [d.path for d in fcache]
    if dir in dirs:
      pathentries = [fe for fe in fcache if fe.path==dir]
      # test current entry against other entries for the same path in the
      # new cache
      for i,pe in enumerate(pathentries):
        notdisjoint = e.segment[0] <= pe.segment[1]
        # if entries match in directory, duration, and are contiguous, append
        if pe.path==os.path.split(e.path)[0]\
             and e.segment.__abs__() == pe.duration\
             and notdisjoint:
          seg = segments.segment(min(pe.segment[0], e.segment[0]),\
                                  max(pe.segment[1], e.segment[1]))
          fcache[fcache.index(pe)].segment = seg
          matched = True
          break

      # if we haven't matched the entry to anything already in the cache add now
      if not matched:
        fe = FrameCacheEntry.from_T050017(e.path)
        fcache.append(fe)

    # if from a new directory add
    else:
      fe = FrameCacheEntry.from_T050017(e.path)
      fcache.append(fe)

  return fcache

# =============================================================================
# Query the LSC Datafind Server
# =============================================================================

def query_datafind_server(url, server=None):

  # try querying the ligo_data_find server
  if not server:
    server = _find_datafind_server()
  if re.search(':', server):
    port = int(server.split(':')[-1])
  else:
    port = None

  cert, key = _get_grid_proxy()
  # if we have a credential then use it when setting up the connection
  if cert and key and port!=80:
    h = httplib.HTTPSConnection(server, key_file=key, cert_file=cert)
  else:
    h = httplib.HTTPConnection(server)

  # query the server
  h.request("GET", url)
  response = h.getresponse()
  _verify_response(response)

  # since status is 200 OK read the types
  body = response.read()
  h.close()
  if body == "":
    return None
  if url.endswith('json'):
    return cjson.decode(body)
  else:
    return body.splitlines()

def _get_grid_proxy():

  """
    Returns paths for X509 certificate and proxy if available.
  """

  # try and get a proxy or certificate
  # FIXME this doesn't check that it is valid, though
  cert = None
  key = None
  try:
    proxy = os.environ['X509_USER_PROXY']
    cert = proxy
    key = proxy
  except:
    try:
      cert = os.environ['X509_USER_CERT']
      key = os.environ['X509_USER_KEY']
    except:
      uid = os.getuid()
      proxy_path = "/tmp/x509up_u%d" % uid
      if os.access(path, os.R_OK):
        cert = proxy_path
        key = proxy_path

  return cert, key

def _verify_response(HTTPresponse):

  """
    Test response of the server to the query and raise the relevant exception
    if necessary.
  """

  if HTTPresponse.status != 200:
    msg = "Server returned code %d: %s" % (HTTPresponse.status,\
                                           HTTPresponse.reason)
    body = HTTPresponse.read()
    msg += body
    raise RuntimeError(msg)

def _find_datafind_server():

  """
    Find the LSC datafind server from the LIGO_DATAFIND_SERVER environment
    variable and raise exception if not found
  """

  var = 'LIGO_DATAFIND_SERVER'
  try:
    server = os.environ[var]
  except KeyError:
    raise RuntimeError("Environment variable %s is not set" % var)

  return server

# =============================================================================
# Read GEO control channels from frames
# =============================================================================

def parse_composite_channels(superchannel, ifo='G1', frame_type='R'):

  """
    Seperate GEO superchannels into seperate channel names.
  """

  # get system name
  tokens = re.split('#', superchannel)
  system = tokens.pop(0)

  channels = ChannelList()

  # parse channels
  while len(tokens) > 1:
    # get subsystem
    subsystem = tokens.pop(0)
    
    # get number of signals
    N = int(tokens.pop(0))

    for i in range(N):
      signal = tokens.pop(0)
      c = '%s:%s-%s_%s' % (ifo, system, subsystem, signal)
      channels.append(Channel(c, frame_type, 1))

  return channels

def separate_composite_data(data, N=1):

  """
    Seperate GEO superchannels data into array of data for each of the composite
    channels.
  """

  out = []
  for i in range(N):
    d = data[i::N]
    out.append(d)

  return out

def get_control_channel(cache, superchannel, channel, start=None, end=None,\
                        ifo='G1', frame_type='R'):

  """
    Get GEO control channel data for a single channel from a given superchannel,
    control channels are sampled at 1Hz.
  """

  # initialise data
  time = numpy.array([])
  data = numpy.array([])

  for frame in cache:
    # check for read access
    if os.access(frame.path, os.R_OK):
      frtime, frdata = fromframefile(frame.path, str(superchannel),\
                                     start=start, end=end)
      channels_list = parse_composite_channels(superchannel)
      chanIndex = channels_list.find(channel)
      channelData = separate_composite_data(frdata,\
                                            len(channels_list))[chanIndex]
      time = numpy.append(time, frtime)
      data = numpy.append(data, channelData)
    else:
      raise RuntimeError("Cannot read frame\n%s" % frame.path)

  return time, data
