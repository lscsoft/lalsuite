"""
This module defines the some useful functions for GSI socket servers.

This file is part of the Grid LSC User Environment (GLUE)

GLUE is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Scott Koranda <skoranda@gravity.phys.uwm.edu>"
from glue import git_version
__date__ = git_version.date
__version__ = git_version.id

import os
import sys
import re
from pyGlobus import security

class Gridmap(object):
  """
  Class to represent the grid-mapfile used for the server. 
  Methods __getitem__ and __setitem__ are added for convenience
  so that an instance can be treated like a dictionary.

  Entries in the grid-mapfile should be on a single line and
  surrounded by double quotes. No mapping to a local user ID 
  is necessary.
  """
  def __init__(self, filePath=None, logger=None):
    """
    Initialize a new instance. If filePath is None then check
    the environment for the variable GRIDMAP and uses the 
    path found as the path to the grid-mapfile to open and parse.
    
    @return: an instance that can be treated like a dictionary
    """
    if not filePath:
      try: 
        path = os.environ["GRIDMAP"]
      except:
        if logger:
          logger.error("Error parsing GRIDMAP from environment: %s" % e)
          logger.error("Authentication to server will fail for all subjects")
    else:   
      path = filePath
      self.path = path
      self.logger = logger


  def __getitem__(self, key):
    f = None
    try:
      f = open(self.path, 'r')
    except Exception as e:
      if self.logger:
        self.logger.error("Error opening grid-mapfile %s: %s" % (self.path, e))
        self.logger.error("Authentication to server will fail for all subjects")
                
    try:
      subjects = []
      lines = [ s.strip() for s in f.readlines()]
      for s in lines:
        a = s.split('"')
        for b in a:
          if b == '': continue
          else:
            subjects.append(b)
            break
    except Exception as e:
      if self.logger:
        self.logger.error("Error parsing grid-mapfile %s: %s" % (self.path, e))
        self.logger.error("Authentication to server may fail for some subjects")

    self.d = {} 
    for s in subjects:
      self.d[s] = 1

    if f:
      f.close()

    if self.d.has_key(key): return 1
    else: return 0


class AuthCallback(object):
  def __init__(self, gridmap, logger, callback=None):
    self.gridmap = gridmap
    self.logger = logger
    self.callback = callback

  def __call__(self, arg, handle, identity, context):
    """
    Callback function for GSI authentication. Each time this 
    function is called a new instance of the GridMap class is created
    that reads the grid-mapfile. If the subject passed in as the identity
    is in the grid-mapfile then 1 is returned and the connection is
    authenticated and allowed to continue. If the subject passed is
    not in grid-mapfile then a 0 is returned and the connection is
    broken.

    This function is supplied to help create an instance of io.AuthData
    which is then supplied to help create an instance of io.TCPIOAttr,
    which is then supplied to help create the instance of 
    io.ThreadingGSITCPSocketServer. See the
    U{pyGlobus<http://www-itg.lbl.gov/gtg/projects/pyGlobus/>} documentation.

    Everything is within a try/except clause since if an exception is raised
    within this function experience shows that the connection will be
    allowed to authenticate!

    @param handle: handle to the io.AuthData instance, not used here
    @param identity: the identity of the party trying to authenticate
    @type identity: string
    @param context: the GSI context, not used here

    @return: 1 to accept and 0 to reject
    """
    try:
      # load gridmap file
      if self.gridmap[identity]:
        if self.logger:
          self.logger.info("Accepted connection from %s" % identity)
        if self.callback:
          self.callback(arg, 1)
        return 1
      else:
        if self.logger:
          self.logger.info("Rejected connection from %s" % identity)
        if self.callback:
          self.callback(arg, 0)
        return 0
    except Exception as e:
      try:
        self.logger.error("Error within authentication callback: %s" % e)
      except Exception as e:
        pass
      return 0


def daemon():
  """
  Common code for making a process into a daemon. 

  Fork once, change directory to /, get a new session, set the umask to 0,
  then fork again.

  @return: Parent never returns. The child gets None.
  """
  # do the first fork
  try:
    pid = os.fork()
    if pid > 0:
      # exit first parent
      sys.exit(0)
  except OSError as e:
    sys.stderr.write("fork #1 failed: %d (%s)\n" % (e.errno, e.strerror))
    sys.exit(1)
        
  # decouple from parent environment
  os.chdir("/")
  os.setsid()
  os.umask(0)
        
  # do the second fork
  try:
    pid = os.fork()
    if pid > 0:
      # exit from second parent
      sys.exit(0)
  except OSError as e:
    sys.stderr.write("fork #2 failed: %d (%s)\n" % (e.errno, e.strerror))
    sys.exit(1)
        
  # send stdout and stderr to /dev/null
  try:
    devnull = open("/dev/null", "w")
    sys.stdout = devnull
    sys.stderr = devnull
  except Exception as e:
    sys.__stderr__.write("Unable to direct to /dev/null: %s\n" % e)
    sys.exit(1)


def checkCredentials():
    """
    Check to make sure that the proper Grid Credentials (a proxy certificate)
    is available in order to authenticate to the remote LSCsegFindServer.
    """
    # verify that we have access to credentials
    try:
      proxyText = security.grid_proxy_info()
    except Exception as e:
      sys.stderr.write("Error verifying credentials: %s\n" % e)
      sys.stderr.write("Run 'grid-proxy-init' to generate a proxy certificate\n")
      return False

    pat = re.compile(r'timeleft : (\d{1,3}):(\d\d):(\d\d)')

    try:
      if isinstance(proxyText, str):
        m = pat.search(proxyText)
      elif isinstance(proxyText, tuple):
        m = pat.search(proxyText[0])
      else:
        raise RuntimeError("bad format for proxyText in checkCredentials")
      hours, minutes, seconds = map(int, m.groups())

    except Exception as e:
      sys.stderr.write("Error parsing proxy information: %s\n" % e)
      return False

    timeleft = seconds + 60 * minutes + 3600 * hours

    if timeleft < 300:
      sys.stderr.write("Less than 5 minutes left for proxy certificate.\n")
      sys.stderr.write("Run 'grid-proxy-init' to generate a new proxy certificate\n")
      return False

    return True
