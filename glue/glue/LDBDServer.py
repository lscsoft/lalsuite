"""
The LDBDServer module provides an API for responding to request from the
LDBDClient by connecting to the DB2 database.

This module requires U{pyGlobus<http://www-itg.lbl.gov/gtg/projects/pyGlobus/>}.


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

from glue import git_version
__date__ = git_version.date
__version__ = git_version.id

import os
import sys
import re
import types
try:
    import pyRXP
except ImportError:
    import pyRXPU as pyRXP
import exceptions
import socket
import SocketServer
import cPickle
from glue import ldbd
try:
  import rlsClient
except:
  pass


def dtd_uri_callback(uri):
  if uri == 'http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt':
    # if the XML file contants a http pointer to the ligolw DTD at CIT then
    # return a local copy to avoid any network problems
    return 'file://localhost' + os.path.join( os.environ["GLUE_PREFIX"],
      'etc/ligolw_dtd.txt' )
  else:
    # otherwise just use the uri in the file
    return uri

def initialize(configuration,log):
  # define the global variables used by the server
  global logger, max_bytes, xmlparser, dbobj, rls
  global dmt_proc_dict, dmt_seg_def_dict, creator_db
  global ldbd_com, db2_com, port

  # initialize the logger
  logger = log
  logger.info("Initializing server module %s" % __name__ )
  
  # initialize the database hash table
  dbobj = ldbd.LIGOMetadataDatabase(configuration['dbname'])
  max_bytes = configuration['max_client_byte_string']

  # initialize the ldbd commands and db2 command restrictions
  port = configuration['port']
  ldbd_com = configuration['ldbd_com'].split(',')
  db2_com = configuration['db2_com'].split(',')

  # create the xml parser
  xmlparser = pyRXP.Parser()

  # use a local copy of the DTD, if one is available
  try:
    GLUE_PREFIX = os.environ["GLUE_PREFIX"]
    xmlparser.eoCB = dtd_uri_callback
    logger.info("Using local DTD in " + 
      'file://localhost' + os.path.join( GLUE_PREFIX, 'etc') )
  except KeyError:
    logger.warning('GLUE_PREFIX not set, unable to use local DTD') 

  # open a connection to the rls server
  rls_server = configuration['rls']
  cert = configuration['certfile']
  key = configuration['keyfile']
  try:
    rls = rlsClient.RlsClient(rls_server,cert,key)
  except:
    rls = None

  # initialize dictionaries for the dmt processes and segments definers
  dmt_proc_dict = {}
  dmt_seg_def_dict = {}
  creator_db = None

def shutdown():
  global logger, max_bytes, xmlparser, dbobj, rls
  global dmt_proc_dict, dmt_seg_def_dict
  logger.info("Shutting down server module %s" % __name__ )
  if rls:
    del rls
  del xmlparser
  del dbobj
  del dmt_proc_dict
  del dmt_seg_def_dict

class ServerHandlerException(exceptions.Exception):
  """
  Class representing exceptions within the ServerHandler class.
  """
  def __init__(self, args=None):
    """
    Initialize an instance.

    @param args: 

    @return: Instance of class ServerHandlerException
    """
    self.args = args
        
class ServerHandler(SocketServer.BaseRequestHandler):
  """
  An instance of this class is created to service each request of the server.
  """
  def handle(self):
    """
    This method does all the work of servicing a request to the server. See
    the documentation for the standard module SocketServer.

    The input from the socket is parsed for the method with the remaining
    strings stripped of null bytes passed to the method in a list.

    There are no parameters. When the instance of the class is created to
    process the request all necessary information is made attributes of the
    class instance.

    @return: None
    """
    global logger
    global max_bytes

    logger.debug("handle method of %s class called" % __name__)

    # mapping of ldbdd RPC protocol names to methods of this class
    methodDict = {
      'PING' : self.ping,
      'QUERY' : self.query,
      'INSERT' : self.insert,
      'INSERTMAP' : self.insertmap,
      'INSERTDMT' : self.insertdmt
    }

    if True:
      # from the socket object create a file object
      self.sfile = self.request.makefile("rw")
      f = self.sfile

      # read all of the input up to limited number of bytes
      input = f.read(size=max_bytes,waitForBytes=2)

      # try 10 more times if we don't have a null byte at the end
      while input[-1] != '\0':
        input += f.read(size=max_bytes,waitForBytes=2)

      # the format should be a method string, followed by a null byte
      # followed by the arguments to the method encoded as null
      # terminated strings

      # check if the last byte is a null byte
      if input[-1] != '\0':
        logger.error("Bad input on socket: %s" % input)
        raise ServerHandlerException("Last byte of input is not null byte")
    #except Exception, e:
    #  logger.error("Error reading input on socket: %s" %  e)
    #  return

    try:
      # parse out the method and arguments 
      stringList = input.split('\0')
      methodString = stringList[0]
      argStringList = stringList[1:-1]
                        
    except Exception, e:
      logger.error("Error parsing method and argument string: %s" % e)

      msg = "ERROR LDBDServer Error: " + \
        "Error parsing method and argument string: %s" % e
      self.__reply__(1, msg)
      return

    # Set read and read/write access to different ports according to their configuration file
    try:
      # ldbd_com lists the allowed methodString for ldbd server based on the port number it is running on
      if methodString in ldbd_com:
        pass
      else:
        msg = '\nTo %s, authorized user please switch to the Apache driven ldbd server through ' % methodString
        msg += '\nsecure connection by specifying procotol "https" in your --segment-url argument'
        msg += '\nFor example, "--segment-url https://segdb.ligo.caltech.edu"'
        logger.error(msg) 
        self.__reply__(1,msg)
        raise ServerHandlerException(msg)
      # list allowed sql commands when methodString is QUERY
      if methodString=='QUERY':
        # get the sql command, for example "SELECT, UPDATE, INSERT"
        # see if the command is in the .ini file "db2_com" variable 
        dbcommand =  argStringList[0].split(' ')[0].upper()
        if dbcommand in db2_com:
          pass
        else:
           msg = 'ldbd server on port %d DO NOT support "%s"' % (port, dbcommand)
           logger.error(msg) 
           self.__reply__(1,msg)
           raise ServerHandlerException(msg)
    except Exception, e:
      logger.error("Error filtering allowed commands: %s" % e) 
      return

                
    try:
      # look up method in dictionary
      method = methodDict[methodString]
    except Exception, e:
      msg = "Error converting method string %s to method call: %s" % \
        (methodString, e)
      logger.error(msg)
                        
      self.__reply__(1, msg)
      return

    try:
      # call the method requested with the rest of strings as input
      result = method(argStringList) 
      self.__reply__( result[0], result[1] )
    except Exception, e:
      logger.error("Error while calling method %s: %s" % (methodString, e))

    return
        
  def __reply__(self, code, msg):
    """
    Format and send a reply back down the socket to the client. The file
    representing the socket is closed at the end of this method.

    @param code: integer representing the error code with 0 for success
                
    @param msg: object to be passed back to the client, either a string
    or a list of items that can be represented by strings
                        
    @return: None
    """
    f = self.sfile
    reply = "%d\0%s\0" % (code, msg)
    f.write(reply)

    # close the file associated with the socket
    f.close()

  def ping(self, arg):
    """
    Bounce back alive statment. Corresponds to the PING method in the
    ldbdd RPC protocol.

    @param arg: list (perhaps empty) of strings representing message sent
      by client to server

    @return: None
    """
    global logger

    logger.debug("Method ping called")
    try:
      hostname = socket.getfqdn()
      msg = "%s at %s is alive" % (__name__, hostname)
    except Exception, e:
      msg = "%s is alive" % __name__

    return (0, msg)

  def query(self, arg):
    """
    Execute an SQL query on the database and return the result as LIGO_LW XML

    @param arg: a text string containing an SQL query to be executed

    @return: None
    """
    global logger
    global xmlparser, dbobj

    # get the query string and log it
    querystr = arg[0]
    logger.debug("Method query called with %s" % querystr)

    # assume failure
    code = 1

    try:
      # create a ligo metadata object
      lwtparser = ldbd.LIGOLwParser()
      ligomd = ldbd.LIGOMetadata(xmlparser,lwtparser,dbobj)

      # execute the query
      rowcount = ligomd.select(querystr)

      # convert the result to xml
      result = ligomd.xml()

      logger.debug("Method query: %d rows returned" % rowcount)
      code = 0
    except Exception, e:
      result = ("Error querying metadata database: %s" % e)
      logger.error(result)

    try:
      del ligomd
      del lwtparser
    except Exception, e:
      logger.error(
        "Error deleting metadata object in method query: %s" % e)

    return (code,result)

  def insert(self, arg):
    msg = '\nTo INSERT, authorized user please switch to the Apache driven ldbd server through '
    msg += '\nsecure connection by specifying procotol "https" in your --segment-url argument'
    msg += '\nFor example, "--segment-url https://segdb.ligo.caltech.edu"'
    logger.error(msg)
    return (1, msg)


  def insertmap(self, arg):
    msg = '\nTo INSERTMAP, authorized user please switch to the Apache driven ldbd server through '
    msg += '\nsecure connection by specifying procotol "https" in your --segment-url argument'
    msg += '\nFor example, "--segment-url https://segdb.ligo.caltech.edu"'
    logger.error(msg)
    return (1, msg)

  def insertdmt(self, arg):
    msg = '\nTo INSERTDMT, authorized user please switch to the Apache driven ldbd server through '
    msg += '\nsecure connection by specifying procotol "https" in your --segment-url argument'
    msg += '\nFor example, "--segment-url https://segdb.ligo.caltech.edu"'
    logger.error(msg)
    return (1, msg)


