"""
The LDBDClient module provides an API for connecting to and making requests of
a LDBDServer.

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

import sys
import os
import exceptions
import types
import re
import cPickle
import xml.parsers.expat

from pyGlobus import io
from pyGlobus import security


def version():
  return __version__


class SimpleLWXMLParser:
  """
  A very simple LIGO_LW XML parser class that reads the only keeps
  tables that do not contain the strings sngl_ or multi_

  The class is not very robust as can have problems if the line
  breaks do not appear in the standard places in the XML file.
  """
  def __init__(self):
    """
    Constructs an instance.

    The private variable ignore_pat determines what tables we ignore.
    """
    self.__p = xml.parsers.expat.ParserCreate()
    self.__in_table = 0
    self.__silent = 0
    self.__ignore_pat = re.compile(r'.*(sngl_|multi_).*', re.IGNORECASE)
    self.__p.StartElementHandler = self.start_element
    self.__p.EndElementHandler = self.end_element

  def __del__(self):
    """
    Destroys an instance by shutting down and deleting the parser.
    """
    self.__p("",1)
    del self.__p

  def start_element(self, name, attrs):
    """
    Callback for start of an XML element. Checks to see if we are
    about to start a table that matches the ignore pattern.

    @param name: the name of the tag being opened
    @type name: string

    @param attrs: a dictionary of the attributes for the tag being opened
    @type attrs: dictionary
    """
    if name.lower() == "table":
      for attr in attrs.keys():
        if attr.lower() == "name":
          if self.__ignore_pat.search(attrs[attr]):
            self.__in_table = 1
        
  def end_element(self, name):
    """
    Callback for the end of an XML element. If the ignore flag is
    set, reset it so we start outputing the table again.

    @param name: the name of the tag being closed
    @type name: string
    """
    if name.lower() == "table":
      if self.__in_table:
        self.__in_table = 0

  def parse_line(self, line):
    """
    For each line we are passed, call the XML parser. Returns the
    line if we are outside one of the ignored tables, otherwise
    returns the empty string.

    @param line: the line of the LIGO_LW XML file to be parsed
    @type line: string

    @return: the line of XML passed in or the null string
    @rtype: string
    """
    self.__p.Parse(line)
    if self.__in_table:
      self.__silent = 1
    if not self.__silent:
      ret = line
    else:
      ret = ""
    if not self.__in_table:
      self.__silent = 0
    return ret


class LDBDClientException(Exception):
  """Exceptions returned by server"""
  def __init__(self,args=None):
    self.args = args


class LDBDClient(object):
  def __init__(self, host, port, identity):
    """
    Open a connection to a LDBD Server and return an instance of
    class LDBDClient. One of the public methods can then be 
    called to send a request to the server.

    @param host: the host on which the LDBD Server runs
    @type host: string

    @param port: port on which the LDBD Server listens
    @type port: integer

    @param identity: string which the LDBD Server identifies itself
    @type identity: string

    @return: Instance of LDBDClient
    """
    try:
      self.__connect__(host,port,identity)
    except Exception, e:
      raise 

  def __del__(self):
    """
    Disconnect from the LDBD server.

    @return: None
    """
    self.__disconnect__()

  def __connect__(self,host,port,identity):
    """
    Attempt to open a connection to the LDBD Server
    using the 'host' and 'port' and expecting the server
    to identify itself with a corresponding host certificate.

    A IOException is raised if the connection cannot be made,
    but this is caught by the __init__ method above and 
    turned into a LDBDClient exception.
        
    @param host: the host on which the LDBD Server runs
    @type host: string

    @param port: port on which the LDBD Server listens
    @type port: integer

    @param identity: string which the LDBD Server identifies itself
    @type identity: string

    @return: None
    """
    # remove the globus tcp port range environment variable if set
    try:
      port_range = os.environ["GLOBUS_TCP_PORT_RANGE"]
      os.environ["GLOBUS_TCP_PORT_RANGE"] = ""
    except:
      pass
    
    self.host = host
    self.port = port
    self.identity = identity

    # redirect stdout and stderror for now
    try:
      f = open("/dev/null", "w")
      sys.stdout = f
      sys.stderr = f
    except:
      pass

    try:
      # create TCPIOAttr instance
      clientAttr = io.TCPIOAttr()
      authData = io.AuthData()
      soc = io.GSITCPSocket()

      if identity is None:
        # try an unauthenticated connection
        clientAttr.set_authentication_mode(
          io.ioc.GLOBUS_IO_SECURE_AUTHENTICATION_MODE_NONE)
        clientAttr.set_authorization_mode(
          io.ioc.GLOBUS_IO_SECURE_AUTHORIZATION_MODE_NONE, authData)
        clientAttr.set_channel_mode(
          io.ioc.GLOBUS_IO_SECURE_CHANNEL_MODE_CLEAR)
        clientAttr.set_delegation_mode(
          io.ioc.GLOBUS_IO_SECURE_DELEGATION_MODE_NONE)

      else:
        # set authentication mode to be GSSAPI
        clientAttr.set_authentication_mode(
          io.ioc.GLOBUS_IO_SECURE_AUTHENTICATION_MODE_GSSAPI)

        # set expected identity
        authData.set_identity(identity)

        # set authorization, channel, and delegation modes
        clientAttr.set_authorization_mode(
          io.ioc.GLOBUS_IO_SECURE_AUTHORIZATION_MODE_IDENTITY, authData)
        clientAttr.set_channel_mode(
          io.ioc.GLOBUS_IO_SECURE_CHANNEL_MODE_CLEAR)
        clientAttr.set_delegation_mode(
          io.ioc.GLOBUS_IO_SECURE_DELEGATION_MODE_LIMITED_PROXY)

      soc.connect(host, port, clientAttr)
      self.socket = soc
      self.sfile = soc.makefile("rw")

    finally:
      sys.stdout = sys.__stdout__
      sys.stderr = sys.__stderr__
      f.close()


  def __disconnect__(self):
    """
    Disconnect from the LDBD Server.

    @return: None
    """
    try:
      self.socket.shutdown(2)
    except:
      pass

  def __response__(self):
    """
    Read the response sent back by the LDBD Server. Parse out the
    return code with 0 for success and non-zero for error, and then
    the list of strings representing the returned result(s).

    @return: tuple containing the integer error code and the list of 
    strings representing the output from the server
    """
    f = self.sfile
       
    response = ""

    # Read in 512 byte chunks until there is nothing left to read.
    # This blocks until the socket is ready for reading and until
    # 512 bytes are received. If the message is less then 512 bytes
    # this will block until the server closes the socket. Since
    # the server always shuts down the socket after sending its
    # reply this should continue to work for now.
    while 1: 
      input = f.read(size = 512, waitForBytes = 512)
      response += input

      if len(input) < 512: break

    # the response from the server must always end in a null byte
    try:
      if response[-1] != '\0':
        msg = "Bad server reponse format. Contact server administrator."
        raise LDBDClientException, msg
    except:
      msg = "Connection refused. The server may be down or you may not have" + \
        "authorization to access this server. Contact server administrator."
      raise LDBDClientException, msg

    # delete the last \0 before splitting into strings
    response = response[0:-1]

    try:
      stringList = response.split('\0')
      code = int(stringList[0])
      output = stringList[1:]
    except Exception, e:
      msg = "Error parsing response from server : %s" % e
      try:
        f.close()
      except:
        pass
      raise LDBDClientException, msg


    f.close()

    return code, output

  def ping(self):
    """
    Ping the LDBD Server and return any message received back as a string.

    @return: message received (may be empty) from LDBD Server as a string
    """

    msg = "PING\0"
    self.sfile.write(msg)

    ret, output = self.__response__()
    reply = str(output[0])

    if ret:
      msg = "Error pinging server %d:%s" % (ret, reply)
      raise LDBDClientException, msg

    return reply

  def query(self,sql):
    """
    Execute an SQL query on the server and fetch the resulting XML file
    back.

    @return: message received (may be empty) from LDBD Server as a string
    """

    msg = "QUERY\0" + sql + "\0"
    self.sfile.write(msg)

    ret, output = self.__response__()
    reply = str(output[0])

    if ret:
      msg = "Error executing query on server %d:%s" % (ret, reply)
      raise LDBDClientException, msg

    return reply

  def insert(self,xmltext):
    """
    Insert the LIGO_LW metadata in the xmltext string into the database.

    @return: message received (may be empty) from LDBD Server as a string
    """

    msg = "INSERT\0" + xmltext + "\0"
    self.sfile.write(msg)

    ret, output = self.__response__()
    reply = str(output[0])

    if ret:
      msg = "Error executing insert on server %d:%s" % (ret, reply)
      raise LDBDClientException, msg

    return reply

  def insertmap(self,xmltext,lfnpfn_dict):
    """
    Insert the LIGO_LW metadata in the xmltext string into the database.

    @return: message received (may be empty) from LDBD Server as a string
    """

    pmsg = cPickle.dumps(lfnpfn_dict)

    msg = "INSERTMAP\0" + xmltext + "\0" + pmsg + "\0"
    self.sfile.write(msg)

    ret, output = self.__response__()
    reply = str(output[0])

    if ret:
      msg = "Error executing insert on server %d:%s" % (ret, reply)
      raise LDBDClientException, msg

    return reply

  def insertdmt(self,xmltext):
    """
    Insert the LIGO_LW metadata in the xmltext string into the database.

    @return: message received (may be empty) from LDBD Server as a string
    """

    msg = "INSERTDMT\0" + xmltext + "\0"
    self.sfile.write(msg)

    ret, output = self.__response__()
    reply = str(output[0])

    if ret:
      msg = "Error executing insert on server %d:%s" % (ret, reply)
      raise LDBDClientException, msg

    return reply
