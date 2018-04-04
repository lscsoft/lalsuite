"""
The LDBDWServer module provides an API for responding to request from the LDBDWClient by connecting to the DB2 database.

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
import cPickle
import logging
import logging.handlers
import ConfigParser
import M2Crypto
import selector
import cjson
import simplejson

from glue import ldbd

# default configuration
configuration = {
  'gridmap' : '/etc/ldbd/grid-mapfile',
  'gridmap_insert' : '/etc/ldbd/grid-mapfile',
  'gridmap_insertmap' : '/etc/ldbd/grid-mapfile',
  'gridmap_insertdmt' : '/etc/ldbd/grid-mapfile',
  'dbname' : 'ldbd_tst',
  'dbuser' : 'grid',
  'dbpasswd' : '',
  'logfile' : '/var/log/ldbd/ldbdd.log',
  'logmaxbytes' : 1024 * 1024 * 1,
  'max_client_byte_string' : 1048576,
  'logbackupcount' : 5,
  'loglevel' : 'DEBUG'
  }

# grab configuration from file
myConfigParser = ConfigParser.ConfigParser()

try:
    ldbdserver_ini = os.environ.get('LDBD_CONFIG_DIR', '/usr3/ldbd/etc')
    config_file = os.path.join(ldbdserver_ini, "ldbdserver.ini")
    myConfigParser.read(config_file)
except Exception, e:
    sys.stder.write("Error: unable to read configuration file : %s\n" % config_file)
    sys.exit(1)

for k in configuration.keys():
    try:
        value = myConfigParser.get('ldbdd',k)
    except ConfigParser.NoOptionError, e:
        sys.stderr.write("Error: missing configuration option %s: %s\n" % (k, e))
        sys.exit(1)
    try:
        configuration[k] = eval(value)
    except Exception, e:
        configuration[k] = value

# set up logging
logFilePath = configuration['logfile']
logMaxBytes = configuration['logmaxbytes']
logBackupCount = configuration['logbackupcount']
logLevel = configuration['loglevel']

# set up logging for the root
rootlogger = logging.getLogger('LDBD')

handler = logging.handlers.RotatingFileHandler(logFilePath, 'a', logMaxBytes, logBackupCount)
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)

rootlogger.addHandler(handler)
rootlogger.setLevel(eval("logging.%s" % logLevel))

# create logger to use with correct prefix
logger = logging.getLogger('LDBD.Server')
logger.info("LDBD server started")

# initialize the database hash table
dbobj = ldbd.LIGOMetadataDatabase(configuration['dbname'])
max_bytes = configuration['max_client_byte_string']

# create the xml parser
xmlparser = pyRXP.Parser()

def dtd_uri_callback(uri):
  if uri == 'http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt':
    # if the XML file contants a http pointer to the ligolw DTD at CIT then
    # return a local copy to avoid any network problems
    return 'file://localhost' + os.path.join( os.environ["GLUE_PREFIX"],
      'etc/ligolw_dtd.txt' )
  else:
    # otherwise just use the uri in the file
    return uri


# use a local copy of the DTD, if one is available
try:
    GLUE_PREFIX = os.environ["GLUE_PREFIX"]
    xmlparser.eoCB = dtd_uri_callback
    logger.info("Using local DTD in " + 
      'file://localhost' + os.path.join( GLUE_PREFIX, 'etc') )
except KeyError:
    logger.warning('GLUE_PREFIX not set, unable to use local DTD') 

# initialize dictionaries for the dmt processes and segments definers
dmt_proc_dict = {}
dmt_seg_def_dict = {}
creator_db = None

class Server(object):
  """
  """
  def __init__(self):
    """
    """

    # create a selector to use for dispatching
    mySelector = selector.Selector()
    self.mySelector = mySelector

    # for convenience
    global configuration
    self.configuration = configuration

    # define dispatches
    mySelector.add('/ldbd/ping[.{format}]', POST=self.ping, GET=self.ping)
    mySelector.add('/ldbd/query[.{format}]', POST=self.query)
    mySelector.add('/ldbd/insert[.{format}]', POST=self.insert)
    mySelector.add('/ldbd/insertmap[.{format}]', POST=self.insertmap)
    mySelector.add('/ldbd/insertdmt[.{format}]', POST=self.insertdmt)

  def __call__(self, environ, start_response):
    """
    """

    # necessary for version 0.8.11 of selector because it expects to
    # use PATH_INFO but mod_wsgi is using the more modern and standard
    # SCRIPT_NAME
    environ['PATH_INFO'] = environ['SCRIPT_NAME'] + environ['PATH_INFO']

    # use selector to choose and execute the correct function 
    return self.mySelector(environ, start_response)

  def checkAuthorizationGridMap(self, environ, mapfile = None):
        """
        """
        # first make sure the client was verified
        verification = environ['SSL_CLIENT_VERIFY']
        if verification != 'SUCCESS':
            return (False, None)

        # grab the client cert
        clientCert = M2Crypto.X509.load_cert_string(environ['SSL_CLIENT_CERT'])
        clientSubject = clientCert.get_subject().__str__()
        logger.debug("The client cert subject is %s" % clientSubject)

        # see if the client cert is a proxy
        try:
            clientCert.get_ext("proxyCertInfo")
            clientCertProxy = True
            logger.debug("The client cert is a proxy")
        except LookupError:
            clientCertProxy = False
            logger.debug("The client cert is not a proxy")

        if not clientCertProxy:
            subject = clientSubject

        else:
            # find the number of certificates in the client chain
            max = None
            for e in environ.keys():
                m = re.search('SSL_CLIENT_CERT_CHAIN_(\d+)', e)
                if m:
                    logger.debug("Found SSL_CLIENT_CERT_CHAIN_%s" % m.group(1))
                    n = int(m.group(1))
                    if n > max: max = n

            certChainLength = max + 1

            logger.debug("There are %d certs in the chain" % certChainLength)
                
            # grab each certificate in the chain
            chainCerts = [ M2Crypto.X509.load_cert_string(environ['SSL_CLIENT_CERT_CHAIN_%d' % i]) for i in range(certChainLength) ]

            # walk the chain and find the first end entity certificate
            logger.debug("Walking the cert chain now...")
            for c in chainCerts:
                s = c.get_subject().__str__()
                logger.debug("Chain cert subject is %s" % s)
                try:
                    c.get_ext("proxyCertInfo")
                    logger.debug("Chain cert %s is a proxy" % s)
                except LookupError:
                    logger.debug("Chain cert %s is not a proxy" % s)
                    break
            subject = s

        logger.debug("Authorizing against %s" % subject)

        # parse the grid-mapfile and see if subject in it
        authorized = False
        if not mapfile:
            mapfile = configuration['gridmap']
        try:
            g = GridMap(mapfile)
            if g.has_key(subject):
                authorized = True
            else:
                authorized = False
        except Exception, e:
            logger.error("Unable to check authorization in grid-mapfile %s: %s" % (mapfile, e))

        return (authorized, subject)

  def ping(self, environ, start_response):
    """
    """
    logger.debug("Method ping called")

    # determine protocol
    try:
        protocol = cjson.decode(environ['wsgi.input'].read())
    except Exception, e:
        if (str(e)).strip() == "empty JSON description":
           logger.debug("No protocol given, request may come from Web GUI")
           protocol = "https"
        else:
           start_response("400 Bad Request", [('Content-type', 'text/plain')])
           msg = "400 Bad Request"
           logger.debug("Error decoding input: %s" % e)
           return [ msg ]


    if protocol == "https": 
      # use generic grid-mapfile for ping operation
      mapfile = self.configuration['gridmap']

      # check authorization
      authorized, subject = self.checkAuthorizationGridMap(environ, mapfile)
      if not authorized:
        start_response("401 Unauthorized", [('Content-type', 'text/plain')])
        msg = "401 Unauthorized\n\nSubject %s is not authorized for method ping" % subject
        logger.info("Subject %s is not authorized for method ping" % subject)
        return [ msg ]
      else:
        logger.info("Subject %s is authorized for method ping" % subject)
    else:
      logger.info("Method ping is being called from internal network")


    try:
      hostname = socket.getfqdn()
      msg = "LDBD server at %s is alive" % hostname
    except Exception, e:
      msg = "LDBD server is alive" 

    format = environ['wsgiorg.routing_args'][1]['format']
    if format == 'json':
        result = cjson.encode(msg)
        contentType = 'application/json'
    else:
        result = msg
        contentType = 'text/plain'

    header = [('Content-Type', contentType)]
    start_response("200 OK", header)

    return [ result ]

  def query(self, environ, start_response):
    """
    """
    global xmlparser, dbobj

    logger.debug("Method query called")

    # determine protocol
    try:
        protocol, querystr = (cjson.decode(environ['wsgi.input'].read())).split(":")
    except Exception, e:
        start_response("400 Bad Request", [('Content-type', 'text/plain')])
        msg = "400 Bad Request"
        logger.debug("Error decoding input: %s" % e)
        return [ msg ]


    if protocol == "https":
      # use generic grid-mapfile for query operation
      mapfile = self.configuration['gridmap']

      # check authorization
      authorized, subject = self.checkAuthorizationGridMap(environ, mapfile)
      if not authorized:
        start_response("401 Unauthorized", [('Content-type', 'text/plain')])
        msg = "401 Unauthorized\n\nSubject %s is not authorized for method query" % subject
        logger.info("Subject %s is not authorized for method query" % subject)
        return [ msg ]
      else:
        logger.info("Subject %s is authorized for method query" % subject)
    else:
      logger.info("Method QUERY is being called from internal network")


    # pick off the format input type
    format = environ['wsgiorg.routing_args'][1]['format']

    # for a POST the format must be JSON
    if format != 'json':
        start_response("400 Bad Request", [('Content-type', 'text/plain')])
        msg = "400 Bad Request\n\nformat must be 'json' for POST operation"
        return [ msg ]

    logger.debug("Method query called with '%s'" % querystr)

    try:
      # create a ligo metadata object
      lwtparser = ldbd.LIGOLwParser()
      ligomd = ldbd.LIGOMetadata(xmlparser,lwtparser,dbobj)

      # execute the query
      rowcount = ligomd.select(querystr)

      # convert the result to xml
      result = ligomd.xml()

      logger.debug("Method query: %d rows returned" % rowcount)
    except Exception, e:
      start_response("500 Internal Server Error", [('Content-type', 'text/plain')])
      msg = "500 Internal Server Error\n\n%s" % e
      logger.error(msg)
      return [ msg ]

    try:
      del ligomd
      del lwtparser
    except Exception, e:
      logger.error("Error deleting metadata object in method query: %s" % e)

    # encode the result
    result = cjson.encode(result)
    
    # return the result
    header = [('Content-Type', 'application/json')]
    start_response("200 OK", header)

    return [ result ]

  def insert(self, environ, start_response):
    """
    """
    global xmlparser, dbobj

    logger.debug("Method insert called")

    # use specific grid-mapfile for insert operation
    mapfile = self.configuration['gridmap_insert']

    # check authorization
    authorized, subject = self.checkAuthorizationGridMap(environ, mapfile)
    if not authorized:
        start_response("401 Unauthorized", [('Content-type', 'text/plain')])
        msg = "401 Unauthorized\n\nSubject %s is not authorized for method insert" % subject
        logger.info("Subject %s is not authorized for method insert" % subject)
        return [ msg ]
    else:
        logger.info("Subject %s is authorized for method insert" % subject)

    # pick off the format input type
    format = environ['wsgiorg.routing_args'][1]['format']

    # for a POST the format must be JSON
    if format != 'json':
        start_response("400 Bad Request", [('Content-type', 'text/plain')])
        msg = "400 Bad Request\n\nformat must be 'json' for POST operation"
        return [ msg ]

    # read the incoming payload
    try:
        #import simplejson  (moved to top)
        wsgiIn=environ['wsgi.input'].read()
        inputString=simplejson.loads(wsgiIn)
        #inputString = cjson.decode(environ['wsgi.input'].read())
    except Exception, e:
        start_response("400 Bad Request", [('Content-type', 'text/plain')])
        msg = "400 Bad Request"
        logger.debug("Error decoding input: %s" % e)
        return [ msg ]

    logger.debug("Method insert called with '%s'" % inputString)

    try:
      # create a ligo metadata object
      lwtparser = ldbd.LIGOLwParser()
      ligomd = ldbd.LIGOMetadata(xmlparser,lwtparser,dbobj)

      # parse the input string into a metadata object
      ligomd.parse(inputString)

      # add a gridcert table to this request containing the users dn
      ligomd.set_dn(subject)

      # insert the metadata into the database
      result = str(ligomd.insert())

      logger.info("Method insert: %s rows affected by insert" % result)
    except Exception, e:
      start_response("500 Internal Server Error", [('Content-type', 'text/plain')])
      msg = "500 Internal Server Error\n\n%s" % e
      logger.error(msg)
      return [ msg ]

    try:
      del ligomd
      del lwtparser
    except Exception, e:
      logger.error("Error deleting metadata object in method query: %s" % e)

    # encode the result
    result = cjson.encode(result)
    
    # return the result
    header = [('Content-Type', 'application/json')]
    start_response("200 OK", header)

    return [ result ]
    
  def insertmap(self, environ, start_response):
    """
    """
    logger.debug("Method insertmap called")

    # use specific grid-mapfile for insert operation
    mapfile = self.configuration['gridmap_insertmap']

    # check authorization
    authorized, subject = self.checkAuthorizationGridMap(environ, mapfile)
    if not authorized:
        start_response("401 Unauthorized", [('Content-type', 'text/plain')])
        msg = "401 Unauthorized\n\nSubject %s is not authorized for method insertmap" % subject
        logger.info("Subject %s is not authorized for method insertmap" % subject)
        return [ msg ]
    else:
        logger.info("Subject %s is authorized for method insertmap" % subject)

    # pick off the format input type
    format = environ['wsgiorg.routing_args'][1]['format']

    # for a POST the format must be JSON
    if format != 'json':
        start_response("400 Bad Request", [('Content-type', 'text/plain')])
        msg = "400 Bad Request\n\nformat must be 'json' for POST operation"
        return [ msg ]

    start_response("500 Internal Server Error", [('Content-type', 'text/plain')])
    msg = "500 Internal Server Error\n\nserver is not initialized for RLS connections"
    logger.error(msg)
    return [ msg ]

  def insertdmt(self, environ, start_response):
    """
    """
    global xmlparser, dbobj
    global dmt_proc_dict, dmt_seg_def_dict, creator_db
    proc_key = {}
    known_proc = {}
    seg_def_key = {}

    logger.debug( "Method dmtinsert called." )
    logger.debug( "Known processes %s, " % str(dmt_proc_dict) )
    logger.debug( "Known segment_definers %s" % str(dmt_seg_def_dict) )

    # use specific grid-mapfile for insertdmt operation
    mapfile = self.configuration['gridmap_insertdmt']

    # check authorization
    authorized, subject = self.checkAuthorizationGridMap(environ, mapfile)
    if not authorized:
        start_response("401 Unauthorized", [('Content-type', 'text/plain')])
        msg = "401 Unauthorized\n\nSubject %s is not authorized for method insertdmt" % subject
        logger.info("Subject %s is not authorized for method insertdmt" % subject)
        return [ msg ]
    else:
        logger.info("Subject %s is authorized for method insertdmt" % subject)

    # pick off the format input type
    format = environ['wsgiorg.routing_args'][1]['format']

    # for a POST the format must be JSON
    if format != 'json':
        start_response("400 Bad Request", [('Content-type', 'text/plain')])
        msg = "400 Bad Request\n\nformat must be 'json' for POST operation"
        return [ msg ]

    # read the incoming payload
    try:
        inputString = cjson.decode(environ['wsgi.input'].read())
    except Exception, e:
        start_response("400 Bad Request", [('Content-type', 'text/plain')])
        msg = "400 Bad Request"
        logger.debug("Error decoding input: %s" % e)
        return [ msg ]

    logger.debug("Method insertdmt called with '%s'" % inputString)

    try:
      # create a ligo metadata object
      lwtparser = ldbd.LIGOLwParser()
      ligomd = ldbd.LIGOMetadata(xmlparser,lwtparser,dbobj)

      # parse the input string into a metadata object
      logger.debug("parsing xml data")
      ligomd.parse(inputString)

      # store the users dn in the process table
      ligomd.set_dn(subject)

      # determine the local creator_db number
      if creator_db is None:
        sql = "SELECT DEFAULT FROM SYSCAT.COLUMNS WHERE "
        sql += "TABNAME = 'PROCESS' AND COLNAME = 'CREATOR_DB'"
        ligomd.curs.execute(sql)
        creator_db = ligomd.curs.fetchone()[0]

      # determine the locations of columns we need in the process table
      process_cols = ligomd.table['process']['orderedcol']
      node_col = process_cols.index('node')
      prog_col = process_cols.index('program')
      upid_col = process_cols.index('unix_procid')
      start_col = process_cols.index('start_time')
      end_col = process_cols.index('end_time')
      pid_col = process_cols.index('process_id')

      # determine and remove known entries from the process table
      rmv_idx = []
      for row_idx,row in enumerate(ligomd.table['process']['stream']):
        uniq_proc = (row[node_col],row[prog_col],row[upid_col],row[start_col])
        try:
          proc_key[str(row[pid_col])] = dmt_proc_dict[uniq_proc]
          known_proc[str(dmt_proc_dict[uniq_proc])] = row[end_col]
          rmv_idx.append(row_idx)
        except KeyError:
          # we know nothing about this process, so query the database
          sql = "SELECT BLOB(process_id) FROM process WHERE "
          sql += "creator_db = " + str(creator_db) + " AND "
          sql += "node = '" + row[node_col] + "' AND "
          sql += "program = '" + row[prog_col] + "' AND "
          sql += "unix_procid = " + str(row[upid_col]) + " AND "
          sql += "start_time = " + str(row[start_col])
          ligomd.curs.execute(sql)
          db_proc_ids = ligomd.curs.fetchall()
          if len(db_proc_ids) == 0:
            # this is a new process with no existing entry
            dmt_proc_dict[uniq_proc] = row[pid_col]
          elif len(db_proc_ids) == 1:
            # the process_id exists in the database so use that insted
            dmt_proc_dict[uniq_proc] = db_proc_ids[0][0]
            proc_key[str(row[pid_col])] = dmt_proc_dict[uniq_proc]
            known_proc[str(dmt_proc_dict[uniq_proc])] = row[end_col]
            rmv_idx.append(row_idx)
          else:
            # multiple entries for this process, needs human assistance
            raise ServerHandlerException("multiple entries for dmt process")

      # delete the duplicate processs rows and clear the table if necessary
      newstream = []
      for row_idx,row in enumerate(ligomd.table['process']['stream']):
        try:
          rmv_idx.index(row_idx)
        except ValueError:
          newstream.append(row)
      ligomd.table['process']['stream'] = newstream
      if len(ligomd.table['process']['stream']) == 0:
        del ligomd.table['process']

      # delete the duplicate process_params rows and clear the table if necessary
      # (the DMT does not write a process_params table, so check for one first)
      if ligomd.table.has_key('process_params'):
        ppid_col = ligomd.table['process_params']['orderedcol'].index('process_id')
        newstream = []
        for row_idx,row in enumerate(ligomd.table['process_params']['stream']):
          # if the process_id in this row is known, delete (i.e. don't copy) it
          try:
            proc_key[str(row[ppid_col])]
          except KeyError:
            newstream.append(row)
        ligomd.table['process_params']['stream'] = newstream
        if len(ligomd.table['process_params']['stream']) == 0:
          del ligomd.table['process_params']

      # turn the known process_id binary for this insert into ascii
      for pid in known_proc.keys():
        pid_str = "x'"
        for ch in pid:
          pid_str += "%02x" % ord(ch)
        pid_str += "'"
        known_proc[pid] = (pid_str, known_proc[pid])

      # determine the locations of columns we need in the segment_definer table
      seg_def_cols = ligomd.table['segment_definer']['orderedcol']
      ifos_col = seg_def_cols.index('ifos')
      name_col = seg_def_cols.index('name')
      vers_col = seg_def_cols.index('version')
      sdid_col = seg_def_cols.index('segment_def_id')

      # determine and remove known entries in the segment_definer table
      rmv_idx = []
      for row_idx,row in enumerate(ligomd.table['segment_definer']['stream']):
        uniq_def = (row[ifos_col],row[name_col],row[vers_col])
        try:
          seg_def_key[str(row[sdid_col])] = dmt_seg_def_dict[uniq_def]
          rmv_idx.append(row_idx)
        except KeyError:
          # we know nothing about this segment_definer, so query the database
          sql = "SELECT BLOB(segment_def_id) FROM segment_definer WHERE "
          sql += "creator_db = " + str(creator_db) + " AND "
          sql += "ifos = '" + row[ifos_col] + "' AND "
          sql += "name = '" + row[name_col] + "' AND "
          sql += "version = " + str(row[vers_col])
          ligomd.curs.execute(sql)
          db_seg_def_id = ligomd.curs.fetchall()
          if len(db_seg_def_id) == 0:
            # this is a new segment_defintion with no existing entry
            dmt_seg_def_dict[uniq_def] = row[sdid_col]
          else:
            dmt_seg_def_dict[uniq_def] = db_seg_def_id[0][0]
            seg_def_key[str(row[sdid_col])] = dmt_seg_def_dict[uniq_def]
            rmv_idx.append(row_idx)

      # delete the necessary rows. if the table is empty, delete it
      newstream = []
      for row_idx,row in enumerate(ligomd.table['segment_definer']['stream']):
        try:
          rmv_idx.index(row_idx)
        except ValueError:
          newstream.append(row)
      ligomd.table['segment_definer']['stream'] = newstream
      if len(ligomd.table['segment_definer']['stream']) == 0:
        del ligomd.table['segment_definer']

      # now update the values in the xml with the values we know about
      for tabname in ligomd.table.keys():
        table = ligomd.table[tabname]
        if tabname == 'process':
          # we do nothing to the process table
          pass
        elif tabname == 'segment' or tabname == 'segment_summary':
          # we need to update the process_id and the segment_def_id columns
          pid_col = table['orderedcol'].index('process_id')
          sdid_col = table['orderedcol'].index('segment_def_id')
          row_idx = 0
          for row in table['stream']:
            try:
              repl_pid = proc_key[str(row[pid_col])]
            except KeyError:
              repl_pid = row[pid_col]
            try:
              repl_sdid = seg_def_key[str(row[sdid_col])]
            except KeyError:
              repl_sdid = row[sdid_col]
            row = list(row)
            row[pid_col] = repl_pid
            row[sdid_col] = repl_sdid
            table['stream'][row_idx] = tuple(row)
            row_idx += 1
        else:
          # we just need to update the process_id column
          pid_col = table['orderedcol'].index('process_id')
          row_idx = 0
          for row in table['stream']:
            try:
              repl_pid = proc_key[str(row[pid_col])]
              row = list(row)
              row[pid_col] = repl_pid
              table['stream'][row_idx] = tuple(row)
            except KeyError:
              pass
            row_idx += 1

      # insert the metadata into the database
      logger.debug("inserting xml data")
      result = str(ligomd.insert())
      logger.debug("insertion complete")

      # update the end time of known processes in the process table
      for pid in known_proc.keys():
        # first check to see if we are backfilling missing segments
        sql = "SELECT end_time,domain FROM process "
        sql += " WHERE process_id = " + known_proc[pid][0]
        ligomd.curs.execute(sql)
        last_end_time = ligomd.curs.fetchone()

        # check the dn in the row we are about to update matches the users dn
        dn = last_end_time[1].strip()
        if subject != dn:
          msg = "\"%s\" does not match dn in existing row entries: " % subject
          msg += "%s (process_id %s)" % (dn, known_proc[pid][0])
          logger.warn(msg)
        else:
          logger.debug('"%s" updating process_id %s' % (dn, known_proc[pid][0]))

        if int(known_proc[pid][1]) <= int(last_end_time[0]):
          logger.debug("Backfilling missing segments for process_id " +
            known_proc[pid][0] + " not updating end_time")
        else:
          # if we are not backfilling, update the end_time of the process
          sql = "UPDATE process SET end_time = " + str(known_proc[pid][1])
          sql += " WHERE process_id = " + known_proc[pid][0]
          sql += " AND end_time < " + str(known_proc[pid][1])
          ligomd.curs.execute(sql)
      ligomd.dbcon.commit()

      logger.info("Method insert: %s rows affected by insert" % result)

    except Exception, e:
      start_response("500 Internal Server Error", [('Content-type', 'text/plain')])
      msg = "500 Internal Server Error\n\n%s" % e
      logger.error(msg)
      return [ msg ]

    try:
      del ligomd
      del lwtparser
      del known_proc
      del seg_def_key
      del proc_key
    except Exception, e:
      logger.error("Error deleting metadata object in method insertdmt: %s" % e)

    # encode the result
    result = cjson.encode(result)
    
    # return the result
    header = [('Content-Type', 'application/json')]
    start_response("200 OK", header)

    return [ result ]

class GridMapError(exceptions.Exception):
    """
    Raised for errors in GridMap class
    """
    pass

class GridMap(dict):
    """
    """
    def __init__(self, path=None):
        """
        """
        self.path = path

        # initialize the base class
        dict.__init__(self)

        if not path:
            msg = "No path for grid-mapfile"
            raise GridMapError(msg)

    def parse(self):
        """
        """
        # clear any existing entries
        self.clear()

        # parse the grid-mapfile
        try:
            f = open(self.path, 'r')
        except Exception, e:
            msg = "Unable to open %s for reading: %s" % (self.path, e)
            raise GridMapError(msg)

        for line in f:
            s = line.strip().split('"')
            if len(s) == 3: subject = s[1]
            elif len(s) == 2: subject = s[1]
            elif len(s) == 1: subject = s[0]
            else:
                msg = "Error parsing line %s" % line
                raise GridMapError(msg)

            dict.__setitem__(self, subject, 1)

        f.close()

    def __getitem__(self, key):
        """
        """
        self.parse()

        return dict.__getitem__(self, key)

    def has_key(self, key):
        """
        """
        self.parse()

        return dict.has_key(self, key)

    def keys(self):
        """
        """
        self.parse()

        return dict.keys(self)

