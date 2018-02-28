"""
The LDRdataFindClient module provides an API for connecting to
and making requests of a LDRdataFindServer.

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
__version__ = git_version.id

import sys
import os
import exceptions

from pyGlobus import io
from pyGlobus import security


def version():
        return __version__


class LSCdataFindClientException(exceptions.Exception):
        """
        Exceptions raised by the classes and methods in this client
        will be instances of this class.
        """
        def __init__(self, args=None):
                """
                Create an instance of this class, ie. an LSCdataFindException.

                @param args:

                @return: Instance of class LSCdataFindException
                """
                self.args = args


class LDRdataFindClientException(exceptions.Exception):
        """
        Exceptions raised by the classes and methods in this module
        will be instances of this class.
        """
        def __init__(self, args=None):
                """
                Create an instance of this class, ie. an LDRdataFindClientException.

                @param args:

                @return: Instance of class LDRdataFindClientException
                """
                self.args = args


class pfnlist(list):
        pass


class lfnlist(list):
        pass


class LDRMetadataQuery(object):
        """
        """
        def __init__(self, attrQuery=None):
                """
                """
                self.attrQuery = attrQuery
                self.regex = None
                self.offset = 0
                self.limit = None
                self.sortAttr = None
                self.sortOrder = None

        def set_query(self, query):
                self.attrQuery = query

        def set_regex_filter(self, filter):
                self.regex = filter

        def set_offset(self, offset):
                self.offset = offset

        def set_limit(self, limit):
                self.limit = limit

        def set_sort_attribute(self, attribute):
                self.sortAttr = attribute

        def set_sort_order(self, sortOrder):
                self.sortOrder = sortOrder

        def __verify(self):
                # we should insure here that
                # - attrQuery is a string and well behaved
                # - regex is a string and well behaved
                # - offset is an integer >= 0
                # - limit is None or integer > 0
                # - sortAttr is None or string
                # - sortOrder is one of "ASC" or "DESC"
                pass

        def __str__(self):
                self._LDRMetadataQuery__verify()

                if not self.limit: 
                        limit = "-1"
                else:
                        limit = str(self.limit)

                if not self.sortAttr:
                        sortAttr = "NONE"
                else:
                        sortAttr = self.sortAttr

                if not self.regex:
                        regex = "NONE"
                else:
                        regex = self.regex

                if not self.sortOrder:
                        sortOrder = "ASC"
                else:
                        sortOrder = self.sortOrder
                        

                s = "%s\0%s\0%d\0%s\0%s\0%s\0" % (
                        self.attrQuery,
                        regex,
                        self.offset,
                        limit,
                        sortAttr,
                        sortOrder
                        )

                return s


class LDRdataFindClient(object):
        """
        Class that represents a client interacting with a LDRdataFindServer. It is expected
        that clients for interacting with a LDRdataFindServer will inherit from this base
        class.
        """
        def __init__(self, host, port, noproxy, disablehostauth):
                """
                Open a connection to a LDRdataFindServer and return an instance of
                class LDRdataFindClient. One of the public methods can then be 
                called to send a request to the server.

                @param host: the host on which the LDRdataFindServer runs
                @type host: string

                @param port: port on which the LDRdataFindServer listens
                @type port: integer

                @pararm noproxy: whether or not to connect with a proxy
                                 certificate.
                @type noproxy: boolean


                @return: Instance of LDRdataFindClient
                """
                # Connection to be performed just before
                # query is sent to server. So just store 
                # the necessary connection information
                self.host = host
                self.port = port
                self.noproxy = noproxy
                self.disablehostauth = disablehostauth


        def __del__(self):
                """
                Disconnect from the LDRdataFindServer.

                @return: None
                """
                self.__disconnect__()


        def __connect__(self):
                """
                Attempt to open a connection to the LDRdataFindServer
                using the 'host' and 'port' and expecting the server
                to identify itself with a corresponding host certificate.

                A IOException is raised if the connection cannot be made,
                but this is caught by the __init__ method above and 
                turned into a LDRdataFindClient exception.
        
                @param host: the host on which the LDRdataFindServer runs
                @type host: string

                @param port: port on which the LDRdataFindServer listens
                @type port: integer

                @param noproxy: whether or not to attempt connection with
                                a proxy certificate.
                @type noproxy: boolean

                @return: None
                """

                # remove the globus tcp port range environment variable if set
                try:
                        port_range = os.environ["GLOBUS_TCP_PORT_RANGE"]
                        os.environ["GLOBUS_TCP_PORT_RANGE"] = ""
                except:
                        pass

                

                # redirect stdout and stderror for now
                try:
                        f = open("/dev/null", "w")
                        sys.stdout = f
                        sys.stderr = f
                except:
                        pass

                try:
                        if self.noproxy:
                                # Try connecting without the proxy
                                # create TCPIOAttr instance and set authentication mode to be NONE
                                myAttr = io.TCPIOAttr()
                                myAttr.set_authentication_mode(io.ioc.GLOBUS_IO_SECURE_AUTHENTICATION_MODE_NONE)

                                authData = io.AuthData()
                                # set authorization, channel, and delegation modes
                                myAttr.set_authorization_mode(io.ioc.GLOBUS_IO_SECURE_AUTHORIZATION_MODE_NONE,authData)
                                myAttr.set_channel_mode(io.ioc.GLOBUS_IO_SECURE_CHANNEL_MODE_CLEAR)
                                myAttr.set_delegation_mode(io.ioc.GLOBUS_IO_SECURE_DELEGATION_MODE_NONE)
                        else:
                                # Try normal GSI-authenticated connection
                                # create TCPIOAttr instance and set authentication mode to be NONE
                                myAttr = io.TCPIOAttr()
                                myAttr.set_authentication_mode(io.ioc.GLOBUS_IO_SECURE_AUTHENTICATION_MODE_GSSAPI)

                                authData = io.AuthData()
                                # set authorization, channel, and delegation modes
                                if self.disablehostauth:
                                        # Disable the GSI hostname checking
                                        myAttr.set_authorization_mode(io.ioc.GLOBUS_IO_SECURE_AUTHORIZATION_MODE_NONE,authData)
                                        myAttr.set_channel_mode(io.ioc.GLOBUS_IO_SECURE_CHANNEL_MODE_CLEAR)
                                        myAttr.set_delegation_mode(io.ioc.GLOBUS_IO_SECURE_DELEGATION_MODE_NONE)
                                else:
                                        # Use the normal GSI hostname checking
                                        myAttr.set_authorization_mode(io.ioc.GLOBUS_IO_SECURE_AUTHORIZATION_MODE_HOST,authData)
                                        myAttr.set_channel_mode(io.ioc.GLOBUS_IO_SECURE_CHANNEL_MODE_CLEAR)
                                        myAttr.set_delegation_mode(io.ioc.GLOBUS_IO_SECURE_DELEGATION_MODE_LIMITED_PROXY)

                        # fi

                        # create socket instance and attempt to connect
                        s = io.GSITCPSocket()
                        s.connect(self.host, self.port, myAttr)
                        self.socket = s
                        self.sfile = s.makefile("rw")


                finally:
                        sys.stdout = sys.__stdout__
                        sys.stderr = sys.__stderr__
                        f.close()

        def __disconnect__(self):
                """
                Disconnect from the LDRdataFindServer.

                @return: None
                """
                try:
                        self.socket.shutdown(2)
                except:
                        pass

        def __response__(self):
                """
                Read the response sent back by the LDRdataFindServer. Parse out the
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
                                raise LDRdataFindClientException(msg)
                except:
                        msg = "Connection refused. The server may be down or you may not have authorization to access this server. Contact server administrator."
                        raise LDRdataFindClientException(msg)

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
                        raise LDRdataFindClientException(msg)

                f.close()

                return code, output

        def ping(self):
                """
                Ping the LDRdataFindServer and return any message received back as
                a string.

                @return: message received (may be empty) from LDRdataFindServer as a string
                """

                msg = "PING\0"
                self.__connect__()
                self.sfile.write(msg)

                ret, output = self.__response__()
                reply = str(output[0])

                if ret:
                        msg = "Error pinging server %d:%s" % (ret, reply)
                        raise LDRdataFindClientException(msg)

                return reply


        def distinctAttrValues(self, attr):
                """
                Query LDRdataFindServer metadata tables for the distince values of an attribute
                and return the values as a list.

                Note that the values will always be returned as strings. Any conversion
                should be done by the client.

                @param attr: name of attribute for which to find distinct values

                @return: list of strings representing the distinct values that the attribute
                        has in the metadata tables
                
                """

                # check argument
                if not isinstance(attr, str):
                        msg = "Argument 'attr' must be a string"
                        raise LDRdataFindClientException(msg)
        
                msg = "DISTINCT\0%s\0" % str(attr)
                self.__connect__()
                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for distinct values of attributes: %s" % str(output)
                        raise LDRdataFindClientException(msg)

                return output
        
        def pfnQuery(self, lfn):
                """
                Query LDRdataFindServer to find the PFN(s) associated with a LFN and return them
                as a list. More then one PFN may be returned.

                @param lfn: the LFN with which to query the LDRdataFindServer

                @return: list of strings representing the PFN(s) for the LFN
                """

                # check argument
                if not isinstance(lfn, str):
                        msg = "Argument 'lfn' must be a string"
                        raise LDRdataFindClientException(msg)

                msg = "LFNPFN\0%s\0" % str(lfn)
                self.__connect__()
                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for PFN with LFN %s : %s" % (str(lfn), str(output))
                        raise LDRdataFindClientException(msg)

                return output
                
        def timeQuery(self, mytype, site, start, end, strict):
                """
                Query LDRdataFindServer for time ranges for a particular frameType.
                Optionally supprts gpsStart and gpsEnd, these will be processed server-side.
                
                @param mytype: frame type
                @param site: observatory
                @param start: gps start time
                @param stop: gps end time
                @param strict: strict query flag
                
                @return: output of query.
                """
                if not mytype:
                        msg = "A frame type, --type, must be specified."
                        raise LSCdataFindClientException(msg)
                if int(end) < int(start):
                        msg = "Supplied start time, %s, is later than supplied end time, %s" % (start,end)
                        raise LSCdataFindClientException(msg)

                
                #Construct RPC
                rpc = "SEGMENT\0type\0%s\0observatory\0%s\0" % (mytype,site)
                
                if start:
                        rpc += "gpsStart\0%s\0" % (start,)
                if end:
                        rpc += "gpsEnd\0%s\0" % (end,)
                if strict:
                        rpc += "strict\0"
                
                self.__connect__()
                self.sfile.write(rpc)
                ret,output = self.__response__()
                
                if ret:
                        msg = "Error querying LDRdataFindServer for times of frame type, %s: %s" % (str(mytype),str(output))
                        raise LSCdataFindClientException(msg)
                return output

        def lfnQueryWithMetadata(self, queryList):
                """
                Query LDRdataFindServer to find the LFN(s) with the appropriate metadata values.

                @param queryList: list of instances of the LDRMetadataQuery class, each describing
                                  a query to be done
                        
                @return: list of strings representing the union of all LFN(s) found from the queries
                """

                # check arguments
                if not isinstance(queryList, list):
                        msg = "Argument must be a list of instances of LDRMetadataQuery"
                        raise LDRdataFindClientException(msg)
                for query in queryList:
                        if not isinstance(query, LDRMetadataQuery):
                                msg = "Argument must be an instance of LDRMetadataQuery"
                                raise LDRdataFindClientException(msg)

                # prepare the messange to send down the socket
                msg = "METALFN\0"
                for q in queryList:
                        msg += "%s" % str(q)
                self.__connect__()
                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for LFNs: %s" % (str(output[0],))
                        raise LDRdataFindClientException(msg)


                return output
                
        def pfnQueryWithMetadata(self, queryList):
                """
                Query LDRdataFindServer to find the PFNs(s) for LFN(s) with the appropriate 
                metadata values.

                @param queryList: list of instances of the LDRMetadataQuery class, each describing
                                  a query to be done

                @return: list of strings representing the PFN(s) found
                """

                # check arguments
                if not isinstance(queryList, list):
                        msg = "Argument must be a list of instances of LDRMetadataQuery"
                        raise LDRdataFindClientException(msg)
                for query in queryList:
                        if not isinstance(query, LDRMetadataQuery):
                                msg = "Argument must be an instance of LDRMetadataQuery"
                                raise LDRdataFindClientException(msg)

                msg = "METAPFN\0" 
                for q in queryList:
                        msg += "%s" % str(q)
                self.__connect__()
                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for PFNs: %s" % (str(output[0],))
                        raise LDRdataFindClientException(msg)

                return output

        def pfnQueryWithMetadataRegExp(self, sql, rexp, offset=None, number=None):
                """
                Query LDRdataFindServer to find the PFNs(s) for LFN(s) with the appropriate 
                metadata values and return those that match a regular expression.

                @param sql: clause that will be part of a SQL query done to find the LFN(s)
                        with the appropriate metadata values. A an example would be
                        
                        gpsStart <= '777777777' AND gpsEnd >= '666666666' AND instrument = 'H' AND runTag = 'S2'

                @param rexp: regular expression against which to match found PFNs

                @param offset: the offset into the list of matching PFNs at which to begin 
                        returning results

                @param number: the total number of PFNs to return

                @return: list of strings representing the PFN(s) found
                """

                # check arguments
                if not isinstance(sql, str):
                        msg = "Argument 'sql' must be a string"
                        raise LDRdataFindClientException(msg)

                if not isinstance(rexp, str):
                        msg = "Argument 'rexp' must be a string"
                        raise LDRdataFindClientException(msg)

                if offset:
                    if not isinstance(offset, int):
                            msg = "Argument 'offset' must be a positive integer or zero"
                            raise LDRdataFindClientException(msg)
                
                    if offset < 0:
                            msg = "Argument 'offset' must be a positive integer or zero"
                            raise LDRdataFindClientException(msg)
                
                if number:
                    if not isinstance(number, int):
                            msg = "Argument 'number' must be a positive integer"
                            raise LDRdataFindClientException(msg)
                
                    if number <= 0:
                            msg = "Argument 'number' must be a positive integer"
                            raise LDRdataFindClientException(msg)

                msg = "METAREPFN\0%s\0%s\0" % (str(sql), str(rexp))
                if offset:
                        msg += "%d\0" % offset
                else:
                        msg += "0\0"
                if number:
                        msg += "%d\0" % number
                self.__connect__()
                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for PFNs with metadata query %s : %s" % (sql, str(output[0]))
                        raise LDRdataFindClientException(msg)

                

                return output

        def pfnQueryWithMetadataUnion(self, sql1, sql2):
                """
                Query LDRdataFindServer to find the PFNs(s) for LFN(s) with the appropriate 
                metadata values found by the union of two queries.

                @param sql1: clause that will be part of a SQL query done to find the LFN(s)
                        with the appropriate metadata values. A an example would be
                        
                        gpsStart >= '777777777' AND gpsStart <= '888888888' AND instrument = 'H' AND runTag = 'S2'

                @param sql2: second clause with similar form as above. An example would be

                        gpsEnd >= '777777777' AND gpsEnd <= '888888888' AND instrument = 'H' AND runTag = 'S2'

                @return: list of strings representing the PFN(s) found
                """

                # check arguments
                if not isinstance(sql1, str):
                        msg = "Argument 'sql1' must be a string"
                        raise LDRdataFindClientException(msg)

                if not isinstance(sql2, str):
                        msg = "Argument 'sql2' must be a string"
                        raise LDRdataFindClientException(msg)

                msg = "METAPFNUNION\0%s\0%s\0" % (str(sql1), str(sql2))
                self.__connect__()
                self.sfile.write(msg)

                ret, output = self.__response__()

                if ret:
                        msg = "Error querying LDRdataFindServer for PFNs with metadata queryies %s AND %s : %s" % (sql1, sql2, str(output[0]))
                        raise LDRdataFindClientException(msg)

                

                return output


class LSCdataFindClient(LDRdataFindClient):
        """
        Class that represents this client interacting with a LDRdataFindServer in
        order to find LSC data.
        """
        def __init__(self, host, port=30010, noproxy=False, disablehostauth=False ):
                """
                Open a connection to a LDRdataFindServer and return an instance of
                class LDRdataFindClient. One of the public methods can then be 
                called to send a request to the server.

                @param host: the host on which the LDRdataFindServer runs
                @type host: string

                @param port: port on which the LDRdataFindServer listens
                @type port: integer
                
                @param noproxy: whether or not to attempt to use a proxy
                                certificate.
                @type noproxy: boolean

                @param disablehostauth: whether or not to attempt to use GSI
                                        host authentication.
                @type disablehostauth: boolean

                @return: Instance of LSCdataFindClient
                """
                # check arguments
                if not isinstance(host, str):
                        msg = "Argument 'host' must be a string"
                        raise LSCdataFindClientException(msg)

                if not isinstance(port, str):
                        msg = "Argument 'port' must be a positive integer"
                        raise LSCdataFindClientException(msg)

                if port <= 0:
                        msg = "Argument 'port' must be a positive integer"
                        raise LSCdataFindClientException(msg)
                 
                LDRdataFindClient.__init__(self, host, port, noproxy, disablehostauth)

        def __check_gps(self, gpsString):
                """
                Minimal checking on GPS time strings. Raises a LSCdataFindClientException if
                the GPS time string is not 9 digits long.

                @param gpsString: The string representing the 9 digit GPS time.

                @returns: None
                """
                if len(gpsString) != 9:
                        msg = "GPS times must be 9 digits"
                        raise LSCdataFindClientException(msg)

                try:
                        a = int(gpsString)
                except Exception, e:
                        msg = "GPS times must be 9 digits"
                        raise LSCdataFindClientException(msg)


        def ping(self, argDict):
                """
                Ping the LDRdataFindServer and return any response sent back.

                @param argDict: Dictionary of arguments passed to all methods.

                @return: None
                """
                response = LDRdataFindClient.ping(self)
                return response

        
        def showObservatories(self, argDict):
                """
                Query LDRdataFindServer for the distinct values for the 'instrument' attribute
                in the metadata table.

                @param argDict: Dictionary of arguments passed to all methods.

                @return: None
                """

                distinctValueList = LDRdataFindClient.distinctAttrValues(self, "site")
                return distinctValueList
        
        
        def showTimes(self,argDict):
                 """
                 Query LDRdataFind server for gps times for existing data corresponding
                 to the specified frame type.
                 
                 @param argDict: Dictionary of arguments passed to all methods.
                 
                 @return: results of query
                 """
                 start = argDict['start']
                 end   = argDict['end']
                 strict = argDict['strict']
                 site = argDict['observatory']
                 mytype = argDict['type']
                 
                 # Perform basic argument check
                 if start:
                     self.__check_gps(start)
                 if end:
                     self.__check_gps(end)
                 if not site:
                     msg = "An observatory must be specified."
                     raise LDRdataFindClientException(msg)
                 if not mytype:
                     msg = "A type must be specified."
                     raise LDRdataFindClientException(msg)
                 
                 timelist = LDRdataFindClient.timeQuery(self,mytype,site,start,end,strict)
                 return timelist
                 
        
        def showTypes(self, argDict):
                """
                Query LDRdataFindServer for the distinct values for the 'frameType' attribute
                in the metadata table.

                @param argDict: Dictionary of arguments passed to all methods.

                @return: None
                """
                distinctValueList = LDRdataFindClient.distinctAttrValues(self, "frameType")
                return distinctValueList


        def singleFrameFind(self, argDict):
                """
                Query the LDRdataFindServer for URLs for a given file.
                
                @param argDict: Dictionary of arguments passed to all methods.

                @return: None
                """

                lfn = argDict['filename']
                pfnList = LDRdataFindClient.pfnQuery(self, lfn)
                return pfnlist(pfnList)


        def findFrameNames(self, argDict):
                """
                Query the LDRdataFindServer for frame files from a particular observatory,
                with a particular frame type, for a particular range of GPS times.

                
                @param argDict: Dictionary of arguments passed to all methods.

                @return: None
                """
                instrument = argDict['observatory']
                type = argDict['type']
                start = argDict['start']
                end = argDict['end']
                offset = argDict['offset']
                number = argDict['limit']
                strict = argDict['strict']

                
               
                
                # check that combination of command-line arguments is sound
                if (not instrument) or (not type) or (not start) or (not end):
                        msg = """\
Bad combination of command line arguments:
--observatory --type --gps-start-time --gps-end-time must all
be present when searching for groups of files
"""
                        raise LSCdataFindClientException(msg)
                
                self.__check_gps(start)
                self.__check_gps(end)
                
                if int(end) < int(start):
                        msg = "Supplied start time, %s, is later than supplied end time, %s" % (start,end)
                        raise LSCdataFindClientException(msg)

                q1 = LDRMetadataQuery()
                q2 = LDRMetadataQuery()

                if strict:
                        q1.set_query("(gpsStart >= '%s' AND gpsStart < '%s') AND gpsEnd < '%s' AND instrument = '%s' AND frameType = '%s'" % (start, end, end, instrument, type))
                        q1.set_sort_attribute("gpsStart")
                        q1.set_sort_order("ASC")

                        if offset: q1.set_offset(offset)
                        if number: q1.set_limit(number)

                        lfnList = LDRdataFindClient.lfnQueryWithMetadata(self, [q1])
                else:
                        q1.set_query("gpsStart >= '%s' AND gpsStart < '%s' AND instrument = '%s' AND frameType = '%s'" % (start, end, instrument, type))
                        q1.set_sort_attribute("gpsStart")
                        q1.set_sort_order("ASC")

                        q2.set_query("gpsEnd > '%s' AND gpsEnd <= '%s' AND instrument = '%s' AND frameType = '%s'" % (start, end, instrument, type))
                        q2.set_sort_attribute("gpsStart")
                        q2.set_sort_order("ASC")

                        if offset:
                                q1.set_offset(offset)
                                q2.set_offset(offset)

                        if number:
                                q1.set_limit(number)
                                q2.set_limit(number)
                                
                        lfnList = LDRdataFindClient.lfnQueryWithMetadata(self, [q1, q2])

                return lfnlist(lfnList)
                

        def findFrameURLs(self, argDict):
                """
                Query the LDRdataFindServer for the URLs for frame files from a particular
                observatory, with a particular frame type, for a particular range of GPS times.
                """
                instrument = argDict['observatory']
                type = argDict['type']
                start = argDict['start']
                end = argDict['end']
                offset = argDict['offset']
                number = argDict['limit']
                strict = argDict['strict']

                # check that combintation of command-line arguments is sound
                if (not instrument) or (not type) or (not start) or (not end):
                        msg = """\
Bad combination of command line arguments:
--observatory --type --gps-start-time --gps-end-time must all
be present when searching for groups of files
"""
                        raise LSCdataFindClientException(msg)
                
                self.__check_gps(start)
                self.__check_gps(end)
                
                if int(end) < int(start):
                        msg = "Supplied start time, %s, is later than supplied end time, %s" % (start,end)
                        raise LSCdataFindClientException(msg)
                q1 = LDRMetadataQuery()
                q2 = LDRMetadataQuery()

                if strict:
                        q1.set_query("(gpsStart >= '%s' AND gpsStart < '%s') AND gpsEnd < '%s' AND instrument = '%s' AND frameType = '%s'" % (start, end, end, instrument, type))
                        q1.set_sort_attribute("gpsStart")
                        q1.set_sort_order("ASC")

                        if offset: q1.set_offset(offset)
                        if number: q1.set_limit(number)

                        pfnList = LDRdataFindClient.pfnQueryWithMetadata(self, [q1])
                else:
                        q1.set_query("gpsStart >= '%s' AND gpsStart < '%s' AND instrument = '%s' AND frameType = '%s'" % (start, end, instrument, type))
                        q1.set_sort_attribute("gpsStart")
                        q1.set_sort_order("ASC")

                        q2.set_query("gpsEnd > '%s' AND gpsEnd <= '%s' AND instrument = '%s' AND frameType = '%s'" % (start, end, instrument, type))
                        q2.set_sort_attribute("gpsStart")
                        q2.set_sort_order("ASC")

                        if offset:
                                q1.set_offset(offset)
                                q2.set_offset(offset)

                        if number:
                                q1.set_limit(number)
                                q2.set_limit(number)
                                
                        pfnList = LDRdataFindClient.pfnQueryWithMetadata(self, [q1, q2])

                return pfnlist(pfnList)


        def findFrameURLsFilter(self, argDict):
                """
                Query the LDRdataFindServer for the URLs for frame files from a particular
                observatory, with a particular frame type, for a particular range of GPS times,
                and filter the results using either the type of URL or by matching against a
                regular expression.
                """
                instrument = argDict['observatory']
                type = argDict['type']
                start = argDict['start']
                end = argDict['end']
                match = argDict['match']
                urlType = argDict['urlType']
                offset = argDict['offset']
                number = argDict['limit']
                strict = argDict['strict']

                if (not instrument) or (not type) or (not start) or (not end):
                        msg = """\
Bad combination of command line arguments:
--observatory --type --gps-start-time --gps-end-time must all
be present when searching for groups of files
"""
                        raise LSCdataFindClientException(msg)

                if offset and not number:
                        msg = "--limit must be used if --offset is used"
                        raise LSCdataFindClientException(msg)

                self.__check_gps(start)
                self.__check_gps(end)

                if int(end) < int(start):
                        msg = "Supplied start time, %s, is later than supplied end time, %s" % (start,end)
                        raise LSCdataFindClientException(msg)


                # should do sanity check here on the urlType and match that have been passed
                # and check for proper quoting
                if match and urlType:
                        rexp = "^%s.*%s" % (str(urlType), str(match))
                elif match and not urlType:
                        rexp = "%s" % str(match)
                elif not match and urlType:
                        rexp = "^%s" % str(urlType)

                q1 = LDRMetadataQuery()
                q2 = LDRMetadataQuery()

                if strict:
                        q1.set_query("(gpsStart >= '%s' AND gpsStart < '%s') AND gpsEnd < '%s' AND instrument = '%s' AND frameType = '%s'" % (start, end, end, instrument, type))
                        q1.set_sort_attribute("gpsStart")
                        q1.set_sort_order("ASC")
                        q1.set_regex_filter(rexp)

                        if offset: q1.set_offset(offset)
                        if number: q1.set_limit(number)

                        pfnList = LDRdataFindClient.pfnQueryWithMetadata(self, [q1])
                else:
                        q1.set_query("gpsStart >= '%s' AND gpsStart < '%s' AND instrument = '%s' AND frameType = '%s'" % (start, end, instrument, type))
                        q1.set_sort_attribute("gpsStart")
                        q1.set_sort_order("ASC")
                        q1.set_regex_filter(rexp)

                        q2.set_query("gpsEnd > '%s' AND gpsEnd <= '%s' AND instrument = '%s' AND frameType = '%s'" % (start, end, instrument, type))
                        q2.set_sort_attribute("gpsStart")
                        q2.set_sort_order("ASC")
                        q2.set_regex_filter(rexp)

                        if offset:
                                q1.set_offset(offset)
                                q2.set_offset(offset)

                        if number:
                                q1.set_limit(number)
                                q2.set_limit(number)
                                
                        pfnList = LDRdataFindClient.pfnQueryWithMetadata(self, [q1, q2])

                return pfnlist(pfnList)
