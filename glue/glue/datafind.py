# -*- coding: utf-8 -*-
# Copyright (C) 2012  Scott Koranda, Duncan Macleod
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

"""The client library for the LIGO Data Replicator (LDR) service.

The DataFind service allows users to query for the location of
Gravitational-Wave Frame (GWF) files containing data from the current
LIGO and Virgo gravitational-wave detectors.

This module defines the L{GWDataFindHTTPConnection} and
L{GWDataFindHTTPSConnection} class objects, for connecting to an LDR
server in open and authenticated access modes respectively.
The authenticated L{GWDataFindHTTPSConnection} connection requires users
have a valid X509 certificate that is registered with the server in
question.

A new connection can be opened as follows:

>>> from glue.datafind import GWDataFindHTTPConnection
>>> connection = GWDataFindHTTPConnection(host, port)

and similar for the HTTPS version.

Users on the LIGO Data Grid (LDG) can connect without giving the name of
the relevant host, so long as the C{LIGO_DATAFIND_SERVER} environment
variable is defined:

>>> connection = GWDataFindHTTPConnection()

Queries for frames can be made using the L{find_frame_urls<GWDataFindHTTPConnection.find_frame_urls>} method of the
relevant connection:

>>> cache = connection.find_frame_urls('L', 'L1_R', 1093564816, 1093651216)

By default, the returned L{Cache<glue.lal.Cache>} object includes both C{gsiftp} and local C{file} versions of each frame, but the C{urlfile} keyword argument can be given to return only one of those:

>>> cache = connection.find_frame_urls('L', 'L1_R', 1093564816, 1093651216, urltype='file')

See the documentation for each connection method for more detailed examples.
"""

from __future__ import division

import os
import sys
import time
import calendar
import httplib
import M2Crypto
import re
import unittest

from cjson import decode
from glue import (lal, git_version, segments)

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__credits__ = "Scott Koranda <scott.koranda@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date

_server_env = "LIGO_DATAFIND_SERVER"
_url_prefix = "/LDR/services/data/v1"


class GWDataFindHTTPConnection(httplib.HTTPConnection):
    """Connection to LIGO data replicator service using HTTP.

    @param host: the name of the server with which to connect
    @param port: the port on which to connect
    @param **kwargs: other keyword arguments accepted by
                     L{httplib.HTTPConnection}

    @type  host: L{str}
    @type  port: L{int}
    """
    LIGOTimeGPSType = lal.LIGOTimeGPS
    def __init__(self, host=None, **kwargs):
        """Connect to the LDR host using HTTPS. Default host is
        defined by the %s environment variable.
        """
        if not host:
            host,port = find_server()
            kwargs.setdefault("port", port)
        httplib.HTTPConnection.__init__(self, host, **kwargs)
    __init__.__doc__ %= _server_env

    def _requestresponse(self, method, url, body=None, headers={}):
        """Internal method to perform request and verify reponse.

        @param method: name of the method to use (e.g. 'GET')
        @param url   : remote URL to query

        @type  method: L{str}
        @type  url   : L{str}

        @returns: L{str} response from server query

        @raises RuntimeError: if query is unsuccessful
        """
        try:
            self.request(method, url)
            response = self.getresponse()
        except Exception as e:
            raise RuntimeError("Unable to query server %s: %s\n\n"
                               "Perhaps you need a valid proxy credential?\n"
                               % (self.host, e))
        if response.status != 200:
            raise RuntimeError("Server returned code %d: %s%s"
                               % (response.status, response.reason,
                                  response.read()))
        return response

    def ping(self):
        """Ping the LDR host to test for life

        @raises RuntimeError: when ping fails
        @returns: 0 if ping was successful
        """
        url = "%s/gwf/%s/%s/%s,%s" % (_url_prefix, 'H', 'R', '1', '2')
        self._requestresponse("HEAD", url)
        return 0

    def find_observatories(self, match=None):
        """Query the LDR host for observatories. Use match to
        restrict returned observatories to those matching the
        regular expression.

        Example:

        >>> connection.find_observatories()
        ['AGHLT', 'G', 'GHLTV', 'GHLV', 'GHT', 'H', 'HL', 'HLT',
         'L', 'T', 'V', 'Z']
        >>> connection.find_observatories("H")
        ['H', 'HL', 'HLT']

        @type  match: L{str}
        @param match:
            name to match return observatories against

        @returns: L{list} of observatory prefixes
        """
        url = "%s/gwf.json" % _url_prefix
        response = self._requestresponse("GET", url)
        sitelist = sorted(set(decode(response.read())))
        if match:
            regmatch = re.compile(match)
            sitelist = [site for site in sitelist if regmatch.search(site)]
        return sitelist

    def find_types(self, site=None, match=None):
        """Query the LDR host for frame types. Use site to restrict
        query to given observatory prefix, and use match to restrict
        returned types to those matching the regular expression.

        Example:

        >>> connection.find_types("L", "RDS")
        ['L1_RDS_C01_LX',
         'L1_RDS_C02_LX',
         'L1_RDS_C03_L2',
         'L1_RDS_R_L1',
         'L1_RDS_R_L3',
         'L1_RDS_R_L4',
         'PEM_RDS_A6',
         'RDS_R_L1',
         'RDS_R_L2',
         'RDS_R_L3',
         'TESTPEM_RDS_A6']

        @param  site: single-character name of site to match
        @param match: type-name to match against

        @type  site: L{str}
        @type match: L{str}

        @returns: L{list} of frame types
        """
        if site:
            url = "%s/gwf/%s.json" % (_url_prefix, site[0])
        else:
            url = "%s/gwf/all.json" % _url_prefix
        response = self._requestresponse("GET", url)
        typelist = sorted(set(decode(response.read())))
        if match:
            regmatch = re.compile(match)
            typelist = [type for type in typelist if regmatch.search(type)]
        return typelist

    def find_times(self, site, frametype, gpsstart=None, gpsend=None):
        """Query the LDR for times for which frames are avaliable

        Use gpsstart and gpsend to restrict the returned times to
        this semiopen interval.

        @returns: L{segmentlist<glue.segments.segmentlist>}

        @param site:
            single-character name of site to match
        @param frametype:
            name of frametype to match
        @param gpsstart:
            integer GPS start time of query
        @param gpsend:
            integer GPS end time of query

        @type       site: L{str}
        @type  frametype: L{str}
        @type   gpsstart: L{int}
        @type     gpsend: L{int}
        """
        if gpsstart and gpsend:
            url = ("%s/gwf/%s/%s/segments/%d,%d.json"
                   % (_url_prefix, site, frametype, gpsstart, gpsend))
        else:
            url = ("%s/gwf/%s/%s/segments.json"
                   % (_url_prefix, site, frametype))

        response = self._requestresponse("GET", url)
        segmentlist = decode(response.read())
        return segments.segmentlist(map(segments.segment, segmentlist))

    def find_frame(self, framefile, urltype=None, on_missing="warn"):
        """Query the LDR host for a single framefile

        @returns: L{Cache<glue.lal.Cache>}

        @param frametype:
            name of frametype to match
        @param urltype:
            file scheme to search for (e.g. 'file')
        @param on_missing:
            what to do when the requested frame isn't found, one of:
                - C{'warn'} (default): print a warning,
                - C{'error'}: raise an L{RuntimeError}, or
                - C{'ignore'}: do nothing

        @type  frametype: L{str}
        @type    urltype: L{str}
        @type on_missing: L{str}

        @raises RuntimeError: if given framefile is malformed
        """
        if on_missing not in ("warn", "error", "ignore"):
            raise ValueError("on_missing must be 'warn', 'error', or 'ignore'.")
        framefile = os.path.basename(framefile)
        # parse file name for site, frame type
        try:
            site,frametype,_,_ = framefile.split("-")
        except Exception as e:
            raise RuntimeError("Error parsing filename %s: %s" % (framefile, e))
        url = ("%s/gwf/%s/%s/%s.json"
              % (_url_prefix, site, frametype, framefile))
        response = self._requestresponse("GET", url)
        urllist  = decode(response.read())
        if len(urllist) == 0:
            if on_missing == "warn":
                sys.stderr.write("No files found!\n")
            elif on_missing == "error":
                raise RuntimeError("No files found!")
        # verify urltype is what we want
        cache = lal.Cache(e for e in
                 [lal.CacheEntry.from_T050017(x, coltype=self.LIGOTimeGPSType)
                  for x in urllist] if not urltype or e.scheme == urltype)
        return cache

    def find_latest(self, site, frametype, urltype=None, on_missing="warn"):
        """Query for the most recent framefile of a given type.

        @param  site:
            single-character name of site to match
        @param frametype:
            name of frametype to match
        @param urltype:
            file scheme to search for (e.g. 'file')
        @param on_missing:
            what to do when the requested frame isn't found, one of:
                - C{'warn'} (default): print a warning,
                - C{'error'}: raise an L{RuntimeError}, or
                - C{'ignore'}: do nothing

        @type       site: L{str}
        @type  frametype: L{str}
        @type    urltype: L{str}
        @type on_missing: L{str}

        @returns: L{Cache<glue.lal.Cache>} with one
                  L{entry<glue.lal.CacheEntry>}

        @raises RuntimeError: if given framefile is malformed
        @raises RuntimeError: if no frames are found and C{on_missing='error'}
        """
        if on_missing not in ('warn', 'error', 'ignore'):
            raise ValueError("on_missing must be 'warn', 'error', or 'ignore'.")
        url = "%s/gwf/%s/%s/latest" % (_url_prefix, site, frametype)
        # if a URL type is specified append it to the path
        if urltype:
            url += "/%s" % urltype
        # request JSON output
        url += ".json"
        response = self._requestresponse("GET", url)
        urllist  = decode(response.read())
        if len(urllist) == 0:
            if on_missing == "warn":
                sys.stderr.write("No files found!\n")
            elif on_missing == "error":
                raise RuntimeError("No files found!")
        return lal.Cache([lal.CacheEntry.from_T050017(x,
                         coltype=self.LIGOTimeGPSType) for x in urllist])

    def find_frame_urls(self, site, frametype, gpsstart, gpsend,
                        match=None, urltype=None, on_gaps="warn"):
        """Find the framefiles for the given type in the [start, end) interval
        frame

        @param site:
            single-character name of site to match
        @param frametype:
            name of frametype to match
        @param gpsstart:
            integer GPS start time of query
        @param gpsend:
            integer GPS end time of query
        @param match:
            regular expression to match against
        @param urltype:
            file scheme to search for (e.g. 'file')
        @param on_gaps:
            what to do when the requested frame isn't found, one of:
                - C{'warn'} (default): print a warning,
                - C{'error'}: raise an L{RuntimeError}, or
                - C{'ignore'}: do nothing

        @type       site: L{str}
        @type  frametype: L{str}
        @type   gpsstart: L{int}
        @type     gpsend: L{int}
        @type      match: L{str}
        @type    urltype: L{str}
        @type    on_gaps: L{str}

        @returns: L{Cache<glue.lal.Cache>}

        @raises RuntimeError: if gaps are found and C{on_gaps='error'}
        """
        if on_gaps not in ("warn", "error", "ignore"):
            raise ValueError("on_gaps must be 'warn', 'error', or 'ignore'.")
        url = ("%s/gwf/%s/%s/%s,%s"
               % (_url_prefix, site, frametype, gpsstart, gpsend))
        # if a URL type is specified append it to the path
        if urltype:
            url += "/%s" % urltype
        # request JSON output
        url += ".json"
        # append a regex if input
        if match:
            url += "?match=%s" % match
        # make query
        response = self._requestresponse("GET", url)
        urllist  = decode(response.read())

        out = lal.Cache([lal.CacheEntry.from_T050017(x,
                         coltype=self.LIGOTimeGPSType) for x in urllist])

        if on_gaps == "ignore":
            return out
        else:
            span    = segments.segment(gpsstart, gpsend)
            seglist = segments.segmentlist(e.segment for e in out).coalesce()
            missing = (segments.segmentlist([span]) - seglist).coalesce()
            if span in seglist:
                return out
            else:
                msg = "Missing segments: \n%s" % "\n".join(map(str, missing))
                if on_gaps=="warn":
                    sys.stderr.write("%s\n" % msg)
                    return out
                else:
                    raise RuntimeError(msg)

class GWDataFindHTTPSConnection(httplib.HTTPSConnection, GWDataFindHTTPConnection):
    """Secured connection to LIGO data replicator service using HTTPS.
    """
    def __init__(self, host=None, **kwargs):
        """Connect to the LDR host using HTTPS.

        Default host is defined by the %s environment variable.
        """
        if not host:
            host, port = find_server()
            kwargs.setdefault("port", port)
        httplib.HTTPSConnection.__init__(self, host, **kwargs)
    __init__.__doc__ %= _server_env


def validate_proxy(path):
    """Validate the users X509 proxy certificate

    Tests that the proxy certificate is RFC 3820 compliant and that it
    is valid for at least the next 15 minutes.

    @returns: L{True} if the certificate validates
    @raises RuntimeError: if the certificate cannot be validated
    """
    # load the proxy from path
    try:
        proxy = M2Crypto.X509.load_cert(path)
    except Exception as e:
        msg = "Unable to load proxy from path %s : %s" % (path, e)
        raise RuntimeError(msg)

    # make sure the proxy is RFC 3820 compliant
    try:
        proxy.get_ext("proxyCertInfo")
    except LookupError:
        subject = proxy.get_subject().as_text()
        if re.search(r'.+CN=proxy$', subject):
            raise RuntimeError("Could not find a RFC 3820 compliant proxy "
                               "credential. Please run 'grid-proxy-init -rfc' "
                               "and try again.")

    # attempt to make sure the proxy is still good for more than 15 minutes
    try:
        expireASN1 = proxy.get_not_after().__str__()
        expireGMT  = time.strptime(expireASN1, "%b %d %H:%M:%S %Y %Z")
        expireUTC  = calendar.timegm(expireGMT)
        now = int(time.time())
        secondsLeft = expireUTC - now
    except Exception as e:
        # problem getting or parsing time so just let the client
        # continue and pass the issue along to the server
        secondsLeft = 3600

    if secondsLeft <= 0:
        raise RuntimeError("Your proxy certificate is expired.\n"
                           "Please generate a new proxy certificate and "
                           "try again. ")
    if secondsLeft < (60 * 15):
        raise RuntimeError("Your proxy certificate expires in less than 15 "
                           "minutes.\nPlease generate a new proxy "
                           "certificate and try again.")

    # return True to indicate validated proxy
    return True

def find_credential():
    """Locate the users X509 certificate and key files

    This method uses the C{X509_USER_CERT} and C{X509_USER_KEY} to locate
    valid proxy information. If those are not found, the standard location
    in /tmp/ is searched.

    @raises RuntimeError: if the proxy found via either method cannot
                          be validated
    @raises RuntimeError: if the cert and key files cannot be located
    """

    rfc_proxy_msg = ("Could not find a RFC 3820 compliant proxy credential."
                     "Please run 'grid-proxy-init -rfc' and try again.")

    # use X509_USER_PROXY from environment if set
    if os.environ.has_key('X509_USER_PROXY'):
        filePath = os.environ['X509_USER_PROXY']
        if validate_proxy(filePath):
            return filePath, filePath
        else:
            raise RuntimeError(rfc_proxy_msg)

    # use X509_USER_CERT and X509_USER_KEY if set
    if (os.environ.has_key('X509_USER_CERT') and
        os.environ.has_key('X509_USER_KEY')):
        certFile = os.environ['X509_USER_CERT']
        keyFile = os.environ['X509_USER_KEY']
        return certFile, keyFile

    # search for proxy file on disk
    uid = os.getuid()
    path = "/tmp/x509up_u%d" % uid

    if os.access(path, os.R_OK):
        if validate_proxy(path):
            return path, path
        else:
            raise RuntimeError(rfc_proxy_msg)

    # if we get here could not find a credential
    raise RuntimeError(rfc_proxy_msg)

def find_server():
    """Find the default server host from the environment

    This method uses the C{LIGO_DATAFIND_SERVER} variable to construct
    a C{(host, port)} tuple.

    @returns: C{(host, port)}: the L{str} host name and L{int} port number

    @raises RuntimeError: if the C{LIGO_DATAFIND_SERVER} environment variable
                          is not set
    """

    if os.environ.has_key(_server_env):
        host = os.environ[_server_env]
        port = None
        if re.search(':', host):
            host, port = host.split(':', 1)
            if port:
                port = int(port)
        return host, port
    else:
        raise RuntimeError("Environment variable %s is not set" % _server_env)


class TestLDR(unittest.TestCase):
    """Small suite of test functions.

    Probably won't work if you're not on an LDAS
    machine...
    """
    def test_HTTPConnection(self):
        h = GWDataFindHTTPConnection()
        h.close()

    def test_HTTPSConnection(self):
        h = GWDataFindHTTPSConnection()
        h.close()

    def test_ping(self):
        h = GWDataFindHTTPConnection()
        h.ping()
        h.close()

    def test_latest(self):
        h = GWDataFindHTTPConnection()
        h.find_latest("L", "R")
        h.close()

    def test_find_observatories(self):
        h = GWDataFindHTTPConnection()
        h.find_observatories()
        h.close()

    def test_find_times(self):
        h = GWDataFindHTTPConnection()
        h.find_times("L", "R")
        h.close()

    def test_find_frame_urls(self):
        h = GWDataFindHTTPConnection()
        h.find_frame_urls("L", "R", 1000000000, 1000001000, on_gaps="ignore")
        h.close()

if __name__ == "__main__":
    unittest.main()
