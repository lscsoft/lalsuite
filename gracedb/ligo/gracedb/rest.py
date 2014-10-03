
# Copyright (C) 2012  LIGO Scientific Collaboration
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


import httplib, socket
import mimetypes
import urllib
import os, sys
import json
from urlparse import urlparse

DEFAULT_SERVICE_URL = "https://gracedb.ligo.org/api/"

#-----------------------------------------------------------------
# Exception(s)

class HTTPError(Exception):
    def __init__(self, status, reason, message):
        self.status = status
        self.reason = reason
        self.message = message
        Exception.__init__(self, (status, reason+" / "+message))

#-----------------------------------------------------------------
# HTTP/S Proxy classses
# Taken from: http://code.activestate.com/recipes/456195/

class ProxyHTTPConnection(httplib.HTTPConnection):

    _ports = {'http' : 80, 'https' : 443}

    def request(self, method, url, body=None, headers={}):
        #request is called before connect, so can interpret url and get
        #real host/port to be used to make CONNECT request to proxy
        o = urlparse(url)
        proto = o.scheme
        port = o.port
        host = o.hostname
        if proto is None:
            raise ValueError, "unknown URL type: %s" % url
        if port is None:
            try:
                port = self._ports[proto]
            except KeyError:
                raise ValueError, "unknown protocol for: %s" % url
        self._real_host = host
        self._real_port = port
        httplib.HTTPConnection.request(self, method, url, body, headers)


    def connect(self):
        httplib.HTTPConnection.connect(self)
        #send proxy CONNECT request
        self.send("CONNECT %s:%d HTTP/1.0\r\n\r\n" % (self._real_host, self._real_port))
        #expect a HTTP/1.0 200 Connection established
        response = self.response_class(self.sock, strict=self.strict, method=self._method)
        (version, code, message) = response._read_status()
        #probably here we can handle auth requests...
        if code != 200:
            #proxy returned and error, abort connection, and raise exception
            self.close()
            raise socket.error, "Proxy connection failed: %d %s" % (code, message.strip())
        #eat up header block from proxy....
        while True:
            #should not use directly fp probablu
            line = response.fp.readline()
            if line == '\r\n': break

class ProxyHTTPSConnection(ProxyHTTPConnection):

    default_port = 443

    def __init__(self, host, port = None, key_file = None, cert_file = None, strict = None):
        ProxyHTTPConnection.__init__(self, host, port)
        self.key_file = key_file
        self.cert_file = cert_file

    def connect(self):
        ProxyHTTPConnection.connect(self)
        #make the sock ssl-aware
        ssl = socket.ssl(self.sock, self.key_file, self.cert_file)
        self.sock = httplib.FakeSocket(self.sock, ssl)

#-----------------------------------------------------------------
# Generic GSI REST

class GsiRest(object):
    """docstring for GracedbRest"""
    def __init__(self, 
            url=DEFAULT_SERVICE_URL, 
            proxy_host=None, 
            proxy_port=3128, 
            cred=None):
        if not cred:
            cred = findUserCredentials()
        if isinstance(cred, (list, tuple)):
            self.cert, self.key = cred
        elif cred:
            self.cert = self.key = cred
        self.connection = None

        if proxy_host:
            self.connector = lambda: ProxyHTTPSConnection(
                    proxy_host, proxy_port, key_file=self.key, cert_file=self.cert)
        else:
            o = urlparse(url)
            port = o.port
            host = o.hostname
            port = port or 443
            self.connector = lambda: httplib.HTTPSConnection(
                    host, port, key_file=self.key, cert_file=self.cert)

    def getConnection(self):
        return self.connector()

    def request(self, method, url, body=None, headers=None, priming_url=None):
        # Bug in Python (versions < 2.7.1 (?))
        # http://bugs.python.org/issue11898
        # if the URL is unicode and the body of a request is binary,
        # the POST/PUT action fails because it tries to concatenate
        # the two which fails due to encoding problems.
        # Workaround is to cast all URLs to str.
        # This is probably bad in general,
        # but for our purposes, today, this will do.
        url = url and str(url)
        priming_url = priming_url and str(priming_url)

        conn = self.getConnection()
        if priming_url:
            conn.request("GET", priming_url, headers={'connection' : 'keep-alive'})
            response = conn.getresponse()
            if response.status != 200:
                response = self.adjustResponse(response)
            else:
                # Throw away the response and make sure to read the body.
                response = response.read()
        conn.request(method, url, body, headers or {})
        response = conn.getresponse()
        return self.adjustResponse(response)

    def adjustResponse(self, response):
#       XXX WRONG.
        if response.status >= 400:
            raise HTTPError(response.status, response.reason, response.read())
        response.json = lambda: json.loads(response.read())
        return response

    def get(self, url, headers=None):
        return self.request("GET", url, headers=headers)

    def head(self, url, headers=None):
        return self.request("HEAD", url, headers=headers)

    def delete(self, url, headers=None):
        return self.request("DELETE", url, headers=headers)

    def options(self, url, headers=None):
        return self.request("OPTIONS", url, headers=headers)

    def post(self, *args, **kwargs):
        return self.post_or_put("POST", *args, **kwargs)

    def put(self, *args, **kwargs):
        return self.post_or_put("PUT", *args, **kwargs)

    def post_or_put(self, method, url, body=None, headers=None, files=None):
        headers = headers or {}
        if not files:
            # Simple urlencoded body
            if isinstance(body, dict):
#           XXX What about the headers in the params?
                if 'content-type' not in headers:
                    headers['content-type'] = "application/json"
                body = json.dumps(body)
        else:
            body = body or {}
            if isinstance(body, dict):
                body = body.items()
            content_type, body = encode_multipart_formdata(body, files)
#           XXX What about the headers in the params?
            headers = {
                'content-type': content_type,
                'content-length': str(len(body)),
#               'connection': 'keep-alive',
            }
        return self.request(method, url, body, headers)

#------------------------------------------------------------------
# GraceDB
#
# Example Gracedb REST client.

class GraceDb(GsiRest):
    """docstring for GraceDb"""
    def __init__(self, service_url=DEFAULT_SERVICE_URL, 
            proxy_host=None, proxy_port=3128, *args, **kwargs):
        GsiRest.__init__(self, service_url, proxy_host, proxy_port, *args, **kwargs)

        self.service_url = service_url
        self._service_info = None

    @property
    def service_info(self):
        if not self._service_info:
            self._service_info = self.request("GET", self.service_url).json()
        return self._service_info

    @property
    def links(self):
        return self.service_info.get('links')

    @property
    def templates(self):
        return self.service_info.get('templates')

    @property
    def groups(self):
        return self.service_info.get('groups')
    
    @property
    def pipelines(self):
        return self.service_info.get('pipelines')

    @property
    def searches(self):
        return self.service_info.get('searches')

    @property
    def em_groups(self):
        return self.service_info.get('em-groups')

    @property
    def wavebands(self):
        return self.service_info.get('wavebands')

    @property
    def eel_statuses(self):
        return self.service_info.get('eel-statuses')

    @property
    def obs_statuses(self):
        return self.service_info.get('obs-statuses')

    def request(self, method, *args, **kwargs):
        if method.lower() in ['post', 'put']:
            kwargs['priming_url'] = self.service_url
        return GsiRest.request(self, method, *args, **kwargs)

    def _getCode(self, input_value, code_dict):
        """Check if input is valid.
           Return coded version if it is"""
        #  code_dict is dict of { code : descriptive_name }
        if input_value in code_dict.keys():
            # Already coded
            return input_value
        if not input_value in code_dict.values():
            # Not valid
            return None
        return [code
                for code, name in code_dict.items()
                if name == input_value][0]

    # Search and filecontents are optional when creating an event.
    def createEvent(self, group, pipeline, filename, search=None, filecontents=None):
        errors = []
        if group not in self.groups:
            errors += ["bad group"]
        if pipeline not in self.pipelines:
            errors += ["bad pipeline"]
        if search and search not in self.searches:
            errors += ["bad search"]
        if errors:
            # XXX Terrible error messages / weak exception type
            raise Exception(str(errors))
        if filecontents is None:
            if filename == '-':
                filename = 'initial.data'
                filecontents = sys.stdin.read()
            else:
                filecontents = open(filename, 'rb').read() 
        fields = [
                  ('group', group),
                  ('pipeline', pipeline),
                 ]
        if search:
            fields.append(('search', search))
        files = [('eventFile', filename, filecontents)]
        # Python httplib bug?  unicode link
        uri = str(self.links['events'])
        return self.post(uri, fields, files=files)

    def replaceEvent(self, graceid, filename, filecontents=None):
        if filecontents is None:
            # Note: not allowing filename '-' here.  We want the event datafile
            # to be versioned.
            filecontents = open(filename, 'rb').read()
        return self.put(
                self.templates['event-detail-template'].format(graceid=graceid),
                files=[('eventFile', filename, filecontents)])

    def event(self, graceid):
        return self.get(
                self.templates['event-detail-template'].format(graceid=graceid))

    def events(self, query=None, orderby=None, count=None, columns=None):
        uri = self.links['events']
        qdict = {}
        if query:   qdict['query'] = query
        if count:   qdict['count'] = count
        if orderby: qdict['orderby'] = orderby
        if columns: qdict['columns'] = columns
        if qdict:
            uri += "?" + urllib.urlencode(qdict)
        while uri:
            response = self.get(uri).json()
            events = response.get('events',[])
            uri = response.get('links',{}).get('next')
            for event in events:
                 yield event

    def numEvents(self, query=None):
        uri = self.links['events']
        if query:
            uri += "?" + urllib.urlencode({'query': query})
        return self.get(uri).json()['numRows']

    def files(self, graceid, filename=None, raw=False):
        template = self.templates['files-template']
        uri = template.format(graceid=graceid, filename=filename or "")
        return self.get(uri)

    def writeFile(self, graceid, filename, filecontents=None):
        template = self.templates['files-template']
        uri = template.format(graceid=graceid, filename=os.path.basename(filename))
        if filecontents is None:
            if filename == '-':
                filename = 'stdin'
                filecontents = sys.stdin.read()
            else:
                filecontents = open(filename, "rb").read()
        elif isinstance(filecontents, file):
            # XXX Does not scale well.
            filecontents = filecontents.read()
        files = [('upload', os.path.basename(filename), filecontents)]
        return self.put(uri, files=files)

    def logs(self, graceid):
        template = self.templates['event-log-template']
        uri = template.format(graceid=graceid)
        return self.get(uri)

    def writeLog(self, graceid, message, filename=None, filecontents=None, 
            tagname=None, displayName=None):
        template = self.templates['event-log-template']
        uri = template.format(graceid=graceid)
        files = None
        if filename:
            if filecontents is None:
                if filename == '-':
                    filename = 'stdin'
                    filecontents = sys.stdin.read()
                else:
                    filecontents = open(filename, "rb").read()
            elif isinstance(filecontents, file):
                # XXX Does not scale well.
                filecontents = filecontents.read()
            files = [('upload', os.path.basename(filename), filecontents)]

        return self.post(uri, body={'message' : message, 'tagname': tagname, 
            'displayName': displayName}, files=files)

    def eels(self, graceid):
        """
        Return a list of EMBB log entries. 
        """
    
        template = self.templates['embb-event-log-template']
        uri = template.format(graceid=graceid)
        return self.get(uri)

    def writeEel(self, graceid, group, waveband, eel_status, 
            obs_status, **kwargs):
        """
        Write an EMBB event log entry. The following arguments are required:
            graceid, group, waveband, eel_status, obs_status.

        Additional keyword arguments may be passed in to be sent in the POST
        data. Only the following kwargs are recognized:
            ra
            dec
            raWidth
            decWidth
            gpstime
            duration
            comment
            extra_info_dict

        Any other kwargs will be ignored.
        """ 
        # validate facility, waveband, eel_status, and obs_status
        if not group in self.em_groups:
            raise ValueError("group must be one of %s" % self.em_groups)
        
        if not waveband in self.wavebands.keys():
            raise ValueError("waveband must be one of %s" % self.wavebands.keys())

        eel_status = self._getCode(eel_status, self.eel_statuses)
        if not eel_status:
            raise ValueError("EEL status must be one of %s" % self.eel_statuses.values())

        obs_status = self._getCode(obs_status, self.obs_statuses)
        if not obs_status:
            raise ValueError("Observation status must be one of %s" % self.obs_statuses.values())

        template = self.templates['embb-event-log-template']
        uri = template.format(graceid=graceid)

        body = {
            'group' : group,
            'waveband' : waveband,
            'eel_status' : eel_status,
            'obs_status' : obs_status,
        }
        body.update(**kwargs)
        return self.post(uri, body=body)

    def labels(self, graceid, label=""):
        template = self.templates['event-label-template']
        uri = template.format(graceid=graceid, label=label)
        return self.get(uri)

    def writeLabel(self, graceid, label):
        template = self.templates['event-label-template']
        uri = template.format(graceid=graceid, label=label)
        return self.put(uri)

    def removeLabel(self, graceid, label):
        template = self.templates['event-label-template']
        uri = template.format(graceid=graceid, label=label)
        return self.delete(uri)

    def tags(self, graceid, n):
        template = self.templates['taglist-template']
        uri = template.format(graceid=graceid, n=n)
        return self.get(uri)

    def createTag(self, graceid, n, tagname, displayName=None):
        template = self.templates['tag-template']
        uri = template.format(graceid=graceid, n=n, tagname=tagname)
        return self.put(uri, body={'displayName': displayName})

    def deleteTag(self, graceid, n, tagname):
        template = self.templates['tag-template']
        uri = template.format(graceid=graceid, n=n, tagname=tagname)
        return self.delete(uri)

    def ping(self):
        return self.get(self.links['self'])

#-----------------------------------------------------------------
# TBD
# Media Types

# Root

# Collection

# Event

# Log

# File

# Label

#-----------------------------------------------------------------
# HTTP upload encoding
# Taken from http://code.activestate.com/recipes/146306/

def encode_multipart_formdata(fields, files):
    """
    fields is a sequence of (name, value) elements for regular form fields.
    files is a sequence of (name, filename, value) elements for data to be uploaded as files
    Return (content_type, body) ready for httplib.HTTP instance
    """
    BOUNDARY = '----------ThIs_Is_tHe_bouNdaRY_$'
    CRLF = '\r\n'
    L = []
    for (key, value) in fields:
        if value is None: continue
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"' % key)
        L.append('')
        # str(value) in case it is unicode
        L.append(str(value))
    for (key, filename, value) in files:
        if value is None: continue
        L.append('--' + BOUNDARY)
        # str(filename) in case it is unicode
        L.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (key, str(filename)))
        L.append('Content-Type: %s' % get_content_type(filename))
        L.append('')
        L.append(value)
    L.append('--' + BOUNDARY + '--')
    L.append('')
    body = CRLF.join(L)
    content_type = 'multipart/form-data; boundary=%s' % BOUNDARY
    return content_type, body

def get_content_type(filename):
    return mimetypes.guess_type(filename)[0] or 'application/octet-stream'

#-----------------------------------------------------------------
# X509 Credentials

# XXX No longer checking whether this is a pre-RFC proxy.  Is this okay?
def findUserCredentials(warnOnOldProxy=1):

    proxyFile = os.environ.get('X509_USER_PROXY')
    certFile = os.environ.get('X509_USER_CERT')
    keyFile = os.environ.get('X509_USER_KEY')

    if certFile and keyFile:
        return certFile, keyFile

    if proxyFile:
        return proxyFile, proxyFile

    # Try default proxy
    proxyFile = os.path.join('/tmp', "x509up_u%d" % os.getuid())
    if os.path.exists(proxyFile):
        return proxyFile, proxyFile

    # Try default cert/key
    homeDir = os.environ.get('HOME')
    certFile = os.path.join(homeDir, '.globus', 'usercert.pem')
    keyFile = os.path.join(homeDir, '.globus', 'userkey.pem')

    if os.path.exists(certFile) and os.path.exists(keyFile):
        return certFile, keyFile


