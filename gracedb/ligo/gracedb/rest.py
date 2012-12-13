
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


import httplib
import mimetypes
import urllib
import os
import sys
import json

DEFAULT_SERVICE_URL = "https://gracedb.ligo.org/api"

#-----------------------------------------------------------------
# Exception(s)

class HTTPError(Exception):
    def __init__(self, status, reason, message):
        self.status = status
        self.reason = reason
        self.message = message
        Exception.__init__(self, (status, reason))

#-----------------------------------------------------------------
# Generic GSI REST

class GsiRest(object):
    """docstring for GracedbRest"""
    def __init__(self, cred=None):
        if not cred:
            cred = findUserCredentials()
        if isinstance(cred, (list, tuple)):
            self.cert, self.key = cred
        elif cred:
            self.cert = self.key = cred
        self.connection = None

    def getConnection(self, url):
        proto, rest = urllib.splittype(url)
        host, path = urllib.splithost(rest)
        host, port = urllib.splitport(host)
        return httplib.HTTPSConnection(host, port,
                key_file=self.key, cert_file=self.cert)

    def request(self, method, url, body=None, headers=None, priming_url=None):
        conn = self.getConnection(url)
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
                    headers['content-type'] = "application/x-www-form-urlencoded"
                body = urllib.urlencode(body)
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
    def __init__(self, service_url=DEFAULT_SERVICE_URL, *args, **kwargs):
        GsiRest.__init__(self, *args, **kwargs)

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
    def analysis_types(self):
        return self.service_info.get('analysis-types')

    @property
    def groups(self):
        return self.service_info.get('groups')

    def _analysisTypeCode(self, analysis_type):
        """Check if analysis_type is valid.
           Return coded version of type if it is"""
        # types is dict of { code : descriptive_name }
        types = self.analysis_types
        if analysis_type in types.keys():
            # Already coded
            return analysis_type
        if not analysis_type in types.values():
            # Not valid
            return None
        return [code
                for code, name in types.items()
                if name == analysis_type][0]

    def createEvent(self, group, analysis_type, filename, filecontents=None):
        errors = []
        analysis_type_code = self._analysisTypeCode(analysis_type)
        if group not in self.groups:
            errors += ["bad group"]
        if not analysis_type_code:
            errors += ["bad type"]
        if errors:
            # XXX Terrible error messages / weak exception type
            raise Exception(str(errors))
        if filecontents is None:
            filecontents = open(filename, 'rb').read() # XXX or 'rb' ?
        fields = [
                  ('group', group),
                  ('type', analysis_type_code),
                 ]
        files = [('eventFile', filename, filecontents)]
        # Python httplib bug?  unicode link
        uri = str(self.links['events'])
        return self.post(uri, fields, files=files)

    def replaceEvent(self, graceid, filename, filecontents=None):
        if filecontents is None:
            filecontents = open(filename, 'rb').read()
        return self.put(
                self.templates['event-detail-template'].format(graceid=graceid),
                files=[('eventFile', filename, filecontents)])

    def event(self, graceid):
        return self.get(
                self.templates['event-detail-template'].format(graceid=graceid))

    def events(self, query=None, orderby=None, count=None):
        uri = self.links['events']
        qdict = {}
        if query:   qdict['query'] = query
        if count:   qdict['count'] = count
        if orderby: qdict['orderby'] = orderby
        if qdict:
            uri += "?" + urllib.urlencode(qdict)
        while uri:
            response = self.get(uri).json()
            events = response['events']
            uri = response['links'].get('next')
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

    def writeLog(self, graceid, message):
        template = self.templates['event-log-template']
        uri = template.format(graceid=graceid)
        return self.post(uri, body={'message' : message})

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
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"' % key)
        L.append('')
        L.append(value)
    for (key, filename, value) in files:
        L.append('--' + BOUNDARY)
        L.append('Content-Disposition: form-data; name="%s"; filename="%s"' % (key, filename))
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


