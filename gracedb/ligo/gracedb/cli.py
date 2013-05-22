#!/usr/bin/env python

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

import mimetypes, urllib
import os, sys, shutil
import json
from rest import GraceDb

DEFAULT_SERVICE_URL = "https://gracedb.ligo.org/gracedb/cli"

GIT_TAG = 'gracedb-1.11-1'

#-----------------------------------------------------------------
# Util routines

def error(*message):
    message = "".join(message) + "\n"
    sys.stderr.write("ERROR: " + message)

def warning(*message):
    message = "".join(message) + "\n"
    sys.stderr.write("WARNING: " + message)

def output(*message):
    message = "".join(message) + "\n"
    sys.stdout.write(message)

#-----------------------------------------------------------------
# HTTP upload encoding
# Taken from http://code.activestate.com/recipes/146306/

typeCodeMap = {
        "LowMass" : "LM",
        "HighMass" : "HM",
        "GRB"      : "GRB",
        "Ringdown" : "RD",
        "Omega"    : "OM",
        "Q"        : "Q",
        "X"        : "X",
        "CWB"      : "CWB",
        "MBTAOnline": "MBTA",
        "Injection": "HWINJ",
}
validTypes = typeCodeMap.keys()

def post_multipart(h, selector, fields, files):
    """
    Post fields and files to an http host as multipart/form-data.
    fields is a sequence of (name, value) elements for regular form fields.
    files is a sequence of (name, filename, value) elements for data to be uploaded as files
    Return the server's response page.
    """
    content_type, body = encode_multipart_formdata(fields, files)
    #h = httplib.HTTP(host)
    h.putrequest('POST', selector)
    h.putheader('content-type', content_type)
    h.putheader('content-length', str(len(body)))
    h.endheaders()
    h.send(body)
    errcode, errmsg, headers = h.getreply()
    return h.file.read()

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
        # str(value) in case it is unicode
        L.append(str(value))
    for (key, filename, value) in files:
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
# Web Service Client

# XXX Replace with something that looks like this.  
# the methods from GraceDb as necessary.
class Client(GraceDb):
    def __init__(self,
            url=DEFAULT_SERVICE_URL, 
            proxy_host=None, proxy_port=3128,
            credentials=None,
            *args, **kwargs):
        if (url[-1] != '/'):
            url += '/'
        super(Client, self).__init__(self, url, proxy_host, proxy_port, 
                                        credentials, *args, **kwargs)

# XXX break here, branson
    def _connect(self):
        self._conn = self.connector()

    def _send(self, method, httpmethod="POST", **kw):
        try:
            kw['cli'] = 'true'
            kw['cli_version'] = "1"
            kw = urllib.urlencode(kw)
            headers = {'connection' : 'keep-alive'}
            url = "%s/%s" % (self.url, method)
            self._connect()
            self._conn.request(httpmethod, url, kw, headers)
            response = self._conn.getresponse()
            rv = response.read()
        except Exception, e:
            return { 'error':  "client send exception: " + str(e) }

        # XXX ridiculous hack to deal with downloaded ligolw search results.
        if response.getheader('content-type') == 'application/xml':
            return rv

        try:
            # XXX Bad!  Should be JSON conversion
            # XXX Well, I could make it so.  Might as well, right?  But I want
            # the responses to look like they did before.
            return eval(rv)
        except Exception, e:
            if "Authorization Required" in rv:
                return { 'error': 'Credentials not accepted' }
            return {'error': "while parsing:%s\nclient send exception: %s" % (rv, str(e))}

    def _is_absolute_url(self, url):
        _, junk = urllib.splittype(url)
        host, _ = urllib.splithost(junk)
        return host is not None

    def rest(self, resource, method="GET", **kw):
        headers = kw.pop('headers', {})
        headers['connection'] = 'keep-alive'
#       If absolute URL, use, otherwise pre-append the REST service URL.
        if self._is_absolute_url(resource):
            url = resource
        else:
            url = "%s%s" % (self.rest_url, resource)
        kw = urllib.urlencode(kw)
        self._connect()
        self._conn.request(method, url, kw, headers)
        return self._conn.getresponse()

    def _upload(self, method, fields, files,
            alert=False, http_method="POST", rawresponse=False):
        # Do an tiny GET request to get SSL primed.
        # Large initial SSL/POST requests will choke the server.
        #r = self._send('ping', 'GET', ack='priming')
        r = self._send('ping', ack='priming')
        if 'error' in r and r['error']:
            return r
        try:
            fields += [('cli', 'true')]
            fields += [('cli_version', "1")]
            fields += [('alert', str(alert))]
            content_type, body = encode_multipart_formdata(fields, files)
            headers = {
                'content-type': content_type,
                'content-length': str(len(body)),
                'connection': 'keep-alive',
            }
            if self._is_absolute_url(method):
                url = method
            else:
                url = "%s/%s" % (self.url, method)
            self._conn.request(http_method, url, body, headers)
            response = self._conn.getresponse()
            if rawresponse:
                return response
            rv = response.read()
        except Exception, e:
            return { "error" : "client upload exception: " + str(e) }
        # XXX should be JSON conversion
        try:
            return eval(rv)
        except Exception, e:
            return {'error': "while parsing:%s\nclient upload exception: %s" % (rv, str(e))}

    def ping(self, msg=""):
        return self._send('ping', ack=msg)

    def search(self, query, columns=None, ligolw=False):
        terms = { "query" : query }
        if columns:
            terms['columns'] = columns
        if ligolw:
            terms['ligolw'] = 1
        return self._send('search', **terms)

    def log(self, graceid, message, alert=False):
        return self._send('log', graceid=graceid, message=message, alert=alert)

    def label(self, graceid, label, alert=False):
        return self._send('label', graceid=graceid, label=label, alert=alert)

    def create(self, group, analysis_type, filename, filecontents=None):
        if analysis_type in typeCodeMap:
            analysis_type = typeCodeMap[analysis_type]
        if filecontents is None:
            if filename == '-':
                filename = 'initial.data'
                filecontents = sys.stdin.read()
            else:
                filecontents = open(filename, 'r').read()
        fields = [
                  ('group', group),
                  ('type', analysis_type),
                 ]
        files = [('eventFile', filename, filecontents)]
        return self._upload('create', fields, files)

    def replace(self, graceid, filename, filecontents=None):
        from ligo.gracedb.rest import GraceDb
        url = self.rest_url
        if url[-1] != '/':
            url += '/'
        server = GraceDb(url)
        response = server.replaceEvent(graceid, filename, filecontents)
        if response.status == 202:  # Accepted
            return "%s updated" % graceid
        else:
            return {'error': "Bad response %s / %s\n%s" %
                    (response.status, response.reason, response.read())}

    def upload(self, graceid, filename, filecontents=None, comment="", alert=False):
        if filecontents is None:
            if filename == '-':
                filename = 'stdin'
                filecontents = sys.stdin.read()
            else:
                filecontents = open(filename, 'r').read()
        fields = [
            ('graceid', graceid),
            ('comment', comment),
            ('filename', filename)
        ]
        files = [ ('upload', filename, filecontents) ]
        return self._upload('upload', fields, files, alert)

    def listfiles(self, graceid):
        response = self.rest('/event/%s/files/' % graceid)
        return response

    def download(self, graceid, filename, destfile):
        # Check that we *could* write the file before we
        # go to the trouble of getting it.  Also, try not
        # to open a file until we know we have data.
        if not isinstance(destfile, file) and destfile != "-":
            if not os.access(os.path.dirname(os.path.abspath(destfile)), os.W_OK):
                raise IOError("%s: Permission denied" % destfile)
        response = self.rest('/event/%s/files/%s' % (graceid, filename))
        if response.status == 200:
            if not isinstance(destfile, file):
                if destfile == '-':
                    destfile = sys.stdout
                else:
                    destfile = open(destfile, "w")
            shutil.copyfileobj(response, destfile)
            return 0
        else:
            return "Error. (%d) %s" % (response.status, response.reason)

#-----------------------------------------------------------------
# Main 


def main():
    usage ="""%%prog [options] GROUP TYPE EVENTFILE
   where GROUP is one of %(groups)s
         TYPE is one of %(types)s
         EVENTFILE is file containing event data. '-' indicates stdin.

%%prog [options] replace GRACEID EVENTFILE
   where GROUP is one of %(groups)s
         TYPE is one of %(types)s
         EVENTFILE is file containing event data. '-' indicates stdin.

%%prog [options] ping
   Test server connection

%%prog [options] upload GRACEID FILE [COMMENT] [TAG_NAME] [DISP_NAME]
   where GRACEID is the id of an existing candidate event in GraCEDb
         FILE      is the name of the file to upload. '-' indicates stdin.
         COMMENT   is an optional annotation to enter into the log
         TAG_NAME  is the name of a tag to add to the log message
         DISP_NAME is the display name of the tag (use for new tags only)
   Upload FILE to the private data area for a candidate event

%%prog [options] download GRACEID FILE [DESTINATION]
   where GRACEID      is the id of an existing candidate event in GraCEDb
         FILE         is the name of the file previosuly uploaded.
         DESTINATION  is the download destination.  '-' indicates stdout.
                      default is same file name as FILE
    Download FILE from private data area of a candidate event

%%prog [options] log GRACEID COMMENT
   where GRACEID  is the id of an existing candidate event in GraCEDb
         COMMENT  is text that will be entered into the event's log
   Enter a comment into the log for a candidate event

%%prog [options] label GRACEID LABEL
    Label event with GRACEDID with LABEL.  LABEL must already exist.

%%prog [options] tag GRACEID LOG_N TAG_NAME [DISP_NAME]
    Tag an existing log message.
    LOG_N is the number of the log message.

%%prog [options] delete_tag GRACEID LOG_N TAG_NAME
    Remove a tag from a log message.

%%prog [options] search SEARCH PARAMS
    Search paramaters are a list of requirements to be satisfied.  They
    may be GPS times, GPS time ranges, graceids and ranges, group(s),
    analysis type(s), labels, etc.  Note that text is case insensitive
    Example: %%prog search G0100..G0200 mbta LUMIN_GO

Environment Variables:
    GRACEDB_SERVICE_URL   (can be overridden by --service-url)
    HTTP_PROXY            (can be overridden by --proxy)
    X509_USER_PROXY
    X509_USER_CERT
    X509_USER_KEY

Credentials are looked for in this order:
    (1) $(X509_USER_CERT) / $(X509_USER_KEY)
    (2) $(X509_USER_PROXY)
    (3) Default location of grid proxy ( /tmp/x509up_u$(UID) )
    (4) $(HOME)/.globus/usercert.pem / $(HOME)/.globus/userkey.pem

Note that comments can only be 200 characters long.
Longer strings will be truncated.""" % {
        'groups' : 'CBC, Burst, Stochastic, CW',
        'types'  : ", ".join(validTypes),
    }

    from optparse import OptionParser
    op = OptionParser(usage=usage)
    op.add_option("-p", "--proxy", dest="proxy",
                  help="HTTP Proxy", metavar="PROXY[:PORT]")
    op.add_option("-s", "--service-url", dest="service",
                  help="GraCEDb Service URL", metavar="URL")
    op.add_option("-f", "--filename", dest="filename",
                  help="If data is read from stdin, use this as the filename.", metavar="NAME")

    op.add_option("-a", "--alert", dest="alert",
                  help="Send an LV alert (not meaningful for search, implied for create)",
                  action="store_true", default=False
                 )

    op.add_option("-c", "--columns", dest="columns",
                  help="Comma separated list of event attributes to include in results (only meaningful in search)",
                  default=None
                 )

    op.add_option("-l", "--ligolw", dest="ligolw",
                  help="Download ligolw file of combined search results (not meaningful outside of search)",
                  action="store_true", default=False
                 )

    options, args = op.parse_args()

    proxy = options.proxy or os.environ.get('HTTP_PROXY', None)
    service = options.service or \
              os.environ.get('GRACEDB_SERVICE_URL', None) or \
              DEFAULT_SERVICE_URL

    proxyport = None
    if proxy and proxy.find(':') > 0:
        try:
            proxy, proxyport = proxy.split(':')
            proxyport = int(proxyport)
        except:
            op.error("Malformed proxy: '%s'" % proxy)
    if proxyport:
        client = Client(service,
                        proxy_host=proxy,
                        proxy_port=proxyport)
    else:
        client = Client(service, proxy_host=proxy)

    if len(args) < 1:
        op.error("not enough arguments")
    elif args[0] == 'ping':
        # XXX Branson: need a ping method in the REST API.
        # Also need to add alert options through the API?
        msg = " ".join(args[1:]) or "PING"
        response = client.ping(msg)
    elif args[0] == 'upload':
        if len(args) < 3:
            op.error("not enough arguments for upload")
        graceid = args[1]
        filename = args[2]
        # In the previous version, you passed in comment and alert options.
        # XXX What do we need to do now?
        # comment = " ".join(args[3:])
        # XXX How does the response differ in these two cases?
        # response = client.upload(graceid, filename, comment=comment, alert=options.alert)
        # The old one was whatever you get from 
        # httplib.HTTPConnection().getresponse()
        # i.e., an httplib.HTTPResponse instance.
        # You then return whatever you would've gotten from response.read(), which is 
        # just the response body. Since we're posting to /cli/upload, we look to see
        # what the function gracedb.view.upload returns.
        # The body is a plaintext message.
        # ERROR: missing arg(s)
        # ERROR: Event does not exist.
        # OK
        # ERROR: problem creating log engry
        # ERROR: could not save file.
        # Can you translate the API output to look like this.  Check status code and
        # do something with it.
        
        # Actually you should make this a log message.
        response = client.writeFile(graceid, filename)
    elif args[0] == 'download':
        if len(args) not in [2,3,4]:
            op.error("not enough arguments for download")
        graceid = args[1]
        if len(args) == 2:
            # get/print listing.
            response = client.listfiles(graceid)
            if response and response.status == 200:
                for fname in json.loads(response.read()):
                    print(fname)
                exit(0)
            print(response.reason)
            exit(1)
        filename = args[2]
        if len(args) == 4:
            outfile = args[3]
        else:
            outfile = os.path.basename(filename)
        response = client.download(graceid, filename, outfile)
        if response:
            # no response means file saved.  any other response is an error message.
            print response
            exit(1)
        exit(0)
    elif args[0] == 'log':
        if len(args) < 3:
            op.error("not enough arguments for log")
        graceid = args[1]
        message = " ".join(args[2:])
        response = client.log(graceid, message, alert=options.alert)
    elif args[0] == 'tag':
        if len(args) not in [3,4]:
            op.error("wrong number of arguments for tag")
        graceid = args[1]
        logN = args[2]
        tagName = args[3]
        dispName = args[4]
        response = client.createTag(graceid, logN, tagName, dispName)
        # XXX Umm... how should I communicate errors.  Print?
        # Were they coming from the cli views before?
        if response.status==404:
            print "ERROR: %s" % response.read()
        elif response.status==409:
            print "FAILED: %s" % response.read()
        elif response.status==201:
            print "Tag operation successful."
    elif args[0] == 'delete_tag':
        if len(args) != 3:
            op.error("wrong number of arguments for delete_tag")
        graceid = args[1]
        logN = args[2]
        tagName = args[3]
        response = client.delete_tag(graceid, logN, tagName, dispName)
# XXX Branson: so far worked on the label action.
    elif args[0] == 'label':
        if len(args) != 3:
            op.error("wrong number of arguments for label")
        graceid = args[1]
        label = args[2]
        response = client.writeLabel(graceid, label, alert=options.alert)
        rdict = response.json()
        if response.status == 201:
            for key in rdict:
                op.error("%s" % rdict[key])
        elif response.status == 404:
            op.error("ERROR: Event matching query does not exist.")
        elif response.status == 400:
            op.error("ERROR: No such label '%s'" % label)
# XXX Branson crap.  How are you going to make stuff into ligolw?
# Look at views.cli_search.  Still how is the best way to handle this?
# Create xml file on server side?
# Another option would be to overwrite the client search function 
# to be the same as in the old client. Then forget about it.
    elif args[0] == 'search':
        response = client.search(" ".join(args[1:]), options.columns, options.ligolw)
        #response = client.events(" ".join(args[1:]))
    elif args[0] == 'replace':
        if len(args) != 3:
            op.error("wrong number of args for replace")
        graceid = args[1]
        filename = args[2]
        response = client.replace(graceid, filename)
    elif args[0] == 'slot':
        if len(args) in [3,4]:
            url = '/events/%s/slot/%s' % (args[1], args[2])
            if len(args) == 3:
                response = client.rest(url, method='GET')
            else:
                headers = {'content-type' : "application/x-www-form-urlencoded" }
                response = client.rest(url,
                        method='PUT',
                        headers=headers,
                        filename=args[3])
            if response and response.status in [200, 201]:
                print( json.loads(response.read()) ) 
                exit(0)
            else:
                print(response.reason)
                exit(1)
        else:
            op.error("wrong number of args for slot")
    elif len(args) == 3:
        group = args[0]
        type = args[1]
        filename = args[2]

        if type in typeCodeMap:
            type = typeCodeMap[type]
        else:
            error("Type must be one of: ", ", ".join(typeCodeMap.keys()))
            sys.exit(1)
        response = client.create(group, type, filename)
        if not response:
            error("There was a problem.  Did you do grid-proxy-init -rfc?")
            sys.exit(1)
    else:
        op.error("")
        sys.exit(1)

    # Output the response.

    exitCode = 0
    # XXX oddball exception for ligolw query responses.
    if isinstance(response, str):
        print(response)
    else:
        if ('error' in response) and response['error']:
            error(response['error'])
            exitCode = 1
        if ('warning' in response) and response['warning']:
            warning(response['warning'])
        if ('output' in response) and response['output']:
            output(response['output'])

    return exitCode

if __name__ == "__main__":
    try:
        code = main()
    except Exception, e:
        error(str(e))
        sys.exit(1)
    sys.exit(code)
