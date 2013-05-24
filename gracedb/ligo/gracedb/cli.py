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

import os, sys, shutil
import json
from rest import GraceDb

DEFAULT_SERVICE_URL = "https://gracedb.ligo.org/gracedb/api"

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

# NB:  We are not checking against this list anymore.  Instead, we will
# get the list of groups and types from the API Root.  However, this is
# left in so that the docstring will still work even without a connection
# to the API.
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


#-----------------------------------------------------------------
# Web Service Client

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

    def search(self, query, columns=None, ligolw=False):
        terms = { "query" : query }
        if columns:
            terms['columns'] = columns
        if ligolw:
            terms['ligolw'] = 1
        headers = {'connection' : 'keep-alive'}

        # construct url via ugly hack.
        url = self.url[:self.url.rindex('api')] + "cli/search"

        try:
            response = self.post(url,body=terms,headers=headers)
            rv = response.read()
            if "Authorization Required" in rv:
                return { 'error': 'Credentials not accepted' }
        except Exception, e:
            return { 'error':  "client send exception: " + str(e) }

        # XXX ridiculous hack to deal with downloaded ligolw search results.
        if response.getheader('content-type') == 'application/xml':
            return rv
        else:
            try:
                rv = json.loads(rv)  
                return rv
            except Exception, e:
                return {'error': "while parsing:%s\nclient send exception: %s" % (rv, str(e))}


    def download(self, graceid, filename, destfile):
        # Check that we *could* write the file before we
        # go to the trouble of getting it.  Also, try not
        # to open a file until we know we have data.
        if not isinstance(destfile, file) and destfile != "-":
            if not os.access(os.path.dirname(os.path.abspath(destfile)), os.W_OK):
                raise IOError("%s: Permission denied" % destfile)
        response = self.files(graceid, filename)
        if response.status == 200:
            if not isinstance(destfile, file):
                if destfile == '-':
                    destfile = sys.stdout
                else:
                    destfile = open(destfile, "w")
            # XXX Check.  This is a django rest framework response object.
            # We want to get the body. Hence .read() 
            shutil.copyfileobj(response.read(), destfile)
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

%%prog [options] upload GRACEID FILE [COMMENT] 
   where GRACEID is the id of an existing candidate event in GraCEDb
         FILE      is the name of the file to upload. '-' indicates stdin.
         COMMENT   is an optional annotation to enter into the log
   Upload FILE to the private data area for a candidate event. To apply 
   a tag, use the --tag-name option (and --tag-display-name if desired.)

%%prog [options] download GRACEID FILE [DESTINATION]
   where GRACEID      is the id of an existing candidate event in GraCEDb
         FILE         is the name of the file previosuly uploaded.
         DESTINATION  is the download destination.  '-' indicates stdout.
                      default is same file name as FILE
    Download FILE from private data area of a candidate event

%%prog [options] log GRACEID COMMENT
   where GRACEID  is the id of an existing candidate event in GraCEDb
         COMMENT  is text that will be entered into the event's log
   Enter a comment into the log for a candidate event.  To apply a tag,
   use the --tag-name option (and --tag-display-name if desired).

%%prog [options] label GRACEID LABEL
    Label event with GRACEDID with LABEL.  LABEL must already exist.

%%prog [options] tag GRACEID LOG_N TAG_NAME [DISP_NAME]
   where GRACEID   is the id of an existing candidate event in GraCEDb
         LOG_N     is the number of the log message.
         TAG_NAME  is the name of the tag
         DISP_NAME is the tag display name (ignored for existing tags)
    Tag an existing log message. Alternatively, the tag name and 
    display name can be passed in with the --tag-name and 
    --tag-display-name options.

%%prog [options] delete_tag GRACEID LOG_N TAG_NAME
    Remove a tag from a log message. Alternatively, the tag name 
    can be passed in with the --tag-name option.

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
        'groups' : 'CBC, Burst, Stochastic, Coherent, Test, External',
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
    op.add_option("-t", "--tag-name", dest="tagName",
                  help="tag name in database (only used for log, upload, tag, and delete_tag)",
                  action="store_true", default=False
                 )
    op.add_option("-t", "--tag-display-name", dest="tagDispName",
                  help="tag display name (ignored for existing tags)",
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
        # what about the response to this?
        response = client.ping()
    elif args[0] == 'upload':
        if len(args) < 3:
            op.error("not enough arguments for upload")
        graceid = args[1]
        filename = args[2]
        comment = " ".join(args[3:])
        tagName = options.tagName
        tagDispName = options.tagDispName
        # XXX How does the old response look?
        # ERROR: missing arg(s)
        # ERROR: Event does not exist.
        # OK
        # ERROR: problem creating log engry
        # ERROR: could not save file.
        
        # XXX Fix this.
        op.warning("comment ignored: %s" % comment)
        response = client.writeFile(graceid, filename)
        # response = client.writeLog(graceid, filename, comment, 
        #    tagName, tagDispName,  options.alert)
    elif args[0] == 'download':
        if len(args) not in [2,3,4]:
            op.error("not enough arguments for download")
        graceid = args[1]
        if len(args) == 2:
            # get/print listing.
            response = client.files(graceid)
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
        # XXX need alert option in here.
#        response = client.writeLog(graceid, message, options.tagName, options.tagDispName, 
#            alert=options.alert)
        response = client.writeLog(graceid, message, options.tagName, options.tagDispName)
    elif args[0] == 'tag':
        if options.tagName:
            if len(args) != 2:
                op.error("wrong number of arguments for tag")
            tagName = options.tagName
            tagDispName = options.tagDispName
        else:
            if len(args) not in [3,4]:
                op.error("wrong number of arguments for tag")
            tagName = args[3]
            tagDispName = args[3]
        graceid = args[1]
        logN = args[2]
        response = client.createTag(graceid, logN, tagName, tagDispName)
    elif args[0] == 'delete_tag':
        if options.tagName:
            if len(args) != 2:
                op.error("wrong number of arguments for delete_tag")
            tagName = options.tagName
        else:
            if len(args) != 3:
                op.error("wrong number of arguments for delete_tag")
            tagName = args[3]
        graceid = args[1]
        logN = args[2]
        response = client.deleteTag(graceid, logN, tagName)
    elif args[0] == 'label':
        if len(args) != 3:
            op.error("wrong number of arguments for label")
        graceid = args[1]
        label = args[2]
        response = client.writeLabel(graceid, label, alert=options.alert)
    elif args[0] == 'search':
        response = client.search(" ".join(args[1:]), options.columns, options.ligolw)
    elif args[0] == 'replace':
        if len(args) != 3:
            op.error("wrong number of args for replace")
        graceid = args[1]
        filename = args[2]
        response = client.replaceEvent(graceid, filename)
    elif len(args) == 3:
        # Create a new event.
        group = args[0]
        type = args[1]
        filename = args[2]

        # Check that the group and type are known to the API.
        # NB: the dictionary returned by the API has keys and values
        # reversed w.r.t. the typeCodeMap above.
        foundType = False
        for key, value in client.analysis_types:
            if type==str(value):
                type = key
                foundType = True
        if not foundType:
            error("Type must be one of: ", ", ".join(client.analysis_types.values()))
            sys.exit(1)

        foundGroup = True if (unicode(group) in client.groups) else False
        if not foundGroup:
            error("Group must be one of: ", ", ".join(client.groups))
            sys.exit(1)

        response = client.createEvent(group, type, filename)
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
        rv = response.read()
        # XXX If you got this far, you should be JSON.
        responseBody = json.loads(rv)
        status = response.status
        if status >= 400:
            exitCode=1
        if isinstance(responseBody, str):
            output("%d: %s" % (status, responseBody))
        else:
            # XXX Think about what will happen here with ping. Right now just 200.
            output("Server returned %d" % status)
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
