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

import os, sys, shutil, urllib
import json
from ligo.gracedb.rest import GraceDb


DEFAULT_SERVICE_URL = "https://gracedb.ligo.org/gracedb/api"

GIT_TAG = 'gracedb-1.15-1'

DEFAULT_COLUMNS = "graceid,labels,group,pipeline,search,far,gpstime,created,dataurl"
 
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

def defaultAccess(e,a):
    if a.find('.') < 0:
        return str(e.get(a,""))
    rv = e
    attrs = a.split('.')
    while attrs and rv:
        rv = rv.get(attrs[0],"")
        attrs = attrs[1:]
    return str(rv)

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
        self.url = url
        super(Client, self).__init__(url, proxy_host, proxy_port, 
                                        credentials, *args, **kwargs)

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
            shutil.copyfileobj(response, destfile)
            return 0
        else:
            return "Error. (%d) %s" % (response.status, response.reason)

    # Hamstring 'adjustResponse' from the example REST client.
    # We don't want it messing with the response from the server.
    def adjustResponse(self, response):
        response.json = lambda: json.loads(response.read())
        return response

#-----------------------------------------------------------------
# Main 


def main():
    usage ="""%%prog [options] GROUP PIPELINE SEARCH EVENTFILE
   where GROUP is one of %(groups)s
         PIPELINE is one of %(pipelines)s
         SEARCH (optional) is one of %(searches)s
         EVENTFILE is file containing event data. '-' indicates stdin.
    NOTE: the groups, pipelines and searches in this docstring may not
    be up top date. To see an accurate list, do:

    %%prog list groups
    %%prog list pipelines
    %%prog list searches

%%prog [options] replace GRACEID EVENTFILE
   where GROUP is one of %(groups)s
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
        'groups'     : 'CBC, Burst, Stochastic, Coherent, Test, External',
        'pipelines'  : 'MBTAOnline, gstlal, gstlal-spiir, HardwareInjection, Fermi, Swift, CWB, CWB2G',
        'searches'   : 'AllSky, LowMass, HighMass, GRB, Test',
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
                  help="Send an LV alert (deprecated; alerts sent by default)",
                  action="store_true", default=None
                 )

    op.add_option("-c", "--columns", dest="columns",
                  help="Comma separated list of event attributes to include in results (only meaningful in search)",
                  default=DEFAULT_COLUMNS
                 )

    op.add_option("-l", "--ligolw", dest="ligolw",
                  help="Download ligolw file of combined search results (not meaningful outside of search). NOTE: Produces an ERROR if any of the events returned by the search do not have coinc.xml files.",
                  action="store_true", default=False
                 )
    op.add_option("-t", "--tag-name", dest="tagName",
                  help="tag name in database (only used for log, upload, tag, and delete_tag)",
                  default=None
                 )
    op.add_option("-d", "--tag-display-name", dest="tagDispName",
                  help="tag display name (ignored for existing tags)",
                  default=None
                 )

    options, args = op.parse_args()

    try:
        from glue.ligolw import ligolw
        from glue.ligolw import lsctables
        from glue.ligolw import utils
        from glue.ligolw.utils import ligolw_add
    except:
        if options.ligolw:
            error("ligolw modules not found")
            exit(1)
        else:
            pass

    proxy = options.proxy or os.environ.get('HTTP_PROXY', None)
    service = options.service or \
              os.environ.get('GRACEDB_SERVICE_URL', None) or \
              DEFAULT_SERVICE_URL

    if options.alert is not None:
        warning("alert option is deprecated.  Alerts are now sent by default.")

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
        response = client.ping()
        output("Client groups: %s" % client.groups)
        output("Client pipelines: %s" % client.pipelines)
        output("Client searches: %s" % client.searches)
        if response.status==200:
            output("%s: 200 OK" % service)
            exit(0)
    elif args[0] == 'list':
        if args[1] == 'groups':
            output(client.groups)
            exit(0)
        elif args[1] == 'pipelines':
            output(client.pipelines)
            exit(0)
        elif args[1] == 'searches':
            output(client.searches)
            exit(0)
        else:
            output("Unknown list object. Please use 'groups', 'pipelines', or 'searches.'")
            exit(1)
    elif args[0] == 'upload':
        if len(args) < 3:
            op.error("not enough arguments for upload")
        graceid = args[1]
        filename = args[2]
        comment = " ".join(args[3:])
        tagName = options.tagName
        tagDispName = options.tagDispName
        response = client.writeLog(graceid, comment, filename, None,
            tagName, tagDispName)
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
        response = client.writeLog(graceid, message, options.tagName, options.tagDispName)
    elif args[0] == 'tag':
        if options.tagName:
            if len(args) != 2:
                op.error("wrong number of arguments for tag")
            tagName = options.tagName
            tagDispName = options.tagDispName
        else:
            if len(args) not in [4,5]:
                op.error("wrong number of arguments for tag")
            tagName = args[3]
            tagDispName = None
            if len(args)==5:
                tagDispName = args[4]
        graceid = args[1]
        logN = args[2]
        response = client.createTag(graceid, logN, tagName, tagDispName)
    elif args[0] == 'delete_tag':
        error("len args = %s" % len(args))
        error("args = %s" % args)
        if options.tagName:
            if len(args) != 2:
                op.error("wrong number of arguments for delete_tag")
            tagName = options.tagName
        else:
            if len(args) != 4:
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
        response = client.writeLabel(graceid, label)
    elif args[0] == 'search':
        query = " ".join(args[1:])
        
        columns = options.columns
        columns = columns.replace('DEFAULTS',DEFAULT_COLUMNS)
        columns = columns.split(',')

        count = None # XXX Let's just get rid of this?
        orderby = None # XXX Should we implement this?

        events = client.events(query, orderby, count, columns)

        if options.ligolw:
            xmldoc = ligolw.Document()
            for e in events:
                graceid = e['graceid']
                try:
                    r = client.files(graceid, "coinc.xml")
                    utils.load_fileobj(r, xmldoc = xmldoc)
                except:
                    error("Missing coinc.xml for %s. Cannot build ligolw output." % graceid)
                    exit(1)
            ligolw_add.reassign_ids(xmldoc)
            ligolw_add.merge_ligolws(xmldoc)
            ligolw_add.merge_compatible_tables(xmldoc)
            xmldoc.write()
        else:
            accessFun = {
                "labels" : lambda e: \
                    ",".join(e['labels'].keys()),
                "dataurl" : lambda e: e['links']['files'],
            }

            header = "#" + "\t".join(columns)
            output(header)
            for e in events:
                row = [ accessFun.get(column, lambda e: defaultAccess(e,column))(e) for column in columns ]
                row = "\t".join(row)
                output(row)
        
        return 0
    elif args[0] == 'replace':
        if len(args) != 3:
            op.error("wrong number of args for replace")
        graceid = args[1]
        filename = args[2]
        response = client.replaceEvent(graceid, filename)
    elif len(args) in [3,4]:
        # Create a new event.
        group = args[0]
        pipeline = args[1]
        if len(args)==3:
            search = None
            filename = args[2]
        else:
            search = args[2]
            filename = args[3]

        # Check that the group, search, and pipeline are known to the API.
        foundGroup = True if (unicode(group) in client.groups) else False
        if not foundGroup:
            error("Group must be one of: ", ", ".join(client.groups))
            sys.exit(1)

        foundPipeline = True if (unicode(pipeline) in client.pipelines) else False
        if not foundPipeline:
            error("Pipeline must be one of: ", ", ".join(client.pipelines))
            sys.exit(1)

        if search:
            foundSearch = True if (unicode(search) in client.searches) else False
            if not foundSearch:
                error("Search must be one of: ", ", ".join(client.searches))
                sys.exit(1)

        response = client.createEvent(group, pipeline, filename, search)
        if not response:
            error("There was a problem.  Did you do ligo-proxy-init?")
            sys.exit(1)

        # XXX Must output graceid for consistency with earlier client.
        # Therefore, must deal with response here rather than at the end.
        exitCode = 0
        status = response.status
        if status >= 400:
            exitCode=1
        try:
            rv = response.read()
        except:
            rv = response
        try:
            rv = json.loads(rv)
        except:
            pass

        if 'graceid' in rv.keys():
            output(rv['graceid'])
        elif 'error' in rv.keys():
            exitCode=1
            error(rv['error'])

        return exitCode

    else:
        op.error("")
        sys.exit(1)

    # Output the response.
    exitCode = 0
    try:
        rv = response.read()
        status = response.status
    except:
        rv = response

    try:
        responseBody = json.loads(rv)
    except:
        responseBody = rv

    if status >= 400:
        exitCode=1
    if isinstance(responseBody, str):
        output("%d: %s" % (status, responseBody))
    else:
        output("Server returned %d" % status)
        if ('error' in responseBody) and response['error']:
            error(response['error'])
            exitCode = 1
        if ('warning' in responseBody) and response['warning']:
            warning(response['warning'])
        if ('output' in responseBody) and response['output']:
            output(response['output'])

    return exitCode

if __name__ == "__main__":
    code = main()
    sys.exit(code)
