#!/usr/bin/env @PYTHONPROG@

"""
Read a text file containing a list of SFT logical filenames (LFNs).
For each LFN query the local Local Replica Catalog (LRC) to find
if possible a "local" URL. If a local URL is found add it to a list.
If a local URL is not found query the Replica Location Indext (RLI) 
to find a remote site that knows about the file. Then query the remote 
LRC to obtain a URL good for transfer and move the file. After the
file is moved register its new location in the local LRC and add
it to the list.

The input list of LFNs should have each LFN on a seperate line.
It is expected that the input will be the output from QueryMetadataLFN.
"""

versionN = 1.0

import os
import sys
import getopt
import rlsClient
import urlparse
import popen2
import time
import re
import exceptions
import LSCGridFTP
import socket

class QueryException(exceptions.Exception):
        """
        Used to raise exceptions in this script.
        """
        def __init__(self, args=None):
                self.args = args

def queryLRCforURL(lfn, myRLS):
        """
        Query the LRC for PFNs or URLs for a LFN. Return a list of PFNs.
        """


        try:
                exists = myRLS.lrc_exists(lfn, rlsClient.obj_lrc_lfn)
        except rlsClient.RlsClientException, e:
                msg = "Unable to test LFN %s for existence: %s" % (lfn, e)
                raise QueryException, msg

        if  exists:
                try:
                        pfnList = myRLS.lrc_get_pfn(lfn)
                except rlsClient.RlsClientException, e:
                        msg = "Unable to query LRC with LFN % for PFNs: %s" % (lfn, e)
                        raise QueryException, msg
        else:
                pfnList = []
        
        del myRLS
        return pfnList

def queryLRClocalURL(pfn, myRLS):
        """
        Query the LRC to determine if a PFN is a local URL. 

        Return True or False.
        """

        try:
                attrList = myRLS.lrc_attr_value_get(pfn, rlsClient.obj_lrc_pfn, "localURL")
        except rlsClient.RlsClientException, e:
                #msg = "Unable to query for local URL attribute of PFN %s: %s" % (pfn, e)
                #raise QueryException, msg
                #
                # not critical so return False
                return False

        localURL = attrList[0].get_val()
        
        if localURL == 1: localURL = True
        else: localURL = False

        return localURL


def findLFN(lfn, myRLS):
        """
        Query RLI to find remote LRC and then query remote
        LRC to find a suitable URL to use to obtain the file.

        Returns a URL as a string.
        """

        try:
                exists = myRLS.rli_lfn_exists(lfn)
        except rlsClient.RlsClientException, e:
                msg = "Unable to query RLI for existence of LFN %s: %s" % (lfn, e)
                raise QueryException, msg

        if not exists:
                msg = "Cannot find LFN %s in the RLI" % lfn
                raise QueryException, msg
                

        try:
                lrcList = myRLS.rli_get_lrc(lfn)
        except rlsClient.RlsClientException, e:
                msg = "Unable to query RLI for LRCs with LFN %s: %s" % (lfn, e)
                raise QueryException, msg

        url = None

        for lrc in lrcList:
                try:
                        thisRLS = rlsClient.RlsClient(lrc)
                except rlsClient.RlsClientException, e:
                        # cannot connect to this RLS so go on to the next
                        continue

                try:
                        exists = thisRLS.lrc_exists(lfn, rlsClient.obj_lrc_lfn)
                except rlsClient.RlsClientException, e:
                        # this remote LRC doesn't know about the file
                        # so go onto the next
                        continue

                if exists:
                        try:
                               pfnList = queryLRCforURL(lfn, thisRLS)
                        except Exception, e:
                                # some problen querying with this remote RLS 
                                # so go onto the next
                                continue

                        for pfn in pfnList:
                                try:
                                        local = queryLRClocalURL(pfn, thisRLS)

                                except Exception, e:
                                        # some problem querying remote RLS so go
                                        # onto the next one
                                        continue

                                if not local:
                                        # we have found a URL so disconnect and break
                                        del thisRLS
                                        url = pfn
                                        break
                                
                # finished with this RLS so disconnect
                try:
                        if thisRLS: del thisRLS
                except:
                        pass

        if not url:
                msg = "Unable to find URL for LFN %s" % lfn
                raise QueryException, msg

        return url
        


def movePublish(lfn, myRLS):
        """
        Query remote catalogs to find URL for file, move it to
        the local filesystem using GridFTP, then publish the
        local URL in the local catalog as well as a GridFTP
        URL that other sites can use.
        """
        global bucketDir

        if bucketDir: 
                transferDir = os.path.abspath(bucketDir)
        else:
                transferDir = os.getcwd()

        hostname = socket.getfqdn()

        srcURL = findLFN(lfn, myRLS)
        destURL = "file://localhost%s/%s" % (transferDir, lfn)
        gsiftpURL = "gsiftp://%s%s/%s" % (hostname, transferDir, lfn)

        print srcURL, destURL

        try:
                myMover = LSCGridFTP.SimpleGet(verbose=True)
                myMover(srcURL,destURL)
        except Exception, e:
                msg = "Unable to get LFN %s from %s to %s: %s" % (lfn, srcURL, destURL, e)
                raise QueryException, msg

        try:
                if myRLS.lrc_exists(lfn, rlsClient.obj_lrc_lfn):
                        myRLS.lrc_add(lfn, destURL)
                else:
                        myRLS.lrc_create_lfn(lfn, destURL)
        except rlsClient.RlsClientException, e:
                msg = "Unable to publish LFN %s into LRC with URL %s: %s" % (lfn, destURL, e)
                raise QueryException, msg

        try:
                myRLS.lrc_add(lfn, gsiftpURL)
        except rlsClient.RlsClientException, e:
                # failure to publish the gsiftp is not critical
                pass

        try:
                myAttr = rlsClient.RlsAttr()
                myAttr.set_name("localURL")
                myAttr.set_objtype(rlsClient.obj_lrc_pfn)
                myAttr.set_type(rlsClient.attr_type_int)
                myAttr.set_val(1)

                myRLS.lrc_attr_add(destURL, myAttr)

        except rlsClient.RlsClientException, e:
                # the attribute may not be defined so try defining 
                # and then setting again
                try:
                        myRLS.lrc_attr_define(myAttr)
                        myRLS.lrc_attr_add(destURL, myAttr)
                except rlsClient.RlsClientException, e2:
                        msg = "Unable to set localURL attribute on PFN %s: %s %s" % (destURL, e, e2)
                        raise QueryException, msg

        try:
                myAttr.set_val(0)
                myRLS.lrc_attr_add(gsiftpURL, myAttr)
                
        except:
                # not critical failure 
                pass

        return destURL

def verifyLocalPath(pfn):
        """
        Verify that we can read the file pointed to by the local
        URL.
        """
        path = urlparse.urlparse(pfn)[2]

        return os.access(path, os.R_OK)

def unpublishLocal(pfn, myRLS):
        """
        Delete a LFN -> PFN mapping in the catalog.
        """

        try:
                lfn = myRLS.lrc_get_lfn(pfn)
                myRLS.lrc_delete(lfn, pfn)
        except rlsClient.RlsClientException, e:
                # not a critical failure
                pass
                

def usage():
        """
        Print a usage message to stderr.
        """

        msg = """\
NAME
        QueryRLS

SYNOPSIS
        QueryRLS --input=PATH --server=URL [ --bucket=DIRECTORY ]
                [ --output=PATH ] [ --bucket=PATH ]

        QueryRLS --version

        QueryRLS --help

DESCRIPTION
        Query the RLS server to find the LFNs listed in the 
        input file. If a file is not available locally then
        query remote RLS servers to find the file and then
        use GridFTP to move the file. When completed write
        out the list of local paths to the LFNs.

        -i, --input
                path to the file containing LFNs

        -b, --bucket
                directory into which files should be moved if 
                necessary, defaults to current working directory

        -s, --server
                URL for the RLS server

        -o, --output
                file into which to write the list of paths
                for the LFNs

        -V, --version
                print version number of this script to stderr and exit

        -h, --help
                print this usage message

EXAMPLE

$ QueryRLS --input=LFNlist --server=rls://hydra.phys.uwm.edu 
        --output=LFNpaths

"""

        print >>sys.stderr, msg


longopt = [
        "input=",
        "server=",
        "output=",
        "bucket=",
        "version",
        "help"
        ]

shortopt = "i:s:o:vhb:"


try:
        opts, args = getopt.getopt(sys.argv[1:], shortopt, longopt)
except getopt.GetoptError:
        print >>sys.stderr, "Error parsing command line"
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

# set default values
inputPath = None
rlsURL = None
outputPath = None
bucketDir = None


for o, a in opts:
        if o in ("-h", "--help"):
                usage()
                sys.exit(0)
        elif o in ("-i", "--input"):
                inputPath = a
        elif o in ("-b", "--bucket"):
                bucketDir = a
        elif o in ("-s", "--server"):
                rlsURL = a
        elif o in ("-o", "--output"):
                outputPath = a
        elif o in ("-V", "--version"):
                print >>sys.stderr, versionN
                sys.exit(0)

         
# sanity checking on inputs
if not inputPath:
        print >>sys.stderr, "An input file must be specified"
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

inputPath = os.path.abspath(inputPath)
if not os.access(inputPath, os.R_OK):
        print >>sys.stderr, "Unable to read file %s" % inputPath
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

if bucketDir:
        bucketDir = os.path.abspath(bucketDir)
        if not os.path.isdir(bucketDir):
                print >>sys.stderr, "Unable to access directory %s" % bucketDir
                print >>sys.stderr, "Use --help for usage"

if not rlsURL:
        print >>sys.stderr, "A RLS server URL must be specified"
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

urlparse.uses_netloc.insert(0, 'rls')

try:
        hostname = urlparse.urlparse(rlsURL)[1]
except:
        print >>sys.stderr, "Unable to parse RLS URL"
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

if outputPath:
        outputPath = os.path.abspath(outputPath)
        if os.path.isfile(outputPath):
                if not os.access(outputPath, os.W_OK):
                        print >>sys.stderr, "Cannot write into file %s" % outputPath
                        print >>sys.stderr, "Use --help for usage"
                        sys.exit(1)
        else:
                dir = os.path.dirname(outputPath)
                if not os.access(dir, os.W_OK):
                        print >>sys.stderr, "Cannot write into directory %s" % dir
                        print >>sys.stderr, "Use --help for usage"
                        sys.exit(1)

# make sure that GLOBUS_LOCATION is in our environment

try:
        GLOBUS_LOCATION = os.environ["GLOBUS_LOCATION"]
except:
        print >>sys.stderr, "GLOBUS_LOCATION not defined in environment"
        sys.exit(1)


# read in the list of LFNs
try:
        f = open(inputPath, "r")
        lfnList = [ s.strip() for s in f.readlines() ]
        f.close()
        print "Found %d LFNs in file %s" % (len(lfnList), inputPath)
except Exception, e:
        print >>sys.stderr, "Error reading LFNs from %s" % inputPath
        print >>sys.stderr, "Use --help for usage"
        sys.exit(1)

# query LRC for a local URL if there is one, and create a list
# of paths for LFNs with local URLs and a list of LFNs that
# have to be moved


try:
        myRLS = rlsClient.RlsClient(rlsURL)
        print "Connected to local RLS server %s" % rlsURL
except rlsClient.RlsClientException, e:
        print >>sys.stderr, "Unable to connect to RLS: %s" % e
        sys.exit(1)

localPathList = []
lfnMoveList = []

for lfn in lfnList:
        local = False
        pfnList = queryLRCforURL(lfn, myRLS)

        for pfn in pfnList:
                if queryLRClocalURL(pfn, myRLS):
                        print "LFN %s has local URL %s" % (lfn, pfn)
                        if verifyLocalPath(pfn):
                                print "file at local URL %s is readable" % pfn
                                localPathList.append(pfn)
                                local = True
                        else:
                                unpublishLocal(pfn, myRLS)
                else:
                        print "LFN %s has NON-local URL %s" % (lfn, pfn)

        if not local:
                lfnMoveList.append(lfn)

for lfn in lfnMoveList:
        print "Finding, moving, and publishing %s" % lfn
        localurl = movePublish(lfn, myRLS)
        print "LFN %s now at %s" % (lfn, localurl)
        localPathList.append(localurl)
                
# close connection to RLS 
try:
        del myRLS        
except Exception, e:
        pass
                
if not outputPath:
        f = sys.stdout
else:
        path = os.path.abspath(outputPath)
        try:
                f = open(path, "w")
        except Exception, e:
                print >>sys.stderr, "Unable to open file %s for writing: %s" % (path, e)
                sys.exit(1)

for p in localPathList:
        print >>f, p

if f != sys.stdout: f.close()


# exit nicely
sys.exit(0)

