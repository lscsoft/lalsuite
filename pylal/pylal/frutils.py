# Copyright (C) 2009-11  Nickolas Fotopoulos
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

from __future__ import division

__author__ = "Nickolas Fotopoulos <nickolas.fotopoulos@ligo.org>"

from bisect import bisect_right
import httplib
import os
import os.path
import shutil
import sys
import operator

from glue.lal import Cache
from glue.segments import segment, segmentlist
from pylal.metaarray import TimeSeries, TimeSeriesList
from pylal.Fr import frgetvect1d

__all__ = ('__author__', 'FrameCache', "AutoqueryingFrameCache")

class FrameCache(object):
    """
FrameCache is a transparent interface to LSC data. The user provides a LAL-
formatted cache file and the returned FrameCache object allows repeated
queries for channels and time, even across frame files. It also supports
smart, lazy local caching. Limitations: It only works for one-dimensional
time-series data.

Constructor:
    FrameCache(cache_entries=None, scratchdir=None, verbose=False)

Inputs:
    cache is a list of glue.lal.CacheEntry objects or a glue.lal.Cache.
    Data will be retrieved from the frame files described within.

    Scratchdir determines where to locally cache frames.  If None, no
    caching is performed.

Example:
>>> from glue import lal
>>> from pylal import frutils
>>> c = lal.Cache.fromfile(open("test.cache"))
>>> d = frutils.FrameCache(c, scratchdir="/tmp", verbose=True)
>>> data = d.fetch("H1:LSC-STRAIN", 861417967, 861417969)
Copying /Users/nvf/temp/H-H1_RDS_C03_L2-861417967-128.gwf -->
          /tmp/H-H1_RDS_C03_L2-861417967-128.gwf.
>>> print data
[  1.68448009e-16   1.69713183e-16   1.71046196e-16 ...,   1.80974629e-16
   1.80911765e-16   1.80804879e-16] {'dt': 6.103515625e-05, 'segments': [segment(861417967, 861417969)], 'comments': [], 'name': 'H1:LSC-STRAIN'}
>>> exit()
Removing /tmp/H-H1_RDS_C03_L2-861417967-128.gwf.

"""

    def __init__(self, cache_entries=None, scratchdir=None, verbose=False):
        """ Initializes interface to frame data.  See .__class__.__doc__"""

        # Simple initializations
        # Use list of segments vs segmentlist to prevent merging.
        self._verbose = verbose
        self._scratchdir = scratchdir
        self._remotefiles = []               # filename list
        self._remotesegs = segmentlist()     # list of segments
        self._remotecoverage = segmentlist() # coalesced copy of remotesegs

        # if we have a scratchdir, maintain independent lists
        if scratchdir is not None:
            self._cachedfiles = []
            self._cachedsegs = segmentlist()
            self._cachecoverage = segmentlist()
        else:
            self._cachedfiles = self._remotefiles
            self._cachedsegs = self._remotesegs
            self._cachecoverage = self._remotecoverage

        if cache_entries is not None:
          self.add_cache(cache_entries)

    def add_cache(self, cache_entries):
        """
        Add information from some cache entries.
        """
        newentries = [entry for entry in cache_entries \
                    if entry.path not in self._remotefiles]
        newfiles = [entry.path for entry in newentries]

        # We iterate here to prevent the segment and files from getting added in
        # an unsorted manner
        for entry in cache_entries:
            # We already have it, skip it
            if entry.path in self._remotefiles:
                continue
            newseg, newfile = entry.segment, entry.path
            insert_idx = bisect_right(self._remotesegs, newseg)
            self._remotesegs.insert(insert_idx, newseg)
            self._remotefiles.insert(insert_idx, newfile)
            self._remotecoverage |= segmentlist([newseg])

        self._remotecoverage.coalesce()

    def __del__(self):
        """
        Clear cache in local scratch.
        """
        if self._scratchdir is None:
            return
        for f,s in zip(self._cachedfiles, self._cachedsegs):
            self._unfetch(f, s)
        return

    def fetch(self, channel, start, end):
        """
        Retrieve data, caching file locations and the files themselves.
        """
        seg = segment(start, end)

        if not self._query(channel, start, end):
            raise ValueError("%s not found in cache" % repr(segmentlist([seg]) - self._remotecoverage))

        # Need to cache files locally
        # Note: seg *will* be in self._cachecoverage if self.scratchdir is None.
        if seg not in self._cachecoverage:
            for f,s in zip(self._remotefiles, self._remotesegs):
                if seg.intersects(s) and s not in self._cachecoverage:
                    dest = os.path.join(self._scratchdir, os.path.split(f)[-1])
                    if self._verbose:
                        print "Copying %s -->\n          %s." % (f, dest)
                    shutil.copy(f, dest)
                    ind = bisect_right(self._cachedsegs, s)
                    self._cachedfiles.insert(ind, dest)
                    self._cachedsegs.insert(ind, s)
                    self._cachecoverage |= segmentlist([s])
            assert seg in self._cachecoverage

        # Finally, return the cached data
        return self._fetch(channel, start, end)

    def _query(self, channel, start, end):
        "Do we know where the frame file is?"
        return segment(start, end) in self._remotecoverage

    def _fetch(self, channel, start, end, comments=[]):
        """
        Internal method to actually retrieve and return data as TimeSeries,
        assuming that self._framefiles is all set.  Does not check boundaries.
        """
        toreturn = TimeSeriesList([])

        if start==end:
            return toreturn

        # Find first frame
        try:
            #tmp = sorted(zip(self._cachedsegs, self._cachedfiles))
            #self._cachedsegs = segmentlist([seg for seg, frfile in tmp])
            #self._cachedfiles = [frfile for seg, frfile in tmp]
            index = self._cachedsegs.find(start)
        except ValueError:
            print >>sys.stderr, "Couldn't find any frame files to cover",\
                str(start),"to",str(end),"among:"
            print >>sys.stderr, str(self._cachedfiles)
            return toreturn

        # Get frames; an error probably means that the frames didn't cover
        # the whole period of time.  Cleanly handles frames of varying lengths.
        now = start
        while now < end:
            dur = min(end, self._cachedsegs[index][1]) - now
            data, GPS_start, t_low, dt, x_unit, y_unit = \
                frgetvect1d(self._cachedfiles[index], channel, now, dur, 0)
            meta = {"name": channel, "dt": dt,
                "segments": [segment(now, now+dur)], "comments": comments}
            toreturn.append(TimeSeries(data, meta))
            now += dur
            index += 1

        if len(toreturn) == 0:
            print >>sys.stderr, "This print statement should never execute."
            print >>sys.stderr,"Couldn't find all frame files needed to cover",\
                str(start), "to", str(end), "among:"
            print >>sys.stderr, str(self._cachedfiles)

        toreturn = toreturn.merge_list()
        toreturn.metadata.segments.coalesce()

        return toreturn

    def unfetch(self, start, end):
        """
        Removes files from local scratch space based on start, end
        pairs.  Silently ignores non-existent times.  Remove if file end
        is between start and end.  This is biased to prevent cache misses
        for future fetches being in the future.  (Processing frames in
        chronological order)
        """
        if self._scratchdir is None:
            return

        for f,s in zip(self._cachedfiles, self._cachedsegs):
            if start < s[1] <= end:
                self._unfetch(f,s)

    def _unfetch(self, filename, seg):
        """
        Internal method to actually remove a file from cache.
        """
        if self._scratchdir is None:
            return
        if filename not in self._cachedfiles:
            print >>sys.stderr, \
                "Cache inconsistency: Delete request for file not in cache."
            return
        if self._verbose: print "Removing %s." % filename
        os.remove(filename)
        self._cachedfiles.remove(filename)
        self._cachedsegs.remove(seg)
        self._cachecoverage -= segmentlist([seg])
        return

#
# Set up a FrameCache subclass that queries LDR on-the-fly for frame locations;
# Contains many bits stolen from ligo_data_find in Glue.
#

def validateProxy(path):
    """
    Test that the proxy certificate is RFC 3820
    compliant and that it is valid for at least
    the next 15 minutes.
    """
    try:
        import M2Crypto
    except ImportError, e:
        print >> sys.stderr, """
validateProxy requires the M2Crypto module.

On CentOS 5 and other RHEL-based platforms
this package is available from the EPEL
repository by doing

yum install m2crypto

For Debian Lenny this package is available
by doing

apt-get install python-m2crypto

Mac OS X users can find this package in MacPorts.

%s
""" % e
        raise

    # load the proxy from path
    try:
        proxy = M2Crypto.X509.load_cert(path)
    except Exception, e:
        msg = "Unable to load proxy from path %s : %s" % (path, e)
        raise RuntimeError(msg)

    # make sure the proxy is RFC 3820 compliant
    try:
        proxy.get_ext("proxyCertInfo")
    except LookupError:
        rfc_proxy_msg = """\
Could not find a RFC 3820 compliant proxy credential.
Please run 'grid-proxy-init -rfc' and try again.
"""
        raise RuntimeError(rfc_proxy_msg)

    # attempt to make sure the proxy is still good for more than 15 minutes
    import time, calendar
    try:
        expireASN1 = proxy.get_not_after().__str__()
        expireGMT  = time.strptime(expireASN1, "%b %d %H:%M:%S %Y %Z")
        expireUTC  = calendar.timegm(expireGMT)
        now = int(time.time())
        secondsLeft = expireUTC - now
    except Exception, e:
        # problem getting or parsing time so just let the client
        # continue and pass the issue along to the server
        secondsLeft = 3600

    if secondsLeft <= 0:
        msg = """\
Your proxy certificate is expired.

Please generate a new proxy certificate and
try again.
"""
        raise RuntimeError(msg)

    if secondsLeft < (60 * 15):
        msg = """\
Your proxy certificate expires in less than
15 minutes.

Please generate a new proxy certificate and
try again.
"""
        raise RuntimeError(msg)

    # return True to indicate validated proxy
    return True

def findCredential():
    """
    Follow the usual path that GSI libraries would
    follow to find a valid proxy credential but
    also allow an end entity certificate to be used
    along with an unencrypted private key if they
    are pointed to by X509_USER_CERT and X509_USER_KEY
    since we expect this will be the output from
    the eventual ligo-login wrapper around
    kinit and then myproxy-login.
    """
    rfc_proxy_msg = """\
Could not find a RFC 3820 compliant proxy credential.
Please run 'grid-proxy-init -rfc' and try again.
"""

    # use X509_USER_PROXY from environment if set
    if os.environ.has_key('X509_USER_PROXY'):
        filePath = os.environ['X509_USER_PROXY']
        if validateProxy(filePath):
            return filePath, filePath
        else:
            raise RuntimeError(rfc_proxy_msg)

    # use X509_USER_CERT and X509_USER_KEY if set
    if os.environ.has_key('X509_USER_CERT'):
        if os.environ.has_key('X509_USER_KEY'):
            certFile = os.environ['X509_USER_CERT']
            keyFile = os.environ['X509_USER_KEY']
            return certFile, keyFile

    # search for proxy file on disk
    uid = os.getuid()
    path = "/tmp/x509up_u%d" % uid

    if os.access(path, os.R_OK):
        if validateProxy(path):
            return path, path
        else:
            raise RuntimeError(rfc_proxy_msg)

    # if we get here could not find a credential
    raise RuntimeError(rfc_proxy_msg)

def query_LDR(server, port, site, frameType, gpsStart, gpsEnd, urlType=None, noproxy=False):
    """
    Return a list of URLs to frames covering the requested time, as returned
    by the LDR server.
    """
    try:
        import cjson
    except ImportError, e:
        print >> sys.stderr, """
    frutils requires the cjson module.

    On CentOS 5 and other RHEL-based platforms
    this package is available from the EPEL
    repository by doing

    yum install python-cjson

    For Debian Lenny this package is available by doing

    apt-get install python-cjson

    Mac OS X users can find this package in MacPorts.

    %s
    """ % e
        raise

    url = "/LDR/services/data/v1/gwf/%s/%s/%s,%s" % (site, frameType, gpsStart, gpsEnd)
    # if a URL type is specified append it to the path
    if urlType:
        url += "/%s" % urlType

    # request JSON output
    url += ".json"

    # make unauthenticated request
    if noproxy or port == 80:
        h = httplib.HTTPConnection(server, port)
    else:
        certFile, keyFile = findCredential()
        h = httplib.HTTPSConnection(server, key_file = keyFile, cert_file = certFile)

    # query the server
    try:
        h.request("GET", url)
        response = h.getresponse()
    except Exception, e:
        msg = "Unable to query server %s: %s\n\nPerhaps you need a valid proxy credential?\n" % (server, e)
        raise RuntimeError(msg)

    # the server did respond to check the status
    if response.status != 200:
        msg = "Server returned code %d: %s" % (response.status, response.reason)
        body = response.read()
        msg += body
        raise RuntimeError(msg)

    # since status is 200 OK read the URLs
    body = response.read()

    # decode the JSON
    return cjson.decode(body)

class AutoqueryingFrameCache(FrameCache):
    """
This subclass of FrameCache will query ligo_data_find automatically,
so no LAL-cache files are required. Limitation: you'll need one instance
per frame type.

Constructor:
    AutoqueryingFrameCache(frametype, hostPortString=None, scratchdir=None,
        verbose=False)

Inputs:
    frametype is the type of GWF frame you seek (e.g. RDS_R_L1).
    hostPortString is the name of the LDR server and optionally,
        with colon separation, the port (e.g. ldr.ligo.caltech.edu)
    scratchdir determines where to locally cache frames. If None, no
        caching is performed.

Example:
>>> from pylal import frutils
>>> d = frutils.AutoqueryingFrameCache(frametype="H1_RDS_C03_L2", scratchdir="/tmp", verbose=True)
>>> data = d.fetch("H1:LSC-STRAIN", 861417967, 861417969)
Copying /Users/nvf/temp/H-H1_RDS_C03_L2-861417967-128.gwf -->
          /tmp/H-H1_RDS_C03_L2-861417967-128.gwf.
>>> print data
[  1.68448009e-16   1.69713183e-16   1.71046196e-16 ...,   1.80974629e-16
   1.80911765e-16   1.80804879e-16] {'dt': 6.103515625e-05, 'segments': [segment(861417967, 861417969)], 'comments': [], 'name': 'H1:LSC-STRAIN'}
>>> exit()
Removing /tmp/H-H1_RDS_C03_L2-861417967-128.gwf.

Using AutoqueryingFrameCache outside of LDG clusters, using Caltech as a
gateway:
 * Just the first time you do this procedure: "sudo mkdir /data && sudo chown
   albert.einstein /data" (replace albert.einstein with your local username;
   /data may be different for different clusters)
 * Set the LIGO_DATAFIND_SERVER environment variable to ldr.ligo.caltech.edu
   (or the LDR server of the LDG cluster nearest you)
 * Use "sshfs -o ssh_command=gsissh
   albert.einstein@ldas-pcdev1.ligo.caltech.edu:/data /data" (replace
   albert.einstein with your cluster username)
 * Use "umount /data" when you're done. Unmounting cleanly will help prevent
   headaches the next time you want to set this up.
    """
    def __init__(self, frametype, hostPortString=None, scratchdir=None,
        verbose=False):
        FrameCache.__init__(self, None, scratchdir, verbose)

        if not frametype:
            raise ValueError("frametype required")
        self.frametype = frametype

        if hostPortString is None:
            if os.environ.has_key('LIGO_DATAFIND_SERVER'):
                hostPortString = os.environ['LIGO_DATAFIND_SERVER']
            else:
                raise ValueError("no way to determine LIGO_DATAFIND_SERVER")
        if hostPortString.find(':') < 0:
            # no port specified
            self.host = hostPortString
            self.port = None
        else:
            # server and port specified
            self.host, portString = hostPortString.split(':')
            self.port = int(portString)

    def _query(self, channel, start, end):
        "Do we know where the frame file is?"
        if segment(start, end) in self._remotecoverage:
            return True
        urls = query_LDR(self.host, self.port, channel[0], self.frametype, start, end, urlType="file")
        if urls:
            new = Cache.from_urls(urls, coltype=int)
            new.sort(key=operator.attrgetter("segment"))
            self.add_cache(new)
        return segment(start, end) in self._remotecoverage
