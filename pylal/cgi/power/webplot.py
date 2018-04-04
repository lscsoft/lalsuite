import cgi
import cgitb ; cgitb.enable()
import os
import shutil
import sys
import tempfile
import time
import urllib

from glue.lal import CacheEntry
from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw.utils import ligolw_add

from pylal import llwapp
from pylal import SnglBurstUtils
from pylal import ligolw_bucluster
from pylal.date import XLALUTCToGPS, LIGOTimeGPS

import eventdisplay

#
# =============================================================================
#
#                             CGI display request
#
# =============================================================================
#

class PlotDescription(object):
	def __init__(self):
		# set defaults
		now = XLALUTCToGPS(time.gmtime())
		self.segment = segments.segment(now, now + (-1 * 3600))
		self.ratewidth = 60.0
		self.freqwidth = 16.0
		self.band = segments.segment(0.0, 2500.0)
		self.set_instrument("H1")
		self.seglist = segments.segmentlist([self.segment])
		self.filename = None
		self.format = None	# force update
		self.set_format("png")
		self.cluster = 0

		return self

	def __del__(self):
		self.remove_tmpfile()

	def remove_tmpfile(self):
		if self.filename:
			os.remove(self.filename)

	def set_format(self, format):
		if format not in ["eps", "png", "svg", "xml"]:
			raise Exception, "unrecognized format %s" % format
		if self.format == format:
			return
		self.format = format
		if format == "eps":
			self.set_extension("eps")
		elif format == "png":
			self.set_extension("png")
		elif format == "svg":
			self.set_extension("svg")
		elif format == "xml":
			self.set_extension("xml")

	def set_extension(self, extension):
		self.remove_tmpfile()
		handle, self.filename = tempfile.mkstemp("." + extension.strip(), "webplot_")
		os.close(handle)

	def set_instrument(self, instrument):
		# update cache name along with instrument
		self.instrument = str(instrument)

	def parse_form(self):
		# parse CGI form
		form = cgi.FieldStorage()

		start = LIGOTimeGPS(form.getfirst("start", str(self.segment[0])))
		duration = LIGOTimeGPS(form.getfirst("dur", str(abs(self.segment))))

		self.segment = segments.segment(start, start + duration)
		self.ratewidth = float(form.getfirst("ratewidth", str(self.ratewidth)))
		self.freqwidth = float(form.getfirst("freqwidth", str(self.freqwidth)))
		self.band = segments.segment(float(form.getfirst("lofreq", str(self.band[0]))), float(form.getfirst("hifreq", str(self.band[1]))))
		self.set_instrument(form.getfirst("inst", self.instrument))
		self.set_format(form.getfirst("format", self.format))
		self.cluster = int(form.getfirst("cluster", "0"))

		return self

	def trig_segment(self):
		# interval in which triggers must be read in order to
		# produce a plot
		return self.segment


#
# =============================================================================
#
#               How to get a table of triggers within a segment
#
# =============================================================================
#

def element_filter(name, attrs):
	return lsctables.IsTableProperties(lsctables.SnglBurstTable, name, attrs) or lsctables.IsTableProperties(lsctables.SimBurstTable, name, attrs) or lsctables.IsTableProperties(lsctables.SearchSummaryTable, name, attrs)

class ContentHandler(ligolw.PartialLIGOLWContentHandler):
	"""
	ContentHandler for reading only sngl_burst, sim_burst and
	search_summary tables.
	"""
	def __init__(self, doc):
		ligolw.PartialLIGOLWContentHandler.__init__(self, doc, element_filter)

ligolw_add.ContentHandler = ContentHandler


def CacheURLs(cachename, seg):
	"""
	Open a trigger cache, and return a list of URLs for files
	intersecting seg.
	"""
	return [c.url for c in map(CacheEntry, file(cachename)) if c.segment.intersects(seg)]

def gettriggers(plotdesc):
	doc = ligolw_add.ligolw_add(ligolw.Document(), CacheURLs(eventdisplay.cache[plotdesc.instrument], plotdesc.segment), verbose = False, non_lsc_tables_ok = False)
	try:
		plotdesc.seglist = table.get_table(doc, lsctables.SearchSummaryTable.tableName).get_outlist().coalesce()
	except:
		plotdesc.seglist = segments.segmentlist()
	try:
		bursttable = table.get_table(doc, lsctables.SnglBurstTable.tableName)
	except:
		bursttable = lsctables.New(lsctables.SnglBurstTable)
	try:
		simtable = table.get_table(doc, lsctables.SimBurstTable.tableName)
	except:
		simtable = lsctables.New(lsctables.SnglBurstTable)

	# cluster
	if plotdesc.cluster:
		ligolw_bucluster.ClusterSnglBurstTable(bursttable, SnglBurstUtils.CompareSnglBurstByPeakTimeAndFreq, ligolw_bucluster.SnglBurstCluster, SnglBurstUtils.CompareSnglBurstByPeakTime)

	# remove triggers and injections that lie outside the required segment
	bursttable.filterRows(lambda row: row.get_peak() in plotdesc.trig_segment())
	simtable.filterRows(lambda row: row.get_time_geocent() in plotdesc.trig_segment())

	return bursttable, simtable


#
# =============================================================================
#
#                                    Output
#
# =============================================================================
#

def SendImage(plotdesc):
	if plotdesc.format == "png":
		print >>sys.stdout, "Content-Type: image/png\n"
	elif plotdesc.format == "eps":
		print >>sys.stdout, "Content-Type: application/postscript\n"
	elif plotdesc.format == "svg":
		print >>sys.stdout, "Content-Type: image/svg+xml\n"
	elif plotdesc.format == "xml":
		print >>sys.stdout, "Content-Type: text/xml\n"
	shutil.copyfileobj(file(plotdesc.filename), sys.stdout)
