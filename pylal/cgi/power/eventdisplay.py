#!/usr/bin/python

import os

from glue import LSCsegFindClient
from glue import segments
from glue import segmentsUtils
from pylal.date import LIGOTimeGPS

#
# Some info
#

s5start = LIGOTimeGPS(815155213)
cache = {
	"G1": os.getcwd() + "/G1/triggers.cache",
	"G1 Injections": os.getcwd() + "/G1/injections.cache",
	"H1": os.getcwd() + "/H1/triggers.cache",
	"H1 Injections": os.getcwd() + "/H1/injections.cache",
	"H2": os.getcwd() + "/H2/triggers.cache",
	"H2 Injections": os.getcwd() + "/H2/injections.cache",
	"L1": os.getcwd() + "/L1/triggers.cache",
	"L1 Injections": os.getcwd() + "/L1/injections.cache"
}


#
# Trigger file segment lists
#

class TrigSegs(object):
	def __init__(self):
		self.G1 = segmentsUtils.fromlalcache(file(cache["G1"]), coltype = LIGOTimeGPS).coalesce()
		self.H1 = segmentsUtils.fromlalcache(file(cache["H1"]), coltype = LIGOTimeGPS).coalesce()
		self.H2 = segmentsUtils.fromlalcache(file(cache["H2"]), coltype = LIGOTimeGPS).coalesce()
		self.L1 = segmentsUtils.fromlalcache(file(cache["L1"]), coltype = LIGOTimeGPS).coalesce()


#
# Segment querying
#

class SegFindConfig(object):
	def __init__(self, host, port, instrument):
		self.host = host
		self.port = port
		self.instrument = instrument

SegFindConfigH1 = SegFindConfig("ldas.ligo-wa.caltech.edu", None, "H1")
SegFindConfigH2 = SegFindConfig("ldas.ligo-wa.caltech.edu", None, "H2")
SegFindConfigL1 = SegFindConfig("ldas.ligo-la.caltech.edu", None, "L1")

def getsegments(config, types, bounds):
	if config.port:
		client = LSCsegFindClient.LSCsegFind(config.host, config.port)
	else:
		client = LSCsegFindClient.LSCsegFind(config.host)
	list = client.findStateSegments({"interferometer" : config.instrument, "type" : types, "start" : str(int(bounds[0])), "end" : str(int(bounds[1])), "lfns" : False, "strict" : True})
	return segments.segmentlist([segments.segment(*map(LIGOTimeGPS, seg)) for seg in list])

