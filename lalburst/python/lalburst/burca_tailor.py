# Copyright (C) 2007-2014  Kipp Cannon
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


import copy
import itertools
import math
import sys


import lal
from lal import rate


from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import array as ligolw_array
from glue.ligolw import param as ligolw_param
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw.utils import search_summary as ligolw_search_summary
from . import snglcoinc
from .SimBurstUtils import MW_CENTER_J2000_RA_RAD, MW_CENTER_J2000_DEC_RAD


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from git_version import date as __date__
from git_version import version as __version__


#
# =============================================================================
#
#                Excess Power Specific Parameter Distributions
#
# =============================================================================
#


class LnLRDensity(snglcoinc.LnLRDensity):
	def __init__(self, instruments):
		self.densities = {}
		for pair in intertools.combinations(sorted(instruments), 2):
			# FIXME:  hard-coded for directional search
			#dt = 0.02 + snglcoinc.light_travel_time(*pair)
			dt = 0.02
			self.densities["%s_%s_dt" % pair] = rate.BinnedLnDPF(rate.NDBins((rate.ATanBins(-dt, +dt, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))))
			self.densities["%s_%s_dband" % pair] = rate.BinnedLnDPF(rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))))
			self.densities["%s_%s_ddur" % pair] = rate.BinnedLnDPF(rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))))
			self.densities["%s_%s_df" % pair] = rate.BinnedLnDPF(rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))))
			self.densities["%s_%s_dh" % pair] = rate.BinnedLnDPF(rate.NDBins((rate.LinearBins(-2.0, +2.0, 12001), rate.LinearBins(0.0, 2 * math.pi, 61))))

	def __call__(self, **params):
		try:
			interps = self.interps
		except AttributeError:
			self.mkinterps()
			interps = self.interps
		return sum(interps[param](value) for param, value in params.items())

	def __iadd__(self, other):
		if type(self) != type(other) or set(self.densities) != set(other.densities):
			raise TypeError("cannot add %s and %s" % (type(self), type(other)))
		for key, pdf in self.densities.items():
			pdf += other.densities[key]
		del self.interps
		return self

	def increment(self, params, weight = 1.0):
		for param, value in params.items():
			self.densities[param].count[value] += weight

	def copy(self):
		new = type(self)([])
		for key, pdf in self.densities.items():
			new.densities[key] = pdf.copy()
		return new

	def mkinterps(self):
		self.interps = dict((key, pdf.mkinterp()) for key, pdf in self.densities.items())

	def finish(self):
		for key, pdf in self.densities.items():
			rate.filter_array(pdf.array, rate.gaussian_window(11, 5))
			pdf.normalize()
		self.mkinterps()

	def to_xml(self, name):
		xml = super(LnLRDensity, self).to_xml(name)
		instruments =  set(key.split("_", 2)[0] for key in self.densities if key.endswith("_dt"))
		instruments |= set(key.split("_", 2)[1] for key in self.densities if key.endswith("_dt"))
		xml.appendChild(ligolw_param.Param.from_pyvalue("instruments", lsctables.instrumentsproperty.set(instruments)))
		for key, pdf in self.densities.items():
			xml.appendChild(pdf.to_xml(key))
		return xml

	@classmethod
	def from_xml(cls, xml, name):
		xml = cls.get_xml_root(xml, name)
		self = cls(lsctables.instrumentsproperty.get(ligolw_param.get_pyvalue(xml, "instruments")))
		for key in self.densities:
			self.densities[key] = rate.BinnedLnPDF.from_xml(xml, key)
		return self


class BurcaCoincParamsDistributions(snglcoinc.LnLikelihoodRatioMixin):
	ligo_lw_name_suffix = u"excesspower_coincparamsdistributions"

	@ligolw_array.use_in
	@ligolw_param.use_in
	@lsctables.use_in
	class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
		pass

	def __init__(self, instruments):
		self.numerator = LnLRDensity(instruments)
		self.denominator = LnLRDensity(instruments)
		self.candidates = LnLRDensity(instruments)

	def __iadd__(self, other):
		if type(self) != type(other):
			raise TypeError(other)
		self.numerator += other.numerator
		self.denominator += other.denominator
		self.candidates += other.candidates
		return self

	def copy(self):
		new = type(self)([])
		new.numerator = self.numerator.copy()
		new.denominator = self.denominator.copy()
		new.candidates = self.candidates.copy()
		return new

	def finish(self):
		self.numerator.finish()
		self.denominator.finish()
		self.candidates.finish()

	@classmethod
	def get_xml_root(cls, xml, name):
		name = u"%s:%s" % (name, cls.ligo_lw_name_suffix)
		xml = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.Name == name]
		if len(xml) != 1:
			raise ValueError("XML tree must contain exactly one %s element named %s" % (ligolw.LIGO_LW.tagName, name))
		return xml[0]

	@classmethod
	def from_xml(cls, xml, name):
		xml = cls.get_xml_root(xml, name)
		self = cls([])
		self.numerator = LnLRDensity.from_xml(xml, "numerator")
		self.denominator = LnLRDensity.from_xml(xml, "denominator")
		self.candidates = LnLRDensity.from_xml(xml, "candidates")

	def to_xml(self, name):
		xml = ligolw.LIGO_LW({u"Name": u"%s:%s" % (name, self.ligo_lw_name_suffix)})
		xml.appendChild(self.numerator.to_xml("numerator"))
		xml.appendChild(self.denominator.to_xml("denominator"))
		xml.appendChild(self.candidates.to_xml("candidates"))
		return xml

	@classmethod
	def from_filenames(cls, filenames, name, verbose = False):
		"""
		Convenience function to deserialize
		CoincParamsDistributions objects from a collection of XML
		files and return their sum.  The return value is a
		two-element tuple.  The first element is the deserialized
		and summed CoincParamsDistributions object, the second is a
		segmentlistdict indicating the interval of time spanned by
		the out segments in the search_summary rows matching the
		process IDs that were attached to the
		CoincParamsDistributions objects in the XML.
		"""
		self = None
		for n, filename in enumerate(filenames, 1):
			if verbose:
				print >>sys.stderr, "%d/%d:" % (n, len(filenames)),
			xmldoc = ligolw_utils.load_filename(filename, verbose = verbose, contenthandler = cls.contenthandler)
			if self is None:
				self = cls.from_xml(xmldoc, name)
				seglists = lsctables.SearchSummaryTable.get_table(xmldoc).get_out_segmentlistdict(set([self.process_id])).coalesce()
			else:
				other = cls.from_xml(xmldoc, name)
				self += other
				seglists |= lsctables.SearchSummaryTable.get_table(xmldoc).get_out_segmentlistdict(set([other.process_id])).coalesce()
				del other
			xmldoc.unlink()
		return self, seglists


#
# All sky version
#


class EPAllSkyCoincParamsDistributions(BurcaCoincParamsDistributions):
	def ln_lr_from_triggers(self, events, offsetvector):
		#
		# check for coincs that have been vetoed entirely
		#

		if len(events) < 2:
			return None

		params = {}

		# the "time" is the ms_snr squared weighted average of the
		# peak times neglecting light-travel times.  because
		# LIGOTimeGPS objects have overflow problems in this sort
		# of a calculation, the first event's peak time is used as
		# an epoch and the calculations are done w.r.t. that time.

		# FIXME: this time is available as the peak_time in the
		# multi_burst table, and it should be retrieved from that
		# table instead of being recomputed
		events = tuple(events)
		t = events[0].peak
		t += sum(float(event.peak - t) * event.ms_snr**2.0 for event in events) / sum(event.ms_snr**2.0 for event in events)
		gmst = lal.GreenwichMeanSiderealTime(t) % (2 * math.pi)

		for event1, event2 in itertools.combinations(sorted(events, key = lambda x: x.ifo), 2):
			if event1.ifo == event2.ifo:
				# a coincidence is parameterized only by
				# inter-instrument deltas
				continue

			prefix = "%s_%s_" % (event1.ifo, event2.ifo)

			# in each of the following, if the list of events contains
			# more than one event from a given instrument, the smallest
			# deltas are recorded

			dt = float(event1.peak + offsetvector[event1.ifo] - event2.peak - offsetvector[event2.ifo])
			name = "%sdt" % prefix
			if name not in params or abs(params[name][0]) > abs(dt):
				#params[name] = (dt,)
				params[name] = (dt, gmst)

			df = (event1.peak_frequency - event2.peak_frequency) / ((event1.peak_frequency + event2.peak_frequency) / 2)
			name = "%sdf" % prefix
			if name not in params or abs(params[name][0]) > abs(df):
				#params[name] = (df,)
				params[name] = (df, gmst)

			dh = (event1.ms_hrss - event2.ms_hrss) / ((event1.ms_hrss + event2.ms_hrss) / 2)
			name = "%sdh" % prefix
			if name not in params or abs(params[name][0]) > abs(dh):
				#params[name] = (dh,)
				params[name] = (dh, gmst)

			dband = (event1.ms_bandwidth - event2.ms_bandwidth) / ((event1.ms_bandwidth + event2.ms_bandwidth) / 2)
			name = "%sdband" % prefix
			if name not in params or abs(params[name][0]) > abs(dband):
				#params[name] = (dband,)
				params[name] = (dband, gmst)

			ddur = (event1.ms_duration - event2.ms_duration) / ((event1.ms_duration + event2.ms_duration) / 2)
			name = "%sddur" % prefix
			if name not in params or abs(params[name][0]) > abs(ddur):
				#params[name] = (ddur,)
				params[name] = (ddur, gmst)

		return self(**params)


#
# Galactic core coinc params
#


def delay_and_amplitude_correct(event, ra, dec):
	# retrieve station metadata

	detector = lal.cached_detector_by_prefix[event.ifo]

	# delay-correct the event to the geocentre

	delay = lal.TimeDelayFromEarthCenter(detector.location, ra, dec, event.peak)
	event.peak -= delay
	event.period = event.period.shift(-delay)
	try:
		event.ms_peak -= delay
	except AttributeError:
		pass
	try:
		event.ms_period = event.ms_period.shift(-delay)
	except AttributeError:
		pass

	# amplitude-correct the event using the polarization-averaged
	# antenna response

	fp, fc = lal.ComputeDetAMResponse(detector.response, ra, dec, 0, lal.GreenwichMeanSiderealTime(event.peak))
	mean_response = math.sqrt(fp**2 + fc**2)
	event.amplitude /= mean_response
	event.ms_hrss /= mean_response

	# done

	return event


class EPGalacticCoreCoincParamsDistributions(BurcaCoincParamsDistributions):
	def ln_lr_from_triggers(self, events, offsetvector):
		ra, dec = MW_CENTER_J2000_RA_RAD, MW_CENTER_J2000_DEC_RAD
		return EPAllSkyCoincParamsDistributions.coinc_params([delay_and_amplitude_correct(copy.copy(event), ra, dec) for event in events], offsetvector)


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
#


process_program_name = "lalapps_burca_tailor"


def gen_likelihood_control(coinc_params_distributions, seglists, name = u"lalapps_burca_tailor", comment = u""):
	xmldoc = ligolw.Document()
	node = xmldoc.appendChild(ligolw.LIGO_LW())

	process = ligolw_process.register_to_xmldoc(xmldoc, program = process_program_name, paramdict = {}, version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = comment)
	coinc_params_distributions.process_id = process.process_id
	ligolw_search_summary.append_search_summary(xmldoc, process, ifos = seglists.keys(), inseg = seglists.extent_all(), outseg = seglists.extent_all())

	node.appendChild(coinc_params_distributions.to_xml(name))

	ligolw_process.set_process_end_time(process)

	return xmldoc
