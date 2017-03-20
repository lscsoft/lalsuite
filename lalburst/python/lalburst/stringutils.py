# Copyright (C) 2009--2014  Kipp Cannon
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


import itertools
import math
import scipy.stats
import sys


import lal


from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import array as ligolw_array
from glue.ligolw import param as ligolw_param
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from glue.offsetvector import offsetvector
from . import git_version
from pylal import rate
from pylal import snglcoinc


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                             Likelihood Machinery
#
# =============================================================================
#


#
# Make a look-up table of time-of-arrival triangulators
#


def triangulators(timing_uncertainties):
	"""
	Return a dictionary of snglcoinc.TOATriangulator objects
	initialized for a variety of instrument combinations.
	timing_uncertainties is a dictionary of instrument->$\\Delta t$
	pairs.  The return value is a dictionary of (instrument
	tuple)->TOATrangulator mappings.  The instrument names in each
	tuple are sorted in alphabetical order, and the triangulators are
	constructed with the instruments in that order (the the
	documentation for snglcoinc.TOATriangulator for more information).

	Example:

	>>> x = triangulators({"H1": 0.005, "L1": 0.005, "V1": 0.005})

	constructs a dictionary of triangulators for every combination of
	two or more instruments that can be constructed from those three.

	The program lalapps_string_plot_binj can be used to measure the
	timing uncertainties for the instruments in a search.
	"""
	allinstruments = sorted(timing_uncertainties.keys())

	triangulators = {}
	for n in range(2, len(allinstruments) + 1):
		for instruments in itertools.combinations(allinstruments, n):
			triangulators[instruments] = snglcoinc.TOATriangulator([lal.cached_detector_by_prefix[instrument].location for instrument in instruments], [timing_uncertainties[instrument] for instrument in instruments])

	return triangulators


#
# Parameter distributions
#


class LnLRDensity(snglcoinc.LnLRDensity):
	# FIXME:  the interps dictionary maintained here should be
	# eliminated in favour of an internal mechanism within the PDFs
	# themselves that performs the interpolation on-the-fly, without
	# requiring an intermediate object to be created
	def __init__(self, instruments):
		self.densities = {}
		for instrument in instruments:
			self.densities["%s_snr2_chi2" % instrument] = rate.BinnedLnPDF(rate.NDBins((rate.ATanLogarithmicBins(10, 1e7, 801), rate.ATanLogarithmicBins(.1, 1e4, 801))))
		for pair in itertools.combinations(sorted(instruments), 2):
			dt = 0.005 + snglcoinc.light_travel_time(instrument1, instrument2)	# seconds
			self.densities["%s_%s_dt" % pair] = rate.BinnedLnPDF(rate.NDBins((rate.ATanBins(-dt, +dt, 801),)))
			self.densities["%s_%s_dA" % pair] = rate.BinnedLnPDF(rate.NDBins((rate.ATanBins(-0.5, +0.5, 801),)))
			self.densities["%s_%s_df" % pair] = rate.BinnedLnPDF(rate.NDBins((rate.ATanBins(-0.2, +0.2, 501),)))
		# only non-negative rss timing residual bins will be used
		# but we want a binning that's linear at the origin so
		# instead of inventing a new one we just use atan bins that
		# are symmetric about 0
		self.densities["instrumentgroup,rss_timing_residual"] = rate.BinnedLnPDF(rate.NDBins((snglcoinc.InstrumentBins(names = instruments), rate.ATanBins(-0.02, +0.02, 1001))))

	def __call__(self, params):
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
			if key.endswith("_snr2_chi2"):
				rate.filter_array(pdf.array, rate.gaussian_window(11, 11, sigma = 20))
			elif key.endswith("_dt") or key.endswith("_dA") or key.endswith("_df"):
				rate.filter_array(pdf.array, rate.gaussian_window(11, sigma = 20))
			elif key.startswith("instrumentgroup"):
				# instrument group filter is a no-op
				pass
			else:
				# shouldn't get here
				raise Exception
			pdf.normalize()
		self.mkinterps()

	def to_xml(self, name):
		xml = super(LnLRDensity, self).to_xml(name)
		instruments = set(key.split("_", 1)[0] for key in self.densities if key.endswith("_snr2_chi2"))
		xml.appendChild(ligolw_param.Param.from_pyvalue("instruments", lsctables.ifos_from_instrument_set(instruments)))
		for key, pdf in self.densities.items():
			xml.appendChild(pdf.to_xml(key))
		return xml

	@classmethod
	def from_xml(cls, name):
		xml = cls.get_xml_root(xml, name)
		self = cls(lsctables.instrument_set_from_ifos(ligolw_param.get_pyvalue(xml, "instruments")))
		for key in self.densities:
			self.densities[key] = rate.BinnedLnPDF.from_xml(xml, key)
		return self


class StringCoincParamsDistributions(snglcoinc.LnLikelihoodRatioMixin):
	ligo_lw_name_suffix = u"stringcusp_coincparamsdistributions"

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

	def copy(self):
		new = type(self)([])
		new.numerator = self.numerator.copy()
		new.denominator = self.denominator.copy()
		new.candidates = self.candidates.copy()
		return new

	@staticmethod
	def coinc_params(events, offsetvector, triangulators):
		#
		# check for coincs that have been vetoed entirely
		#

		if len(events) < 2:
			return None

		#
		# Initialize the parameter dictionary, sort the events by
		# instrument name (the multi-instrument parameters are defined for
		# the instruments in this order and the triangulators are
		# constructed this way too), and retrieve the sorted instrument
		# names
		#

		params = {}
		events = tuple(sorted(events, key = lambda event: event.ifo))
		instruments = tuple(event.ifo for event in events)

		#
		# zero-instrument parameters
		#

		ignored, ignored, ignored, rss_timing_residual = triangulators[instruments](tuple(event.peak + offsetvector[event.ifo] for event in events))
		# FIXME:  rss_timing_residual is forced to 0 to disable this
		# feature.  all the code to compute it properly is still here and
		# given suitable initializations, the distribution data is still
		# two-dimensional and has a suitable filter applied to it, but all
		# events are forced into the RSS_{\Delta t} = 0 bin, in effect
		# removing that dimension from the data.  We can look at this again
		# sometime in the future if we're curious why it didn't help.  Just
		# delete the next line and you're back in business.
		rss_timing_residual = 0.0
		params["instrumentgroup,rss_timing_residual"] = (frozenset(instruments), rss_timing_residual)

		#
		# one-instrument parameters
		#

		for event in events:
			params["%s_snr2_chi2" % evemt.ifo] = (event.snr**2.0, event.chisq / event.chisq_dof)

		#
		# two-instrument parameters.  note that events are sorted by
		# instrument
		#

		for event1, event2 in itertools.combinations(events, 2):
			assert event1.ifo != event2.ifo

			prefix = "%s_%s_" % (event1.ifo, event2.ifo)

			dt = float((event1.peak + offsetvector[event1.ifo]) - (event2.peak + offsetvector[event2.ifo]))
			params["%sdt" % prefix] = (dt,)

			dA = math.log10(abs(event1.amplitude / event2.amplitude))
			params["%sdA" % prefix] = (dA,)

			# f_cut = central_freq + bandwidth/2
			f_cut1 = event1.central_freq + event1.bandwidth / 2
			f_cut2 = event2.central_freq + event2.bandwidth / 2
			df = float((math.log10(f_cut1) - math.log10(f_cut2)) / (math.log10(f_cut1) + math.log10(f_cut2)))
			params["%sdf" % prefix] = (df,)

		#
		# done
		#

		return params

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
		return self

	def to_xml(self, name):
		xml = ligolw.LIGO_LW({u"Name": u"%s:%s" % (name, self.ligo_lw_name_suffix)})
		xml.appendChild(self.numerator.to_xml("numerator"))
		xml.appendChild(self.denominator.to_xml("denominator"))
		xml.appendChild(self.candidates.to_xml("candidates"))
		return xml

	def add_slidelessbackground(self, database, experiments, param_func_args = ()):
		# FIXME:  this needs to be taught how to not slide H1 and
		# H2 with respect to each other

		# segment lists
		seglists = database.seglists - database.vetoseglists

		# construct the event list dictionary.  remove vetoed
		# events from the lists and save event peak times so they
		# can be restored later
		eventlists = {}
		orig_peak_times = {}
		for event in database.sngl_burst_table:
			if event.peak in seglists[event.ifo]:
				try:
					eventlists[event.ifo].append(event)
				except KeyError:
					eventlists[event.ifo] = [event]
				orig_peak_times[event] = event.peak

		# parse the --thresholds H1,L1=... command-line options from burca
		delta_t = [float(threshold.split("=")[-1]) for threshold in ligolw_process.get_process_params(database.xmldoc, "lalapps_burca", "--thresholds")]
		if not all(delta_t[0] == threshold for threshold in delta_t[1:]):
			raise ValueError("\Delta t is not unique in lalapps_burca arguments")
		delta_t = delta_t.pop()

		# construct the coinc generator.  note that H1+H2-only
		# coincs are forbidden, which is affected here by removing
		# that instrument combination from the object's internal
		# .rates dictionary
		coinc_generator = snglcoinc.CoincSynthesizer(eventlists, seglists, delta_t)
		if frozenset(("H1", "H2")) in coinc_generator.rates:
			del coinc_generator.rates[frozenset(("H1", "H2"))]

		# build a dictionary of time-of-arrival generators
		toa_generator = dict((instruments, coinc_generator.plausible_toas(instruments)) for instruments in coinc_generator.rates.keys())

		# how many coincs?  the expected number is obtained by
		# multiplying the total zero-lag time for which at least
		# two instruments were on by the sum of the rates for all
		# coincs to get the mean number of coincs per zero-lag
		# observation time, and multiplying that by the number of
		# experiments the background should simulate to get the
		# mean number of background events to simulate.  the actual
		# number simulated is a Poisson-distributed RV with that
		# mean.
		n_coincs, = scipy.stats.poisson.rvs(float(abs(segmentsUtils.vote(seglists.values(), 2))) * sum(coinc_generator.rates.values()) * experiments)

		# generate synthetic background coincs
		zero_lag_offset_vector = offsetvector.fromkeys(seglists, 0.0)
		for n, events in enumerate(coinc_generator.coincs(lsctables.SnglBurst.get_peak)):
			# n = 1 on 2nd iteration, so placing this condition
			# where it is in the loop causes the correct number
			# of events to be added to the background
			if n >= n_coincs:
				break
			# assign fake peak times
			toas = toa_generator[frozenset(event.ifo for event in events)].next()
			for event in events:
				event.peak = toas[event.ifo]
			# compute coincidence parameters
			self.denominator.increment(self.coinc_params(events, zero_lag_offset_vector, *param_func_args))

		# restore original peak times
		for event, peak_time in orig_peak_times.iteritems():
			event.peak = peak_time


#
# I/O
#


def load_likelihood_data(filenames, verbose = False):
	coinc_params = None
	seglists = None
	for n, filename in enumerate(filenames, 1):
		if verbose:
			print >>sys.stderr, "%d/%d:" % (n, len(filenames)),
		xmldoc = ligolw_utils.load_filename(filename, verbose = verbose, contenthandler = StringCoincParamsDistributions.contenthandler)
		this_coinc_params = StringCoincParamsDistributions.from_xml(xmldoc, u"string_cusp_likelihood")
		this_seglists = lsctables.SearchSummaryTable.get_table(xmldoc).get_out_segmentlistdict(set([this_coinc_params.process_id])).coalesce()
		xmldoc.unlink()
		if coinc_params is None:
			coinc_params = this_coinc_params
		else:
			coinc_params += this_coinc_params
		if seglists is None:
			seglists = this_seglists
		else:
			seglists |= this_seglists
	return coinc_params, seglists


#
# =============================================================================
#
#                                   Livetime
#
# =============================================================================
#


def time_slides_livetime(seglists, time_slides, min_instruments, verbose = False, clip = None):
	"""
	seglists is a segmentlistdict of times when each of a set of
	instruments were on, time_slides is a sequence of
	instrument-->offset dictionaries, each vector of offsets in the
	sequence is applied to the segmentlists and the total time during
	which at least min_instruments were on is summed and returned.  If
	clip is not None, after each offset vector is applied to seglists
	the result is intersected with clip before computing the livetime.
	If verbose is True then progress reports are printed to stderr.
	"""
	livetime = 0.0
	seglists = seglists.copy()	# don't modify original
	N = len(time_slides)
	if verbose:
		print >>sys.stderr, "computing the live time for %d time slides:" % N
	for n, time_slide in enumerate(time_slides):
		if verbose:
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
		seglists.offsets.update(time_slide)
		if clip is None:
			livetime += float(abs(segmentsUtils.vote(seglists.values(), min_instruments)))
		else:
			livetime += float(abs(segmentsUtils.vote((seglists & clip).values(), min_instruments)))
	if verbose:
		print >>sys.stderr, "\t100.0%"
	return livetime


def time_slides_livetime_for_instrument_combo(seglists, time_slides, instruments, verbose = False, clip = None):
	"""
	like time_slides_livetime() except computes the time for which
	exactly the instruments given by the sequence instruments were on
	(and nothing else).
	"""
	livetime = 0.0
	# segments for instruments that must be on
	onseglists = seglists.copy(keys = instruments)
	# segments for instruments that must be off
	offseglists = seglists.copy(keys = set(seglists) - set(instruments))
	N = len(time_slides)
	if verbose:
		print >>sys.stderr, "computing the live time for %s in %d time slides:" % (", ".join(instruments), N)
	for n, time_slide in enumerate(time_slides):
		if verbose:
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
		onseglists.offsets.update(time_slide)
		offseglists.offsets.update(time_slide)
		if clip is None:
			livetime += float(abs(onseglists.intersection(onseglists.keys()) - offseglists.union(offseglists.keys())))
		else:
			livetime += float(abs((onseglists & clip).intersection(onseglists.keys()) - offseglists.union(offseglists.keys())))
	if verbose:
		print >>sys.stderr, "\t100.0%"
	return livetime


#
# =============================================================================
#
#                              Database Utilities
#
# =============================================================================
#


def create_recovered_likelihood_table(connection, coinc_def_id):
	"""
	Create a temporary table named "recovered_likelihood" containing
	two columns:  "simulation_id", the simulation_id of an injection,
	and "likelihood", the highest likelihood ratio at which that
	injection was recovered by a coincidence of type coinc_def_id.
	"""
	cursor = connection.cursor()
	cursor.execute("""
CREATE TEMPORARY TABLE recovered_likelihood (simulation_id TEXT PRIMARY KEY, likelihood REAL)
	""")
	cursor.execute("""
INSERT OR REPLACE INTO
	recovered_likelihood
SELECT
	sim_burst.simulation_id AS simulation_id,
	MAX(coinc_event.likelihood) AS likelihood
FROM
	sim_burst
	JOIN coinc_event_map AS a ON (
		a.table_name == "sim_burst"
		AND a.event_id == sim_burst.simulation_id
	)
	JOIN coinc_event_map AS b ON (
		b.coinc_event_id == a.coinc_event_id
	)
	JOIN coinc_event ON (
		b.table_name == "coinc_event"
		AND b.event_id == coinc_event.coinc_event_id
	)
WHERE
	coinc_event.coinc_def_id == ?
GROUP BY
	sim_burst.simulation_id
	""", (coinc_def_id,))
	cursor.close()


def create_sim_burst_best_string_sngl_map(connection, coinc_def_id):
	"""
	Construct a sim_burst --> best matching coinc_event mapping.
	"""
	connection.cursor().execute("""
CREATE TEMPORARY TABLE
	sim_burst_best_string_sngl_map
AS
	SELECT
		sim_burst.simulation_id AS simulation_id,
		(
			SELECT
				sngl_burst.event_id
			FROM
				coinc_event_map AS a
				JOIN coinc_event_map AS b ON (
					b.coinc_event_id == a.coinc_event_id
				)
				JOIN coinc_event ON (
					coinc_event.coinc_event_id == a.coinc_event_id
				)
				JOIN sngl_burst ON (
					b.table_name == 'sngl_burst'
					AND b.event_id == sngl_burst.event_id
				)
			WHERE
				a.table_name == 'sim_burst'
				AND a.event_id == sim_burst.simulation_id
				AND coinc_event.coinc_def_id == ?
			ORDER BY
				(sngl_burst.chisq / sngl_burst.chisq_dof) / (sngl_burst.snr * sngl_burst.snr)
			LIMIT 1
		) AS event_id
	FROM
		sim_burst
	WHERE
		event_id IS NOT NULL
	""", (coinc_def_id,))


def create_sim_burst_best_string_coinc_map(connection, coinc_def_id):
	"""
	Construct a sim_burst --> best matching coinc_event mapping for
	string cusp injections and coincs.
	"""
	# FIXME:  this hasn't finished being ported from the inspiral code
	connection.cursor().execute("""
CREATE TEMPORARY TABLE
	sim_burst_best_string_coinc_map
AS
	SELECT
		sim_burst.simulation_id AS simulation_id,
		(
			SELECT
				coinc_inspiral.coinc_event_id
			FROM
				coinc_event_map AS a
				JOIN coinc_event_map AS b ON (
					b.coinc_event_id == a.coinc_event_id
				)
				JOIN coinc_inspiral ON (
					b.table_name == 'coinc_event'
					AND b.event_id == coinc_inspiral.coinc_event_id
				)
			WHERE
				a.table_name == 'sim_burst'
				AND a.event_id == sim_burst.simulation_id
				AND coinc_event.coinc_def_id == ?
			ORDER BY
				(sngl_burst.chisq / sngl_burst.chisq_dof) / (sngl_burst.snr * sngl_burst.snr)
			LIMIT 1
		) AS coinc_event_id
	FROM
		sim_burst
	WHERE
		coinc_event_id IS NOT NULL
	""", (coinc_def_id,))
