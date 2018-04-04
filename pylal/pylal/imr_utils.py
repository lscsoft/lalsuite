# Copyright (C) 2012  Chad Hanna
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

import sys
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue import segments
from glue import segmentsUtils
from glue.ligolw import table
from pylal import db_thinca_rings
from pylal import rate
import numpy
import math
import copy
from glue.ligolw.utils import search_summary as ligolw_search_summary
from glue.ligolw.utils import segments as ligolw_segments
from glue.ligolw.utils import process
from lalsimulation import SimInspiralTaylorF2ReducedSpinComputeChi, SimIMRPhenomBComputeChi

try:
	import sqlite3
except ImportError:
	# pre 2.5.x
	from pysqlite2 import dbapi2 as sqlite3


def allowed_analysis_table_names():
	return (dbtables.lsctables.MultiBurstTable.tableName, dbtables.lsctables.CoincInspiralTable.tableName, dbtables.lsctables.CoincRingdownTable.tableName)


def make_sim_inspiral_row_from_columns_in_db(connection):
	"""
	get the unique mapping of a sim inspiral row from columns in this
	database
	"""
	return lsctables.table.get_table(dbtables.get_xml(connection), lsctables.SimInspiralTable.tableName).row_from_cols


def time_within_segments(geocent_end_time, geocent_end_time_ns, zero_lag_segments = None):
	"""
	Return True if injection was made in the given segmentlist, if no
	segments just return True
	"""
	if zero_lag_segments is None:
		return True
	else:
		return lsctables.LIGOTimeGPS(geocent_end_time, geocent_end_time_ns) in zero_lag_segments


def get_min_far_inspiral_injections(connection, segments = None, table_name = "coinc_inspiral"):
	"""
	This function returns the found injections from a database and the
	minimum far associated with them as tuple of the form (far, sim). It also tells
	you all of the injections that should have been injected.  Subtracting the two
	outputs	should tell you the missed injections
	"""

	if table_name == dbtables.lsctables.CoincInspiralTable.tableName:
		found_query = 'SELECT sim_inspiral.*, coinc_inspiral.combined_far FROM sim_inspiral JOIN coinc_event_map AS mapA ON mapA.event_id == sim_inspiral.simulation_id JOIN coinc_event_map AS mapB ON mapB.coinc_event_id == mapA.coinc_event_id JOIN coinc_inspiral ON coinc_inspiral.coinc_event_id == mapB.event_id JOIN coinc_event on coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id WHERE mapA.table_name = "sim_inspiral" AND mapB.table_name = "coinc_event" AND injection_in_segments(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)'

	elif table_name == dbtables.lsctables.CoincRingdownTable.tableName:
		found_query = 'SELECT sim_inspiral.*, coinc_ringdown.false_alarm_rate FROM sim_inspiral JOIN coinc_event_map AS mapA ON mapA.event_id == sim_inspiral.simulation_id JOIN coinc_event_map AS mapB ON mapB.coinc_event_id == mapA.coinc_event_id JOIN coinc_ringdown ON coinc_ringdown.coinc_event_id == mapB.event_id JOIN coinc_event on coinc_event.coinc_event_id == coinc_ringdown.coinc_event_id WHERE mapA.table_name = "sim_inspiral" AND mapB.table_name = "coinc_event" AND injection_in_segments(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)'

	elif table_name == dbtables.lsctables.MultiBurstTable.tableName:
		found_query = 'SELECT sim_inspiral.*, multi_burst.false_alarm_rate FROM sim_inspiral JOIN coinc_event_map AS mapA ON mapA.event_id == sim_inspiral.simulation_id JOIN coinc_event_map AS mapB ON mapB.coinc_event_id == mapA.coinc_event_id JOIN multi_burst ON multi_burst.coinc_event_id == mapB.event_id JOIN coinc_event on coinc_event.coinc_event_id == multi_burst.coinc_event_id WHERE mapA.table_name = "sim_inspiral" AND mapB.table_name = "coinc_event" AND injection_in_segments(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)'

	else:
		raise ValueError("table must be in " + " ".join(allowed_analysis_table_names()))

	def injection_was_made(end_time, end_time_ns, segments = segments):
		return time_within_segments(end_time, end_time_ns, segments)

	# restrict the found injections to only be within certain segments
	connection.create_function("injection_in_segments", 2, injection_was_made)

	# get the mapping of a record returned by the database to a sim
	# inspiral row. Note that this is DB dependent potentially, so always
	# do this!
	make_sim_inspiral = make_sim_inspiral_row_from_columns_in_db(connection)

	found_injections = {}

	for values in connection.cursor().execute(found_query):
		# all but the last column is used to build a sim inspiral object
		sim = make_sim_inspiral(values[:-1])
		far = values[-1]
		# update with the minimum far seen until now
		this_inj = found_injections.setdefault(sim.simulation_id, (far, sim))
		if far < this_inj[0]:
			found_injections[sim.simulation_id] = (far, sim)

	total_query = 'SELECT * FROM sim_inspiral WHERE injection_in_segments(geocent_end_time, geocent_end_time_ns)'

	total_injections = {}
	# Missed injections start as a copy of the found injections
	missed_injections = {}
	for values in connection.cursor().execute(total_query):
		sim = make_sim_inspiral(values)
		total_injections[sim.simulation_id] = sim
		missed_injections[sim.simulation_id] = sim

	# now actually remove the missed injections
	for k in found_injections:
		del missed_injections[k]
		

	return found_injections.values(), total_injections.values(), missed_injections.values()


def get_max_snr_inspiral_injections(connection, segments = None, table_name = "coinc_inspiral"):
	"""
	Like get_min_far_inspiral_injections but uses SNR to rank injections.
	"""

	if table_name == dbtables.lsctables.CoincInspiralTable.tableName:
		found_query = 'SELECT sim_inspiral.*, coinc_inspiral.snr FROM sim_inspiral JOIN coinc_event_map AS mapA ON mapA.event_id == sim_inspiral.simulation_id JOIN coinc_event_map AS mapB ON mapB.coinc_event_id == mapA.coinc_event_id JOIN coinc_inspiral ON coinc_inspiral.coinc_event_id == mapB.event_id JOIN coinc_event on coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id WHERE mapA.table_name = "sim_inspiral" AND mapB.table_name = "coinc_event" AND injection_in_segments(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)'

	elif table_name in allowed_analysis_table_names():
		raise NotImplementedError("get_max_snr_inspiral_injections has not yet implemented querying against the table %s. Please consider submitting a patch. See get_min_far_inspiral_injections for how to construct your query." % table_name)
	else:
		raise ValueError("table must be in " + " ".join(allowed_analysis_table_names()))


	def injection_was_made(end_time, end_time_ns, segments = segments):
		return time_within_segments(end_time, end_time_ns, segments)

	# restrict the found injections to only be within certain segments
	connection.create_function("injection_in_segments", 2, injection_was_made)

	# get the mapping of a record returned by the database to a sim
	# inspiral row. Note that this is DB dependent potentially, so always
	# do this!
	make_sim_inspiral = make_sim_inspiral_row_from_columns_in_db(connection)

	found_injections = {}

	for values in connection.cursor().execute(found_query):
		# all but the last column is used to build a sim inspiral object
		sim = make_sim_inspiral(values[:-1])
		snr = values[-1]
		# update with the minimum far seen until now
		this_inj = found_injections.setdefault(sim.simulation_id, (snr, sim))
		if snr > this_inj[0]:
			found_injections[sim.simulation_id] = (snr, sim)

	total_query = 'SELECT * FROM sim_inspiral WHERE injection_in_segments(geocent_end_time, geocent_end_time_ns)'

	total_injections = {}
	# Missed injections start as a copy of the found injections
	missed_injections = {}
	for values in connection.cursor().execute(total_query):
		sim = make_sim_inspiral(values)
		total_injections[sim.simulation_id] = sim
		missed_injections[sim.simulation_id] = sim

	# now actually remove the missed injections
	for k in found_injections:
		del missed_injections[k]

	return found_injections.values(), total_injections.values(), missed_injections.values()


def get_instruments_from_coinc_event_table(connection):
	"""
	This function returns a list of the instruments analyzed according to the coinc_event_table
	"""
	instruments = []
	for ifos in connection.cursor().execute('SELECT DISTINCT(instruments) FROM coinc_event WHERE instruments!=""'):
		# ignore null columns
		if ifos[0]:
			instruments.append(frozenset(lsctables.instrument_set_from_ifos(ifos[0])))
	return instruments


def get_segments(connection, xmldoc, table_name, live_time_program, veto_segments_name = None, data_segments_name = "datasegments"):
	segs = segments.segmentlistdict()

	if table_name == dbtables.lsctables.CoincInspiralTable.tableName:
		if live_time_program == "gstlal_inspiral":
			segs = ligolw_segments.segmenttable_get_by_name(xmldoc, data_segments_name).coalesce()
			segs &= ligolw_search_summary.segmentlistdict_fromsearchsummary(xmldoc, live_time_program).coalesce()
		elif live_time_program == "thinca":
			segs = db_thinca_rings.get_thinca_zero_lag_segments(connection, program_name = live_time_program).coalesce()
		else:
			raise ValueError("for burst tables livetime program must be one of gstlal_inspiral, thinca")
		if veto_segments_name is not None:
			veto_segs = db_thinca_rings.get_veto_segments(connection, veto_segments_name)
			segs -= veto_segs
		return segs
	elif table_name == dbtables.lsctables.CoincRingdownTable.tableName:
		segs = ligolw_search_summary.segmentlistdict_fromsearchsummary(xmldoc, live_time_program).coalesce()
		if veto_segments_name is not None:
			veto_segs = ligolw_segments.segmenttable_get_by_name(xmldoc, veto_segments_name).coalesce()
			segs -= veto_segs
		return segs
	elif table_name == dbtables.lsctables.MultiBurstTable.tableName:
		if live_time_program == "omega_to_coinc":
			segs = ligolw_search_summary.segmentlistdict_fromsearchsummary(xmldoc, live_time_program).coalesce()
			if veto_segments_name is not None:
				veto_segs = ligolw_segments.segmenttable_get_by_name(xmldoc, veto_segments_name).coalesce()
				segs -= veto_segs
		elif live_time_program == "waveburst":
			segs = db_thinca_rings.get_thinca_zero_lag_segments(connection, program_name = live_time_program).coalesce()
			if veto_segments_name is not None:
				veto_segs = db_thinca_rings.get_veto_segments(connection, veto_segments_name)
				segs -= veto_segs
		else:
			raise ValueError("for burst tables livetime program must be one of omega_to_coinc, waveburst")
		return segs
	else:
		raise ValueError("table must be in " + " ".join(allowed_analysis_table_names()))


def get_event_fars(connection, table_name, segments = None):
	"""
	return the false alarm rate of the most rare zero-lag coinc by instruments
	"""

	def event_in_requested_segments(end_time, end_time_ns, segments = segments):
		return time_within_segments(end_time, end_time_ns, segments)

	connection.create_function("event_in_requested_segments", 2, event_in_requested_segments)

	if table_name == dbtables.lsctables.CoincInspiralTable.tableName:
		query = 'SELECT coinc_event.instruments, coinc_inspiral.combined_far AS combined_far, EXISTS(SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset != 0) FROM coinc_inspiral JOIN coinc_event ON (coinc_inspiral.coinc_event_id == coinc_event.coinc_event_id) WHERE event_in_requested_segments(coinc_inspiral.end_time, coinc_inspiral.end_time_ns);'

	elif table_name == dbtables.lsctables.MultiBurstTable.tableName:
		query = 'SELECT coinc_event.instruments, multi_burst.false_alarm_rate AS combined_far, EXISTS(SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset != 0) FROM multi_burst JOIN coinc_event ON (multi_burst.coinc_event_id == coinc_event.coinc_event_id) WHERE event_in_requested_segments(multi_burst.peak_time, multi_burst.peak_time_ns);'

	elif table_name == dbtables.lsctables.CoincRingdownTable.tableName:
		query = 'SELECT coinc_event.instruments, coinc_ringdown.false_alarm_rate AS combined_far, EXISTS(SELECT * FROM time_slide WHERE time_slide.time_slide_id == coinc_event.time_slide_id AND time_slide.offset != 0) FROM coinc_ringdown JOIN coinc_event ON (coinc_ringdown.coinc_event_id == coinc_event.coinc_event_id) WHERE event_in_requested_segments(coinc_ringdown.start_time, coinc_ringdown.start_time_ns);'

	else:
		raise ValueError("table must be in " + " ".join(allowed_analysis_table_names()))

	for inst, far, ts in connection.cursor().execute(query):
		inst = frozenset(lsctables.instrument_set_from_ifos(inst))
		yield (inst, far, ts)


def compute_search_efficiency_in_bins(found, total, ndbins, sim_to_bins_function = lambda sim: (sim.distance,)):
	"""
	This program creates the search efficiency in the provided ndbins.  The
	first dimension of ndbins must be the distance.  You also must provide a
	function that maps a sim inspiral row to the correct tuple to index the ndbins.
	"""

	input = rate.BinnedRatios(ndbins)

	# increment the numerator with the missed injections
	[input.incnumerator(sim_to_bins_function(sim)) for sim in found]

	# increment the denominator with the total injections
	[input.incdenominator(sim_to_bins_function(sim)) for sim in total]

	# regularize by setting empty bins to zero efficiency
	input.denominator.array[input.numerator.array < 1] = 1e35

	# pull out the efficiency array, it is the ratio
	eff = rate.BinnedArray(rate.NDBins(ndbins), array = input.ratio())

	# compute binomial uncertainties in each bin
	k = input.numerator.array
	N = input.denominator.array
	eff_lo_arr = ( N*(2*k + 1) - numpy.sqrt(4*N*k*(N - k) + N**2) ) / (2*N*(N + 1))
	eff_hi_arr = ( N*(2*k + 1) + numpy.sqrt(4*N*k*(N - k) + N**2) ) / (2*N*(N + 1))

	eff_lo = rate.BinnedArray(rate.NDBins(ndbins), array = eff_lo_arr)
	eff_hi = rate.BinnedArray(rate.NDBins(ndbins), array = eff_hi_arr)

	return eff_lo, eff, eff_hi


def compute_search_volume(eff):
	"""
	Integrate efficiency to get search volume.
	"""
	# get distance bins
	ndbins = eff.bins
	dx = ndbins[0].upper() - ndbins[0].lower()
	r = ndbins[0].centres()

	# we have one less dimension on the output
	vol = rate.BinnedArray(rate.NDBins(ndbins[1:]))

	# integrate efficiency to obtain volume
	vol.array = numpy.trapz(eff.array.T * 4. * numpy.pi * r**2, r, dx)

	return vol


def guess_nd_bins(sims, bin_dict = {"distance": (200, rate.LinearBins)}):
	"""
	Given a dictionary of bin counts and bin objects keyed by sim
	attribute, come up with a sensible NDBins scheme
	"""
	return rate.NDBins([bintup[1](min([getattr(sim, attr) for sim in sims]), max([getattr(sim, attr) for sim in sims]) + sys.float_info.min, bintup[0]) for attr, bintup in bin_dict.items()])


def guess_distance_mass1_mass2_bins_from_sims(sims, mass1bins = 11, mass2bins = 11, distbins = 200):
	"""
	Given a list of the injections, guess at the mass1, mass2 and distance
	bins.
	"""
	return guess_nd_bins(sims, bin_dict = {"distance": (distbins, rate.LinearBins), "mass1": (mass1bins, rate.LinearBins), "mass2": (mass2bins, rate.LinearBins)})


def guess_distance_spin1z_spin2z_bins_from_sims(sims, spin1bins = 11, spin2bins = 11, distbins = 200):
	"""
	Given a list of the injections, guess at the spin1, spin2 and distance
	bins.
	"""
	return guess_nd_bins(sims, bin_dict = {"distance": (distbins, rate.LinearBins), "spin1z": (spin1bins, rate.LinearBins), "spin2z": (spin2bins, rate.LinearBins)})


def guess_distance_effective_spin_parameter_bins_from_sims(sims, chibins = 11, distbins = 200):
	"""
	Given a list of the injections, guess at the chi = (m1*s1z +
	m2*s2z)/(m1+m2) and distance bins.
	"""
	dist_chi_vals = map(sim_to_distance_effective_spin_parameter_bins_function, sims)

	distances = [tup[0] for tup in dist_chi_vals]
	chis = [tup[1] for tup in dist_chi_vals]

	return rate.NDBins([rate.LinearBins(min(distances), max(distances), distbins), rate.LinearBins(min(chis), max(chis), chibins)])


def guess_distance_mass_ratio_bins_from_sims(sims, qbins = 11, distbins = 200):
	"""
	Given a list of the injections, guess at the chi and distance
	bins.
	"""
	dist_mratio_vals = map(sim_to_distance_mass_ratio_bins_function, sims)

	distances = [tup[0] for tup in dist_mratio_vals]
	mratios = [tup[1] for tup in dist_mratio_vals]

	return rate.NDBins([rate.LinearBins(min(distances), max(distances), distbins), rate.LinearBins(min(mratios), max(mratios), qbins)])


def guess_distance_chirp_mass_bins_from_sims(sims, mbins = 11, distbins = 200):
	"""
	Given a list of the injections, guess at the chirp mass and distance
	bins.
	"""
	dist_mchirp_vals = map(sim_to_distance_chirp_mass_bins_function, sims)

	distances = [tup[0] for tup in dist_mchirp_vals]
	mchirps = [tup[1] for tup in dist_mchirp_vals]

	return rate.NDBins([rate.LinearBins(min(distances), max(distances), distbins), rate.LinearBins(min(mchirps), max(mchirps), mbins)])


def guess_distance_total_mass_bins_from_sims(sims, nbins = 11, distbins = 200):
       """
       Given a list of the injections, guess at the mass1, mass2 and distance
       bins. Floor and ceil will be used to round down to the nearest integers.
       """

       total_lo = numpy.floor(min([sim.mass1 + sim.mass2 for sim in sims]))
       total_hi = numpy.ceil(max([sim.mass1 + sim.mass2 for sim in sims]))
       mindist = numpy.floor(min([sim.distance for sim in sims]))
       maxdist = numpy.ceil(max([sim.distance for sim in sims]))

       return rate.NDBins((rate.LinearBins(mindist, maxdist, distbins), rate.LinearBins(total_lo, total_hi, nbins)))


def sim_to_distance_mass1_mass2_bins_function(sim):
	"""
	create a function to map a sim to a distance, mass1, mass2 NDBins based object
	"""

	return (sim.distance, sim.mass1, sim.mass2)


def sim_to_distance_total_mass_bins_function(sim):
       """
       create a function to map a sim to a distance, total mass NDBins based object
       """
       return (sim.distance, sim.mass1 + sim.mass2)


def sim_to_distance_spin1z_spin2z_bins_function(sim):
	"""
	create a function to map a sim to a distance, spin1z, spin2z NDBins based object
	"""

	return (sim.distance, sim.spin1z, sim.spin2z)


def sim_to_distance_effective_spin_parameter_bins_function(sim):
	"""
	Map a sim_inspiral row to a distance, "chi" spin parameter
	bin. For IMR waveforms, "chi" refers to the effective spin,

	   chi = (m1*s1z + m2*s2z)/(m1 + m2)

	where s1z, s2z are the components of the spins along the
	direction of the total angular momentum. For inspiral
	waveforms, "chi" refers to the reduced spin,

	   chi_red = chi_s + delta*chi_a - 76.*eta/113*chi_s,

	where chi_s and chi_a are the symmetric and anti-symmetric
	combinations of the spins, and delta=(m1-m2)/(m1+m2). Some
	waveforms, e.g., SpinTaylorT4, use different coordinate
	conventions and require a coordinate transformation before
	applying these definitions.
	"""

	if sim.waveform.startswith("SpinTaylorT4"):
		chi1 = sim.spin1x * math.sin(sim.inclination) + sim.spin1z * math.cos(sim.inclination)
		chi2 = sim.spin2x * math.sin(sim.inclination) + sim.spin2z * math.cos(sim.inclination)
		chi = SimInspiralTaylorF2ReducedSpinComputeChi(sim.mass1, sim.mass2, chi1, chi2)

	elif sim.waveform.startswith("SpinTaylorT5"):
		chi1 = sim.spin1z
		chi2 = sim.spin2z
		chi = SimInspiralTaylorF2ReducedSpinComputeChi(sim.mass1, sim.mass2, chi1, chi2)

	elif sim.waveform.startswith("IMRPhenomB") or sim.waveform.startswith("IMRPhenomC") or sim.waveform.startswith("SEOBNR"):
		chi = SimIMRPhenomBComputeChi(sim.mass1, sim.mass2, sim.spin1z, sim.spin2z)

	else:
		raise ValueError(sim.waveform)

	return (sim.distance, chi)


def sim_to_distance_mass_ratio_bins_function(sim):
	"""
	create a function to map a sim to a distance, mass ratio NDBins based object
	"""
	# note that if you use symmetrize_sims() below, m2/m1 > 1
	# which just strikes me as more intuitive
	return (sim.distance, sim.mass2/sim.mass1)


def sim_to_distance_chirp_mass_bins_function(sim):
	"""
	create a function to map a sim to a distance, chirp mass NDBins based object
	"""
	return (sim.distance, sim.mchirp)

def symmetrize_sims(sims, col1, col2):
	"""
	symmetrize by two columns that should be symmetric.  For example mass1 and mass2
	"""
	for sim in sims:
		c1 = getattr(sim, col1)
		c2 = getattr(sim, col2)
		if c1 > c2:
			setattr(sim, col1, c2)
			setattr(sim, col2, c1)
	return sims

class DataBaseSummary(object):
	"""
	This class stores summary information gathered across the databases
	"""

	def __init__(self, filelist, live_time_program = None, veto_segments_name = None, data_segments_name = "datasegments", tmp_path = None, verbose = False):

		self.segments = segments.segmentlistdict()
		self.instruments = set()
		self.table_name = None
		self.found_injections_by_instrument_set = {}
		self.missed_injections_by_instrument_set = {}
		self.total_injections_by_instrument_set = {}
		self.zerolag_fars_by_instrument_set = {}
		self.ts_fars_by_instrument_set = {}
		self.numslides = set()

		for f in filelist:
			if verbose:
				print >> sys.stderr, "Gathering stats from: %s...." % (f,)
			working_filename = dbtables.get_connection_filename(f, tmp_path = tmp_path, verbose = verbose)
			connection = sqlite3.connect(working_filename)
			xmldoc = dbtables.get_xml(connection)

			sim = False

			# look for a sim inspiral table.  This is IMR work we have to have one of these :)
			try:
				sim_inspiral_table = table.get_table(xmldoc, dbtables.lsctables.SimInspiralTable.tableName)
				sim = True
			except ValueError:
				pass

			# look for the relevant table for analyses
			for table_name in allowed_analysis_table_names():
				try:
					setattr(self, table_name, table.get_table(xmldoc, table_name))
					if self.table_name is None or self.table_name == table_name:
						self.table_name = table_name
					else:
						raise ValueError("detected more than one table type out of " + " ".join(allowed_analysis_table_names()))
				except ValueError:
					setattr(self, table_name, None)

			# the non simulation databases are where we get information about segments
			if not sim:
				self.numslides.add(connection.cursor().execute('SELECT count(DISTINCT(time_slide_id)) FROM time_slide').fetchone()[0])
				[self.instruments.add(ifos) for ifos in get_instruments_from_coinc_event_table(connection)]
				# save a reference to the segments for this file, needed to figure out the missed and found injections
				self.this_segments = get_segments(connection, xmldoc, self.table_name, live_time_program, veto_segments_name, data_segments_name = data_segments_name)
				# FIXME we don't really have any reason to use playground segments, but I put this here as a reminder
				# self.this_playground_segments = segmentsUtils.S2playground(self.this_segments.extent_all())
				self.segments += self.this_segments

				# get the far thresholds for the loudest events in these databases
				for (instruments_set, far, ts) in get_event_fars(connection, self.table_name):
					if not ts:
						self.zerolag_fars_by_instrument_set.setdefault(instruments_set, []).append(far)
					else:
						self.ts_fars_by_instrument_set.setdefault(instruments_set, []).append(far)
			# get the injections
			else:
				# We need to know the segments in this file to determine which injections are found
				self.this_injection_segments = get_segments(connection, xmldoc, self.table_name, live_time_program, veto_segments_name, data_segments_name = data_segments_name)
				self.this_injection_instruments = []
				distinct_instruments = connection.cursor().execute('SELECT DISTINCT(instruments) FROM coinc_event WHERE instruments!=""').fetchall()
				for instruments, in distinct_instruments:
					instruments_set = frozenset(lsctables.instrument_set_from_ifos(instruments))
					self.this_injection_instruments.append(instruments_set)
					segments_to_consider_for_these_injections = self.this_injection_segments.intersection(instruments_set) - self.this_injection_segments.union(set(self.this_injection_segments.keys()) - instruments_set)
					found, total, missed = get_min_far_inspiral_injections(connection, segments = segments_to_consider_for_these_injections, table_name = self.table_name)
					if verbose:
						print >> sys.stderr, "%s total injections: %d; Found injections %d: Missed injections %d" % (instruments, len(total), len(found), len(missed))
					self.found_injections_by_instrument_set.setdefault(instruments_set, []).extend(found)
					self.total_injections_by_instrument_set.setdefault(instruments_set, []).extend(total)
					self.missed_injections_by_instrument_set.setdefault(instruments_set, []).extend(missed)

			# All done
			dbtables.discard_connection_filename(f, working_filename, verbose = verbose)
		if len(self.numslides) > 1:
			raise ValueError('number of slides differs between input files')
		elif self.numslides:
			self.numslides = min(self.numslides)
		else:
			self.numslides = 0

		# FIXME
		# Things left to do
		# 1) summarize the far threshold over the entire dataset
