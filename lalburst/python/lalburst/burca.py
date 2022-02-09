# Copyright (C) 2006--2011,2013,2015--2019,2021  Kipp Cannon
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


from bisect import bisect_left, bisect_right
import itertools
import math
import sys


from ligo.lw import lsctables
from ligo.segments import PosInfinity
from . import snglcoinc


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from .git_version import date as __date__
from .git_version import version as __version__


#
# =============================================================================
#
#                                 Speed Hacks
#
# =============================================================================
#


class SnglBurst(lsctables.SnglBurst):
	__slots__ = ()

	#
	# compare self's peak time to the LIGOTimeGPS instance other.
	# allows bisection searches by GPS time to find ranges of triggers
	# quickly
	#

	def __lt__(self, other):
		return self.peak < other

	def __le__(self, other):
		return self.peak <= other

	def __eq__(self, other):
		return self.peak == other

	def __ne__(self, other):
		return self.peak != other

	def __ge__(self, other):
		return self.peak >= other

	def __gt__(self, other):
		return self.peak > other


#
# =============================================================================
#
#                          CoincTables Customizations
#
# =============================================================================
#


#
# For use with excess power.
#


ExcessPowerBBCoincDef = lsctables.CoincDef(search = "excesspower", search_coinc_type = 0, description = "sngl_burst<-->sngl_burst coincidences")


class ExcessPowerCoincTables(snglcoinc.CoincTables):
	def __init__(self, xmldoc):
		super(ExcessPowerCoincTables, self).__init__(xmldoc)

		# find the multi_burst table or create one if not found
		try:
			self.multibursttable = lsctables.MultiBurstTable.get_table(xmldoc)
		except ValueError:
			self.multibursttable = lsctables.New(lsctables.MultiBurstTable, ("process_id", "duration", "central_freq", "bandwidth", "snr", "confidence", "amplitude", "coinc_event_id"))
			xmldoc.childNodes[0].appendChild(self.multibursttable)

	def make_multi_burst(self, process_id, coinc_event_id, events, offset_vector):
		multiburst = self.multibursttable.RowType(
			process_id = process_id,
			coinc_event_id = coinc_event_id,
			# populate the false alarm rate with none.
			false_alarm_rate = None
		)

		# snr = root sum of ms_snr squares

		multiburst.snr = math.sqrt(sum(event.ms_snr**2.0 for event in events))

		# peak time = ms_snr squared weighted average of peak times
		# (no attempt to account for inter-site delays).
		# LIGOTimeGPS objects don't like being multiplied by
		# things, so the first event's peak time is used as a
		# reference epoch.

		t = events[0].peak + offset_vector[events[0].ifo]
		multiburst.peak = t + sum(event.ms_snr**2.0 * float(event.peak + offset_vector[event.ifo] - t) for event in events) / multiburst.snr**2.0

		# duration = ms_snr squared weighted average of durations

		multiburst.duration = sum(event.ms_snr**2.0 * event.duration for event in events) / multiburst.snr**2.0

		# central_freq = ms_snr squared weighted average of peak
		# frequencies

		multiburst.central_freq = sum(event.ms_snr**2.0 * event.peak_frequency for event in events) / multiburst.snr**2.0

		# bandwidth = ms_snr squared weighted average of bandwidths

		multiburst.bandwidth = sum(event.ms_snr**2.0 * event.bandwidth for event in events) / multiburst.snr**2.0

		# confidence = minimum of confidences

		multiburst.confidence = min(event.confidence for event in events)

		# "amplitude" = h_rss of event with highest confidence

		multiburst.amplitude = max(events, key = lambda event: event.ms_confidence).ms_hrss

		# done

		return multiburst

	def coinc_rows(self, process_id, time_slide_id, events, table_name):
		coinc, coincmaps = super(ExcessPowerCoincTables, self).coinc_rows(process_id, time_slide_id, events, table_name)
		coinc.insts = (event.ifo for event in events)
		return coinc, coincmaps, self.make_multi_burst(process_id, coinc.coinc_event_id, events, self.time_slide_index[time_slide_id])

	def append_coinc(self, coinc, coincmaps, multiburst):
		coinc = super(ExcessPowerCoincTables, self).append_coinc(coinc, coincmaps)
		multiburst.coinc_event_id = coinc.coinc_event_id
		self.multibursttable.append(multiburst)
		return coinc


#
# For use with string cusp
#


StringCuspBBCoincDef = lsctables.CoincDef(search = "StringCusp", search_coinc_type = 0, description = "sngl_burst<-->sngl_burst coincidences")


class StringCuspCoincTables(snglcoinc.CoincTables):
	def coinc_rows(self, process_id, time_slide_id, events, table_name):
		coinc, coincmaps = super(StringCuspCoincTables, self).coinc_rows(process_id, time_slide_id, events, table_name)
		coinc.insts = (event.ifo for event in events)
		return coinc, coincmaps


#
# =============================================================================
#
#                            Coincidence Generators
#
# =============================================================================
#


#
# For use with excess power coincidence test
#


class ep_coincgen_doubles(snglcoinc.coincgen_doubles):
	class singlesqueue(snglcoinc.coincgen_doubles.singlesqueue):
		@staticmethod
		def event_time(event):
			return event.peak


	class get_coincs(object):
		@staticmethod
		def max_edge_peak_delta(events):
			# the largest difference between any event's peak
			# time and either its start or stop times
			if not events:
				return 0.0
			return max(max(float(event.peak - event.start), float(event.start + event.duration - event.peak)) for event in events)

		@staticmethod
		def comparefunc(a, offseta, b, coinc_window):
			if abs(a.central_freq - b.central_freq) > (a.bandwidth + b.bandwidth) / 2:
				return True

			astart = a.start + offseta
			bstart = b.start
			if astart > bstart + b.duration + coinc_window:
				# a starts after the end of b
				return True

			if bstart > astart + a.duration + coinc_window:
				# b starts after the end of a
				return True

			# time-frequency times intersect
			return False

		def __init__(self, events):
			self.events = events
			# for this instance, replace the method with the
			# pre-computed value
			self.max_edge_peak_delta = self.max_edge_peak_delta(events)

		def __call__(self, event_a, offset_a, coinc_window):
			# event_a's peak time
			peak = event_a.peak

			# difference between event_a's peak and start times
			dt = float(peak - event_a.start)

			# largest difference between event_a's peak time
			# and either its start or stop times
			dt = max(dt, event_a.duration - dt)

			# add our own max_edge_peak_delta and the light
			# travel time between the two instruments (when
			# done, if event_a's peak time differs by more than
			# this much from the peak time of an event in this
			# list then it is *impossible* for them to be
			# coincident)
			dt += self.max_edge_peak_delta + coinc_window

			# apply time shift
			peak += offset_a

			# extract the subset of events from this list that
			# pass coincidence with event_a (use bisection
			# searches for the minimum and maximum allowed peak
			# times to quickly identify a subset of the full
			# list)
			return [event_b for event_b in self.events[bisect_left(self.events, peak - dt) : bisect_right(self.events, peak + dt)] if not self.comparefunc(event_a, offset_a, event_b, coinc_window)]



#
# For use with string coincidence test
#


class string_coincgen_doubles(ep_coincgen_doubles):
	class get_coincs(object):
		def __init__(self, events):
			self.events = events

		def __call__(self, event_a, offset_a, coinc_window):
			peak = event_a.peak + offset_a
			template_id = event_a.template_id
			return [event for event in self.events[bisect_left(self.events, peak - coinc_window) : bisect_right(self.events, peak + coinc_window)] if event.template_id == template_id]


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def burca(
	xmldoc,
	process_id,
	coincgen_doubles,
	CoincTables,
	coinc_definer_row,
	delta_t,
	ntuple_comparefunc = lambda events, offset_vector: False,
	min_instruments = 2,
	incremental = True,
	verbose = False
):
	#
	# prepare the coincidence table interface.
	#

	if verbose:
		print("indexing ...", file=sys.stderr)
	coinc_tables = CoincTables(xmldoc, coinc_definer_row)

	#
	# construct offset vector assembly graph
	#

	time_slide_graph = snglcoinc.TimeSlideGraph(coincgen_doubles, coinc_tables.time_slide_index, delta_t, min_instruments = min_instruments, verbose = verbose)

	#
	# collect events.
	#

	sngl_burst_table = lsctables.SnglBurstTable.get_table(xmldoc)
	if not incremental:
		# normal version:  push everything into the graph, then
		# pull out all coincs in one operation below using the
		# final flush
		for instrument, events in itertools.groupby(sorted(sngl_burst_table, key = lambda row: row.ifo), lambda event: event.ifo):
			time_slide_graph.push(instrument, tuple(events), PosInfinity)
	else:
		# slower diagnostic version.  simulate an online
		# incremental analysis by pushing events into the graph in
		# time order and collecting candidates as we go.  we still
		# do the final flush operation below.
		for instrument, events in itertools.groupby(sorted(sngl_burst_table, key = lambda row: (row.peak, row.ifo)), lambda event: event.ifo):
			events = tuple(events)
			if time_slide_graph.push(instrument, events, max(event.peak for event in events)):
				for node, events in time_slide_graph.pull(coinc_sieve = ntuple_comparefunc):
					coinc_tables.append_coinc(*coinc_tables.coinc_rows(process_id, node.time_slide_id, events, "sngl_burst"))

	#
	# retrieve all remaining coincidences.
	#

	for node, events in time_slide_graph.pull(coinc_sieve = ntuple_comparefunc, flush = True):
		coinc_tables.append_coinc(*coinc_tables.coinc_rows(process_id, node.time_slide_id, events, "sngl_burst"))

	#
	# done
	#

	return xmldoc
