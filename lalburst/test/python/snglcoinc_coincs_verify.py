from bisect import bisect_left, bisect_right
import functools
import itertools
import math
import numpy
from optparse import OptionParser
import os
import random
import sys
from tqdm import tqdm

from lal import LIGOTimeGPS
from lalburst import burca
from lalburst import offsetvector
from lalburst import snglcoinc

import igwn_segments as segments
from igwn_ligolw import ligolw
from igwn_ligolw import lsctables
from igwn_ligolw import utils as ligolw_utils
from igwn_ligolw.utils import process as ligolw_process
from igwn_ligolw.utils import time_slide as ligolw_time_slide


#
# construct a synthetic burst trigger document
#


def make_input_test_document(
	segs = {
		"H1": segments.segment(LIGOTimeGPS(1e6), LIGOTimeGPS(1e6 + 300)),
		"L1": segments.segment(LIGOTimeGPS(1e6), LIGOTimeGPS(1e6 + 300)),
		"V1": segments.segment(LIGOTimeGPS(1e6), LIGOTimeGPS(1e6 + 300))
	},
	mean_rates = {"H1": 1.0, "L1": math.pi, "V1": math.e},
	n_offset_vectors = 3
):
	assert set(mean_rates) == set(segs), "need rates and segments for the same instruments"

	# start an XML document with some process metadata
	xmldoc = ligolw.Document()
	xmldoc.appendChild(ligolw.LIGO_LW())
	process = ligolw_process.register_to_xmldoc(xmldoc, "snglcoinc_coincs_verify", {}, instruments = segs.keys(), comment = "artificial data for coincidence engine test suite")

	# add time slide information
	for i in range(n_offset_vectors):
		ligolw_time_slide.get_time_slide_id(xmldoc, offsetvector.offsetvector((instrument, random.uniform(-28. * n_offset_vectors, +28. * n_offset_vectors)) for instrument in segs), create_new = process)

	# add a sngl_burst table
	snglbursttable = xmldoc.childNodes[-1].appendChild(lsctables.SnglBurstTable.new(["event_id", "ifo", "peak_time", "peak_time_ns", "process:process_id"]))
	snglbursttable.sync_next_id()

	# fill with random events
	for instrument, segment in segs.items():
		for i in range(numpy.random.poisson(mean_rates[instrument] * float(abs(segment)))):
			snglbursttable.append(snglbursttable.RowType(
				event_id = snglbursttable.get_next_id(),
				ifo = instrument,
				peak = segment[0] + random.uniform(0., float(abs(segment))),
				process_id = process.process_id
			))

	# done
	return xmldoc


#
# brute force all-possible-coincs algorithm to generate correct answer
#


def do_brute_force_coinc(xmldoc, delta_t = 0.015, min_instruments = 1):
	# register ourselves in the process metadata
	process = ligolw_process.register_to_xmldoc(xmldoc, "brute_force_coinc", {})

	# construct coinc tables interface for convenience
	coinctables = snglcoinc.CoincTables(xmldoc, burca.StringCuspBBCoincDef)

	# this is normally not the correct way to get the instrument list
	# from a document, but here we only care about instruments that
	# made triggers.  if others were included in the analysis we don't
	# care about them.
	instruments = list(set(row.ifo for row in lsctables.SnglBurstTable.get_table(xmldoc)))
	dt = dict((frozenset(pair), delta_t + snglcoinc.light_travel_time(*pair)) for pair in itertools.combinations(instruments, 2))

	# create look-up table of triggers indexed by instrument, sorted by
	# time
	sngls_by_instrument = dict((instrument, []) for instrument in instruments)
	for row in lsctables.SnglBurstTable.get_table(xmldoc):
		sngls_by_instrument[row.ifo].append(row)
	for instrument in instruments:
		sngls_by_instrument[instrument].sort(key = lambda row: row.peak)

	# iterate over time slides
	for time_slide_id, offsetvector in lsctables.TimeSlideTable.get_table(xmldoc).as_dict().items():
		print("offsets = %s" % str(offsetvector))

		# set of event_id sets to exclude from the coinc set
		# because they've already been found in a higher-order
		# coinc
		exclude = set()

		# iterate over coinc orders from highest to lowest
		# (quadruples, triples, doubles, ...)
		for n in range(len(instruments), min_instruments - 1, -1):
			# iterate over unique choices of that many
			# instruments
			for select_instruments in itertools.combinations(instruments, n):
				offsets = tuple(offsetvector[instrument] for instrument in select_instruments)
				total = functools.reduce(lambda a, b: a * b, map(lambda instrument: len(sngls_by_instrument[instrument]), select_instruments), 1)
				# iterate over every n-way combination of
				# triggers from the selected instruments
				for sngls in tqdm(itertools.product(*map(sngls_by_instrument.__getitem__, select_instruments)), total = total, desc = ", ".join(sorted(select_instruments))):
					# if any Delta t fails coincidence,
					# discard this combination of
					# triggers
					if any(abs(t_a - t_b) > dt[frozenset((instrument_a, instrument_b))] for (instrument_a, t_a), (instrument_b, t_b) in itertools.combinations(((instrument, sngl.peak + offset) for instrument, sngl, offset in zip(select_instruments, sngls, offsets)), 2)):
						continue
					# if this combination of event IDs
					# participated in a higher-order
					# coinc, discard this combination
					# of triggers
					event_ids = frozenset(sngl.event_id for sngl in sngls)
					if event_ids in exclude:
						continue
					# record this coinc
					coinctables.append_coinc(*coinctables.coinc_rows(process.process_id, time_slide_id, sngls, "sngl_burst"))
					exclude |= {frozenset(x) for i in range(n - 1, min_instruments - 1, -1) for x in itertools.combinations(event_ids, i)}

	# done
	process.set_end_time_now()
	return xmldoc


#
# use snglcoinc to construct coincs
#


class coincgen_doubles(snglcoinc.coincgen_doubles):
	class singlesqueue(snglcoinc.coincgen_doubles.singlesqueue):
		@staticmethod
		def event_time(event):
			return event.peak

	class get_coincs(object):
		def __init__(self, events):
			self.events = events
			self.times = tuple(map(coincgen_doubles.singlesqueue.event_time, events))

		def __call__(self, event_a, offset_a, coinc_window):
			peak = event_a.peak + offset_a
			return self.events[bisect_left(self.times, peak - coinc_window) : bisect_right(self.times, peak + coinc_window)]


def do_snglcoinc(xmldoc, delta_t = 0.015, min_instruments = 1):
	# register ourselves in the process metadata
	process = ligolw_process.register_to_xmldoc(xmldoc, "snglcoinc_coinc", {})

	# construct coinc tables interface for convenience
	coinctables = snglcoinc.CoincTables(xmldoc, burca.StringCuspBBCoincDef)

	# construct offset vector assembly graph
	time_slide_graph = snglcoinc.TimeSlideGraph(coincgen_doubles, coinctables.time_slide_index, delta_t, min_instruments = min_instruments, verbose = True)

	# collect coincs using incremental approach to minic online
	# searches
	with tqdm(desc = "snglcoinc", total = len(lsctables.SnglBurstTable.get_table(xmldoc))) as progress:
		for instrument, events in itertools.groupby(sorted(lsctables.SnglBurstTable.get_table(xmldoc), key = lambda row: (row.peak, row.ifo)), lambda event: event.ifo):
			events = tuple(events)
			progress.update(len(events))
			if time_slide_graph.push(instrument, events, max(event.peak for event in events)):
				for node, events in time_slide_graph.pull():
					coinctables.append_coinc(*coinctables.coinc_rows(process.process_id, node.time_slide_id, events, "sngl_burst"))
		for node, events in time_slide_graph.pull(flush = True):
			coinctables.append_coinc(*coinctables.coinc_rows(process.process_id, node.time_slide_id, events, "sngl_burst"))

	# done
	process.set_end_time_now()
	return xmldoc


#
# compare the two coincidence algorithm results
#


class summary(object):
	def __init__(self, xmldoc, algorithm_name):
		process_id, = lsctables.ProcessTable.get_table(xmldoc).get_ids_by_program(algorithm_name)

		coinc_def_id = lsctables.CoincDefTable.get_table(xmldoc).get_coinc_def_id(search = burca.StringCuspBBCoincDef.search, search_coinc_type = burca.StringCuspBBCoincDef.search_coinc_type, create_new = False)

		self.sngls = dict((row.event_id, row) for row in lsctables.SnglBurstTable.get_table(xmldoc))
		if len(lsctables.SnglBurstTable.get_table(xmldoc)) - len(self.sngls):
			raise ValueError("document contains %d duplicate sngl_burst rows" % (len(lsctables.SnglBurstTable.get_table(xmldoc)) - len(self.sngls)))
		self.coincs = dict((row.coinc_event_id, row) for row in lsctables.CoincTable.get_table(xmldoc) if row.coinc_def_id == coinc_def_id and row.process_id == process_id)
		self.offsetvectors = lsctables.TimeSlideTable.get_table(xmldoc).as_dict()

		coinc_event_ids = frozenset(self.coincs)
		coinc_map = [row for row in lsctables.CoincMapTable.get_table(xmldoc) if row.table_name == "sngl_burst" and row.coinc_event_id in coinc_event_ids]
		coinc_map.sort(key = lambda row: row.coinc_event_id)
		self.coinc_map = {}
		for time_slide_id in self.offsetvectors:
			self.coinc_map[time_slide_id] = dict((frozenset(row.event_id for row in rows), coinc_event_id) for coinc_event_id, rows in itertools.groupby(coinc_map, key = lambda row: row.coinc_event_id) if self.coincs[coinc_event_id].time_slide_id == time_slide_id)

	def get_sngls(self, event_ids):
		return tuple(self.sngls[event_id] for event_id in event_ids)

	def min_time(self, time_slide_id, event_ids):
		offsetvector = self.offsetvectors[time_slide_id]
		return min(sngl.peak + offsetvector[sngl.ifo] for sngl in self.get_sngls(event_ids))

	def print_summary(self, time_slide_id, event_ids):
		offsetvector = self.offsetvectors[time_slide_id]
		sngls = self.get_sngls(event_ids)
		peak_times = []
		for sngl in sorted(sngls, key = lambda row: row.ifo):
			peak_times.append(sngl.peak + offsetvector[sngl.ifo])
			print("\tevent ID %d: %s: %s + %g s = %s" % (sngl.event_id, sngl.ifo, sngl.peak, offsetvector[sngl.ifo], sngl.peak + offsetvector[sngl.ifo]))
		print("\tmax Delta t = %g s" % float(max(peak_times) - min(peak_times)))


#
# Main
#


def parse_command_line():
	parser = OptionParser(
	)
	parser.add_option("--mode", metavar = "mode", default = "test", help = "Select operating mode.  Allowed values are \"initialize\" (generate random test document), \"brute-force\" (apply brute-force algorithm to test document to construct correct-answer document), \"test\" (default;  apply snglcoinc algorithm to test document and compare to brute-force result contained in the document).")
	parser.add_option("--delta-t", metavar = "seconds", type = float, default = 0.015, help = "Coincidence window (default = 0.015 s).")
	parser.add_option("--min-instruments", metavar = "count", type = int, default = 1, help = "Require candidates to have at least this many instruments participating (default = 2).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")

	options, filenames = parser.parse_args()

	if options.mode not in ("initialize", "brute-force", "test"):
		raise ValueError("invalide --mode %s" % options.mode)
	if options.delta_t < 0.:
		raise ValueError("--delta-t must be >= 0")
	if options.min_instruments < 1:
		raise ValueError("--min-instruments must be >= 1")
	if not filenames:
		filenames = [os.path.join(os.environ.get("LAL_TEST_SRCDIR", "."), "snglcoinc_coincs_verify_input.xml.gz")]
	elif len(filenames) != 1:
		raise ValueError("only exactly one filename may be provided.  if no filenames are given, a default will be used.")

	return options, filenames[0]


options, filename = parse_command_line()


if options.mode == "initialize":
	#
	# generate a document with random triggers.  save to disk
	#

	xmldoc = make_input_test_document()
	ligolw_utils.write_filename(xmldoc, filename, verbose = options.verbose)

elif options.mode == "brute-force":
	#
	# load document containing triggers.  use all-possible-coincs
	# brute-force algorithm to compute coincs, save to disk
	#

	xmldoc = ligolw_utils.load_filename(filename, verbose = options.verbose)
	do_brute_force_coinc(xmldoc, delta_t = options.delta_t, min_instruments = options.min_instruments)
	ligolw_utils.write_filename(xmldoc, filename, verbose = options.verbose)

elif options.mode == "test":
	#
	# load document containing coincs generated by brute-force
	# algorithm.  recompute coincs using snglcoinc incremental
	# coincidence engine.  compare the two sets of coincs
	#

	xmldoc = ligolw_utils.load_filename(filename, verbose = options.verbose)
	do_snglcoinc(xmldoc, delta_t = options.delta_t, min_instruments = options.min_instruments)

	# intialize the summary objects.  this checks for duplicate singles
	algorithm_names = ("brute_force_coinc", "snglcoinc_coinc")
	summaries = [summary(xmldoc, algorithm_name) for algorithm_name in algorithm_names]

	# check that the time slide vectors are the same or we're wasting
	# our time
	if summaries[0].offsetvectors != summaries[1].offsetvectors:
		raise ValueError("documents do not contain identical offset vectors, or their IDs are not equivalent")

	# check that the candidate lists are identical
	m = 0
	print("\ncoincs in %s that are not in %s:" % algorithm_names)
	for time_slide_id in summaries[0].offsetvectors:
		for n, event_ids in enumerate(sorted(set(summaries[0].coinc_map[time_slide_id]) - set(summaries[1].coinc_map[time_slide_id]), key = lambda event_ids: summaries[0].min_time(time_slide_id, event_ids)), start = m):
			print("%d:" % n)
			summaries[0].print_summary(time_slide_id, event_ids)
			m += 1
	print("\ncoincs in %s that are not in %s:" % tuple(reversed(algorithm_names)))
	for time_slide_id in summaries[0].offsetvectors:
		for n, event_ids in enumerate(sorted(set(summaries[1].coinc_map[time_slide_id]) - set(summaries[0].coinc_map[time_slide_id]), key = lambda event_ids: summaries[1].min_time(time_slide_id, event_ids)), start = m):
			print("%d:" % n)
			summaries[1].print_summary(time_slide_id, event_ids)
			m += 1
	if m:
		raise ValueError("documents do not contain identical candidate sets")
