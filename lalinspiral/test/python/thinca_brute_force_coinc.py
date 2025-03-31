##python
import itertools
from optparse import OptionParser
import sys
from tqdm import tqdm


from lalburst import snglcoinc
from lalinspiral.thinca import InspiralCoincDef


from igwn_ligolw import lsctables
from igwn_ligolw import utils as ligolw_utils
from igwn_ligolw.utils import process as ligolw_process


def parse_command_line():
	parser = OptionParser(
	)
	parser.add_option("--delta-t", metavar = "seconds", type = float, default = 0.005, help = "Coincidence window (default = 0.005 s).")
	parser.add_option("--min-instruments", metavar = "count", type = int, default = 2, help = "Require candidates to have at least this many instruments participating (default = 2).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")

	options, filenames = parser.parse_args()

	process_params = options.__dict__.copy()

	if options.delta_t < 0.:
		raise ValueError("--delta-t must be >= 0")
	if options.min_instruments < 1:
		raise ValueError("--min-instruments must be >= 1")

	return options, process_params, filenames


options, process_params, filenames = parse_command_line()


for filename in filenames:
	xmldoc = ligolw_utils.load_filename(filename, verbose = options.verbose)

	process = ligolw_process.register_to_xmldoc(xmldoc, "brute_force_coinc", process_params)

	# don't bother with coinc_inspiral table
	coinctables = snglcoinc.CoincTables(xmldoc, InspiralCoincDef)

	# this is normally not the correct way to get the instrument list
	# from a document, but here we only care about instruments that
	# made triggers.  if others were included in the analysis we don't
	# care about them.
	instruments = list(set(row.ifo for row in lsctables.SnglInspiralTable.get_table(xmldoc)))
	dt = dict((frozenset(pair), options.delta_t + snglcoinc.light_travel_time(*pair)) for pair in itertools.combinations(instruments, 2))

	sngls = dict.fromkeys(row.template_id for row in lsctables.SnglInspiralTable.get_table(xmldoc))
	for template_id in sngls:
		sngls[template_id] = dict((instrument, []) for instrument in instruments)
	for row in lsctables.SnglInspiralTable.get_table(xmldoc):
		sngls[row.template_id][row.ifo].append(row)
	for template_id in sngls:
		for instrument in instruments:
			sngls[template_id][instrument].sort(key = lambda row: row.end)

	offsetvectors = lsctables.TimeSlideTable.get_table(xmldoc).as_dict()

	# iterate over instrument-->trigger list dictionaries, one for each
	# template

	for sngls_by_instrument in tqdm(sngls.values()):
		for time_slide_id, offsetvector in offsetvectors.items():
			# set of event_id sets to exclude from the coinc
			# set because they've already been found in a
			# higher-order coinc
			exclude = set()

			for n in range(len(instruments), options.min_instruments - 1, -1):
				for select_instruments in itertools.combinations(instruments, n):
					offsets = tuple(offsetvector[instrument] for instrument in select_instruments)
					for sngls in itertools.product(*map(sngls_by_instrument.__getitem__, select_instruments)):
						event_ids = frozenset(sngl.event_id for sngl in sngls)
						if event_ids in exclude:
							continue
						ends = tuple((instrument, sngl.end + offset) for instrument, sngl, offset in zip(select_instruments, sngls, offsets))
						if any(abs(t_a - t_b) > dt[frozenset((instrument_a, instrument_b))] for (instrument_a, t_a), (instrument_b, t_b) in itertools.combinations(ends, 2)):
							continue
						coinctables.append_coinc(*coinctables.coinc_rows(process.process_id, time_slide_id, sngls, "sngl_inspiral"))
						exclude |= {frozenset(x) for i in range(n - 1, options.min_instruments - 1, -1) for x in itertools.combinations(event_ids, i)}

	process.set_end_time_now()

	ligolw_utils.write_filename(xmldoc, filename, verbose = options.verbose)
