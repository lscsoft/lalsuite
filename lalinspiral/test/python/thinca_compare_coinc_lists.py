import itertools
from optparse import OptionParser
import sys


from lalinspiral.thinca import InspiralCoincDef


from ligo.lw import ligolw
from ligo.lw import lsctables
from ligo.lw import utils as ligolw_utils


@lsctables.use_in
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	pass


def parse_command_line():
	parser = OptionParser(
	)
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")

	options, filenames = parser.parse_args()

	if len(filenames) != 2:
		raise ValueError("exactly two filenames are required")

	return options, filenames


options, filenames = parse_command_line()


class index(object):
	def __init__(self, xmldoc):
		coinc_def_id = lsctables.CoincDefTable.get_table(xmldoc).get_coinc_def_id(search = InspiralCoincDef.search, search_coinc_type = InspiralCoincDef.search_coinc_type, create_new = False)

		self.sngls = dict((row.event_id, row) for row in lsctables.SnglInspiralTable.get_table(xmldoc))
		print("document contains %d duplicate sngl_inspiral rows" % (len(lsctables.SnglInspiralTable.get_table(xmldoc)) - len(self.sngls)))
		self.coincs = dict((row.coinc_event_id, row) for row in lsctables.CoincTable.get_table(xmldoc) if row.coinc_def_id == coinc_def_id)
		self.offsetvectors = lsctables.TimeSlideTable.get_table(xmldoc).as_dict()

		coinc_event_ids = frozenset(self.coincs)
		coinc_map = [row for row in lsctables.CoincMapTable.get_table(xmldoc) if row.table_name == "sngl_inspiral" and row.coinc_event_id in coinc_event_ids]
		coinc_map.sort(key = lambda row: row.coinc_event_id)
		print("document is missing %d sngl_inspiral rows that appear in coincs" % (len(set(row.event_id for row in coinc_map) - set(self.sngls))))
		print("document has %d sngl_inspiral rows that do not appear in coincs" % (len(set(self.sngls) - set(row.event_id for row in coinc_map))))
		self.coinc_map = {}
		for time_slide_id in self.offsetvectors:
			self.coinc_map[time_slide_id] = dict((frozenset(row.event_id for row in rows), coinc_event_id) for coinc_event_id, rows in itertools.groupby(coinc_map, key = lambda row: row.coinc_event_id) if self.coincs[coinc_event_id].time_slide_id == time_slide_id)

	def get_sngls(self, event_ids):
		return tuple(self.sngls[event_id] for event_id in event_ids)

	def min_time(self, time_slide_id, event_ids):
		offsetvector = self.offsetvectors[time_slide_id]
		return min(sngl.end + offsetvector[sngl.ifo] for sngl in self.get_sngls(event_ids))

	def print_summary(self, time_slide_id, event_ids):
		offsetvector = self.offsetvectors[time_slide_id]
		sngls = self.get_sngls(event_ids)
		end_times = []
		for sngl in sorted(sngls, key = lambda row: row.ifo):
			end_times.append(sngl.end + offsetvector[sngl.ifo])
			print("\tevent ID %d: %s: %s + %g s = %s" % (sngl.event_id, sngl.ifo, sngl.end, offsetvector[sngl.ifo], sngl.end + offsetvector[sngl.ifo]))
		print("\tmax Delta t = %g s" % float(max(end_times) - min(end_times)))


indexes = [index(ligolw_utils.load_filename(filename, contenthandler = LIGOLWContentHandler, verbose = options.verbose)) for filename in filenames]

if indexes[0].offsetvectors != indexes[1].offsetvectors:
	raise ValueError("documents do not contain identical offset vectors, or their IDs are not equivalent")


print("\ncoincs in document 1 that are not in document 2:")
n = 0
for time_slide_id in indexes[0].offsetvectors:
	for n, event_ids in enumerate(sorted(set(indexes[0].coinc_map[time_slide_id]) - set(indexes[1].coinc_map[time_slide_id]), key = lambda event_ids: indexes[0].min_time(time_slide_id, event_ids)), start = n):
		print("%d:" % n)
		indexes[0].print_summary(time_slide_id, event_ids)

print("\ncoincs in document 2 that are not in document 1:")
n = 0
for time_slide_id in indexes[0].offsetvectors:
	for n, event_ids in enumerate(sorted(set(indexes[1].coinc_map[time_slide_id]) - set(indexes[0].coinc_map[time_slide_id]), key = lambda event_ids: indexes[1].min_time(time_slide_id, event_ids)), start = n):
		print("%d:" % n)
		indexes[1].print_summary(time_slide_id, event_ids)
