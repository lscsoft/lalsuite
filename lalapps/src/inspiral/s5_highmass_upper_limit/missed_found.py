import sys
import sqlite3
from optparse import OptionParser
from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue.ligolw import utils
from glue.ligolw import table
from glue import segmentsUtils

from pylal import db_thinca_rings
from pylal import llwapp
from pylal import rate
from pylal import SimInspiralUtils
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

import matplotlib
matplotlib.use('Agg')
import pylab

from pylal import git_version
__author__ = "Chad Hanna <channa@ligo.caltech.edu>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

lsctables.LIGOTimeGPS = LIGOTimeGPS

class Summary(object):
	def __init__(self, opts, flist):
		self.segments = segments.segmentlistdict()
		self.non_inj_fnames = []
		self.inj_fnames = []
		self.found = {}
		self.missed = {}
		self.opts = opts
		self.veto_segments = segments.segmentlistdict()
		self.zero_lag_segments = {}
		self.instruments = []
		self.livetime = {}
		self.multi_burst_table = None
		self.coinc_inspiral_table = None

		for f in flist:
			if opts.verbose: print >> sys.stderr, "Gathering stats from: %s...." % (f,)
			working_filename = dbtables.get_connection_filename(f, tmp_path=opts.tmp_space, verbose = opts.verbose)
			connection = sqlite3.connect(working_filename)
			dbtables.DBTable_set_connection(connection)
			xmldoc = dbtables.get_xml(connection)

			# look for a sim table
			try:
				sim_inspiral_table = dbtables.lsctables.SimInspiralTable.get_table(xmldoc)
				self.inj_fnames.append(f)
				sim = True
			except ValueError:
				self.non_inj_fnames.append(f)
				sim = False

			# FIGURE OUT IF IT IS A BURST OR INSPIRAL RUN
			try:
				self.multi_burst_table = dbtables.lsctables.MultiBurstTable.get_table(xmldoc)
			except ValueError:
				self.multi_burst_table = None
			try:
				self.coinc_inspiral_table = dbtables.lsctables.CoincInspiralTable.get_table(xmldoc)
			except ValueError:
				self.coinc_inspiral_table = None
			if self.multi_burst_table and self.coinc_inspiral_table:
				print >>sys.stderr, "both burst and inspiral tables found.  Aborting"
				raise ValueError

			if not sim:
				self.get_instruments(connection)
				self.segments += self.get_segments(connection,xmldoc)
				#FIXME, don't assume veto segments are the same in every file!
				self.veto_segments = self.get_veto_segments(connection)

			dbtables.discard_connection_filename(f, working_filename, verbose = opts.verbose)
			dbtables.DBTable_set_connection(None)

		# remove redundant instruments
		self.instruments = list(set(self.instruments))
		# FIXME Do these have to be done by instruments?
		self.segments -= self.veto_segments

		# segments and livetime by instruments
		for i in self.instruments:
			self.zero_lag_segments[i] = self.segments.intersection(i) - self.segments.union(set(self.segments.keys()) - i)
			self.livetime[i] = float(abs(self.zero_lag_segments[i]))

	def get_segments(self, connection,xmldoc):
		if self.coinc_inspiral_table:
			segments = db_thinca_rings.get_thinca_zero_lag_segments(connection, program_name = self.opts.live_time_program)
		if self.multi_burst_table:
			#FIXME CWB case of rings not handled
			segments = llwapp.segmentlistdict_fromsearchsummary(xmldoc, self.opts.live_time_program).coalesce()
		return segments

	def get_veto_segments(self, connection):
		if self.coinc_inspiral_table:
			if self.opts.veto_segments_name is not None: return db_thinca_rings.get_veto_segments(connection, opts.self.veto_segments_name)
		# FIXME BURST CASE VETOS NOT HANDLED
		else: return segments.segmentlistdict()

	def get_instruments(self, connection):
		for i in connection.cursor().execute('SELECT DISTINCT(instruments) FROM coinc_event'):
			if i[0]: self.instruments.append(frozenset(lsctables.instrument_set_from_ifos(i[0])))

	def get_injections(self, instruments, FAR=float("inf")):
		injfnames = self.inj_fnames
		zero_lag_segments = self.zero_lag_segments[instruments]
		verbose = self.opts.verbose
		found = []
		missed = []
		print >>sys.stderr, ""
		for cnt, f in enumerate(injfnames):
			print >>sys.stderr, "getting injections below FAR: " + str(FAR) + ":\t%.1f%%\r" % (100.0 * cnt / len(injfnames),),
			working_filename = dbtables.get_connection_filename(f, tmp_path = opts.tmp_space, verbose = verbose)
			connection = sqlite3.connect(working_filename)
			dbtables.DBTable_set_connection(connection)
			xmldoc = dbtables.get_xml(connection)
			# DON'T BOTHER CONTINUING IF THE INSTRUMENTS OF INTEREST ARE NOT HERE
			instruments_in_this_file = []
			for i in connection.cursor().execute('SELECT DISTINCT(instruments) FROM coinc_event'):
				if i[0]: instruments_in_this_file.append(frozenset(lsctables.instrument_set_from_ifos(i[0])))
			if instruments not in instruments_in_this_file:
				connection.close()
				dbtables.discard_connection_filename(f, working_filename, verbose = verbose)
				dbtables.DBTable_set_connection(None)
				continue

			# WORK OUT CORRECT SEGMENTS FOR THIS FILE WHERE WE SHOULD SEE INJECTIONS
			segments = self.get_segments(connection, xmldoc)
			segments -= self.veto_segments
			#print thincasegments
			zero_lag_segments  = segments.intersection(instruments) - segments.union(set(segments.keys()) - instruments)
			###############

			# DEFINE THE INJECTION WAS MADE FUNCTION
			def injection_was_made(geocent_end_time, geocent_end_time_ns, zero_lag_segments = zero_lag_segments):
				"""
				return True if injection was made in the given segmentlist
				"""
				return lsctables.LIGOTimeGPS(geocent_end_time, geocent_end_time_ns) in zero_lag_segments

			connection.create_function("injection_was_made", 2, injection_was_made)
			make_sim_inspiral = lsctables.SimInspiralTable.get_table(dbtables.get_xml(connection)).row_from_cols

			# INSPIRAL
			if self.coinc_inspiral_table:
				for values in connection.cursor().execute("""
SELECT
  sim_inspiral.*,
  -- true if injection matched a coinc below the false alarm rate threshold
  EXISTS (
    SELECT
      *
    FROM
      coinc_event_map AS mapa
      JOIN coinc_event_map AS mapb ON (
        mapa.coinc_event_id == mapb.coinc_event_id
      )
      JOIN coinc_inspiral ON (
        mapb.table_name == "coinc_event"
        AND mapb.event_id == coinc_inspiral.coinc_event_id
      )
    WHERE
      mapa.table_name == "sim_inspiral"
      AND mapa.event_id == sim_inspiral.simulation_id
      AND coinc_inspiral.combined_far < ?
  )
FROM
  sim_inspiral
WHERE
  -- only interested in injections that were injected
  injection_was_made(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
				""", (FAR,)):
					sim = make_sim_inspiral(values)
					if values[-1]:
						found.append(sim)
					else:
						missed.append(sim)

			# BURSTS
			if self.multi_burst_table:
				for values in connection.cursor().execute("""
SELECT
  sim_inspiral.*,
  -- true if injection matched a coinc below the false alarm rate threshold
  EXISTS (
    SELECT
      *
    FROM
      coinc_event_map AS mapa
      JOIN coinc_event_map AS mapb ON (
        mapa.coinc_event_id == mapb.coinc_event_id
      )
      JOIN multi_burst ON (
        mapb.table_name == "coinc_event"
        AND mapb.event_id == multi_burst.coinc_event_id
      )
    WHERE
      mapa.table_name == "sim_inspiral"
      AND mapa.event_id == sim_inspiral.simulation_id
      AND multi_burst.false_alarm_rate < ?
  )
FROM
  sim_inspiral
WHERE
  -- only interested in injections that were injected
  injection_was_made(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
				""", (FAR,)):
					sim = make_sim_inspiral(values)
					if values[-1]:
						found.append(sim)
					else:
						missed.append(sim)
			# done
			dbtables.discard_connection_filename(f, working_filename, verbose = verbose)
			dbtables.DBTable_set_connection(None)

			print >>sys.stderr, "\nFound = %d Missed = %d" % (len(found), len(missed))
		return found, missed

	def set_instruments_to_calculate(self):
		if not self.opts.instruments: return self.instruments
		if self.opts.instruments in self.instruments:
			return frozenset(lsctables.instrument_set_from_ifos(i[0]))
		else:
			print >> sys.stderr, "Instruments %s do not exist in DB, nothing will be calculated" % (str(frozenset(lsctables.instrument_set_from_ifos(i[0]))),)
		return []


def parse_command_line():
	parser = OptionParser(version = "%prog CVS $Id$", usage = "%prog [options] [file ...]", description = "%prog computes mass/mass upperlimit")
	parser.add_option("--instruments", metavar = "name[,name,...]", help = "Set the list of instruments.  Required.  Example \"H1,H2,L1\"")
	parser.add_option("--live-time-program", default = "thinca", metavar = "name", help = "Set the name of the program whose rings will be extracted from the search_summary table.  Default = \"thinca\".")
	parser.add_option("--far", help = "FAR to use for injection finding instead of loudest event")
	parser.add_option("--veto-segments-name", default = "vetoes", help = "Set the name of the veto segments to use from the XML document.")
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("--verbose", action = "store_true", help = "Be verbose.")

	opts, filenames = parser.parse_args()

	if opts.instruments: opts.instruments = lsctables.instrument_set_from_ifos(opts.instruments)
	if not filenames:
		print >>sys.stderr, "must specify at least one database file"
		sys.exit(1)
	return opts, filenames

opts, filenames = parse_command_line()

summ = Summary(opts, filenames)

for instruments in summ.set_instruments_to_calculate():
	found, missed = summ.get_injections(instruments)
	pylab.semilogy([f.mchirp for f in found], [f.eff_dist_l for f in found],'.')
	pylab.semilogy([m.mchirp for m in missed], [m.eff_dist_l for m in missed],'k.')
	pylab.savefig('test.png')


