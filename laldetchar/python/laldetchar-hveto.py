#
#  Copyright (C) 2013 Chris Pankow
#
#  This program is free software; ynu can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Fnundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
# Based on the reference implementation in MATLAB by Josh Smith, et al.
__doc__ = """trigger based hierarchical veto"""


import sys
import itertools
from collections import defaultdict
from optparse import OptionParser, OptionGroup

import numpy

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

from glue.lal import Cache
from glue.segments import segment, segmentlist, segmentlistdict
from glue.segmentsUtils import fromsegwizard
from glue.ligolw import ligolw, utils, lsctables, ilwd, table
from glue.ligolw.utils import process
from glue.ligolw.utils import segments as ligolw_segments

import laldetchar
import lal
# FIXME: Will this be needed with Karl's bindings?
# Yes... lalburst doesn't import the wrapped version of SnglBurst, so trying to
# get its attributes will fail until they are defined with this
# But that's been fixed, so we should try it without
import lalmetaio

# FIXME: gridspec...
#import breakout

from laldetchar import hveto
from laldetchar.hveto import glib_utils as gl
from laldetchar.hveto import plot_utils as pu

#
# Utilities
#

def process_options():
	"""
	Process options and check for required values.
	"""
	opt = OptionParser()
	opt.add_option( "-c", "--input-cache", help="Read triggers from the files in this cache." )
	opt.add_option( "-v", "--verbose", action="store_true", help="Be verbose." )
	veto_settings = OptionGroup( opt, "HVeto settings" )
	veto_settings.add_option( "-i", "--instrument", help="Instrument against which to veto. Required." )
	veto_settings.add_option( "-r", "--reference-channel", help="Channel against which to veto. Required." )
	veto_settings.add_option( "-t", "--reference-triggers", help="File path to load reference triggers. Required." )
	veto_settings.add_option( "-s", "--significance-threshold", type=float, default=15, help="Significance below which to terminate the rounds. Default is 15." )
	veto_settings.add_option( "--snr-thresh", action="append", help="Add an SNR threshold to use in veto round. Can be given multiple times for different values. WARNING: This will override the default settings, *not* append to them." )
	veto_settings.add_option( "--time-window", action="append", help="Add a time window to use in veto round. Can be given multiple times for different values. WARNING: This will override the default settings, *not* append to them." )
	veto_settings.add_option( "-S", "--min-ref-snr", type=float, default=8, help="Minimum SNR threshold to load a trigger in the reference channel." )
	# FIXME: Strictly speaking the ignore list is required because I'm not
	# sure what the	 function will do with out one?
	veto_settings.add_option( "-I", "--ignore-list", help="Text file, one channel per line with a list of channels to ignore when loading triggers." )
	# FIXME:
	#veto_settings.add_option( "-C", "--ignore-channel", action="append", help="Ignore these channels. Given several times, will ignore several channels. Do not prepend instrument. E.g. -C LSC-DARM_CTRL." )
	veto_settings.add_option( "--write-coinc", action="store_true", default=False, help="If set, output table will include coinc tables indicating which triggers were coincided in the process of execution." )
	opt.add_option_group( veto_settings )

	livetime_settings = OptionGroup( opt, "livetime settings" )
	# FIXME:
	#livetime_settings.add_option( "-L", "--livetime-definer", action="append", help="Name of segment definer entry from which to draw live segments. See '-l' option. If none is indicated, use all segments. Provide several times for deveral different segment definers." )
	livetime_settings.add_option( "-l", "--livetime-segments", help="File from which to parse livetime segments. Will assume, at first, a LIGOLW XML file with valid segment definer and segment tables. If this fails, will try segwizard format. Required." )
	livetime_settings.add_option( "--segment-definer", help="In tandem with --livetime-segments will retrieve segments with this definer. If none is provided, all segments will be used. Note: this option is REQUIRED if querying a databse. Example: H1:DMT-SCIENCE:1 (version is required)" )
	livetime_settings.add_option( "--segment-database", help="Query this URL for segments. Takes precedence over providing a file." )
	livetime_settings.add_option( "--gps-start", type=int, help="GPS start of analysis." )
	livetime_settings.add_option( "--gps-end", type=int, help="GPS end of analysis." )
	opt.add_option_group( livetime_settings )

	opts, args = opt.parse_args()
	if opts.instrument is None:
		print >>sys.stderr, "Instrument must be indicated."
		exit()
	if opts.reference_channel is None:
		print >>sys.stderr, "Reference channel must be indicated."
		exit()
	if opts.reference_triggers is None:
		print >>sys.stderr, "Reference triggers must be present."
		exit()
	if (opts.livetime_segments or opts.segment_database) is None:
		print >>sys.stderr, "Must provide livetime segments file or segment database location."
		exit()
	if opts.segment_database and (opts.segment_definer is None):
		print >>sys.stderr, "Must provide definer for segment database querying."
		exit()
	if len(args) == 0 and opts.input_cache is None:
		print >>sys.stderr, "Must provide input arguments or set --input-cache."
		exit()
	if opts.input_cache is not None:
		with open(opts.input_cache) as cache:
			c = Cache.fromfile(cache)
			args.extend( c.pfnlist() )
	if opts.ignore_list is None:
		print >>sys.stderr, "Must provide a channel ignore list."
		exit()

	return opts, args

#
# Default settings
#
DEFAULT_SNR_THRESH = [10.0, 12.0, 15.0, 20.0, 40.0, 100.0, 300.0]
DEFAULT_TIME_WINDOWS = [0.1, 0.2, 0.4, 0.8, 1.0]

# Begin the main script here
opt, args = process_options()

verbose = opt.verbose

#
# Thresholds
#

# SNR thresholds
snr_thresh = map(float, opt.snr_thresh or DEFAULT_SNR_THRESH)
# Minimum SNR of auxiliary channel to process
min_aux_chan_snr = snr_thresh[0]

# Window in which to check for coincidences
twind = map(float, opt.time_window or DEFAULT_TIME_WINDOWS)

# Threshold to stop the rounds
sig_stop_thresh = opt.significance_threshold

# Type of coincidence:
# 0 -- Record all coincidences
# 1 -- Record one or no coincidence per reference trigger
# NOTE: I don't provide this as an option yet. The original algorithm only used
# type 1, and type 0 may not even work after all the tooling I've done. Fair
# warning.
coinc_type = 1

# Channels to ignore
# FIXME: Seg faults will occur if this is a bad fileame
ignore_list = opt.ignore_list

#
# Reference channel loading
#

# Reference channel name
ref_ifo, ref_chan = opt.instrument, opt.reference_channel
# Reference channel filename
# FIXME: We shouldn't require the user to separate out the reference channel
ref_chan_fname = opt.reference_triggers

# Minimum SNR of reference channel to process
min_ref_chan_snr = opt.min_ref_snr

# Read in reference channel trigger files
trig_seq = laldetchar.PopulateTrigSequenceFromFile( None, ref_chan_fname, min_ref_chan_snr, ignore_list )

#
# Auxiliary channel loading
#

# Read in auxiliary channel trigger files
pstr = ""
for i, fname in enumerate(args):
	sys.stderr.write( "\b"*len(pstr) )
	sys.stderr.write( " "*len(pstr) )
	sys.stderr.write( "\b"*len(pstr) )
	pstr = "%d/%d %s (%d triggers)" % (i+1, len(args), fname, laldetchar.GetGSequenceLength( trig_seq ) )
	sys.stderr.write( pstr )
	trig_seq = laldetchar.PopulateTrigSequenceFromFile( trig_seq, fname, min_aux_chan_snr, ignore_list )
print ""

#
# Livetime segments
#

# Livetime and live segments
# Query from a segment database
if opt.segment_database:
	livesegs = hveto.query_segments_db( opt.segment_database, opt.gps_start, opt.gps_end, opt.segment_definer )
# Name of livetime file to process
elif opt.livetime_segments:
	# Get segments from an XML file
	try:
		livesegs = hveto.query_segments_xml( opt.livetime_segments, opt.gps_start, opt.gps_end, opt.segment_definer)
	# Okay, try segwizard
	except:
		livesegs = fromsegwizard( open(opt.livetime_segments) ).coalesce()
if not opt.gps_start and not opt.gps_end:
	runseg = segment(livesegs[0][0], livesegs[-1][1])
else:
	runseg = segment(opt.gps_start, opt.gps_end)
runlive = livetime = float(abs(livesegs))
print "Live time before vetoes: %f, covering %2.2f%% of %s" % (livetime, 100*livetime/abs(runseg), str(runseg))
if runlive == 0:
	exit()
lalsegl = gl.lalseg_from_seglist( livesegs )

chanlist = sorted(gl.get_chan_list( trig_seq ))
print "Channels present:\n\t" + "\n\t".join(chanlist)

# Preprocessing. Remove triggers which are outside the livetime segments or
# below our minimum SNR threshold
print "Length before prune: %d" % laldetchar.GetGSequenceLength( trig_seq )
laldetchar.DetCharPruneTrigs( trig_seq, lalsegl, min_ref_chan_snr, None );
total_ref_trigs = laldetchar.GetGSequenceLength( trig_seq )
print "Length after prune: %d" % total_ref_trigs

# Make a histogram of the reference channel triggers, and save it to compare
# against the final result
# FIXME: This is *slow*
"""
#preveto_snr = [sb.snr for sb in gl.seq_to_sbtable( trig_seq ) if sb.channel == ref_chan ]
preveto_snr = gl.seq_to_sbtable( trig_seq, get_col="snr" )
max_snr = numpy.log10( max(preveto_snr) )
snr_bins = numpy.logspace( numpy.log10(min_ref_chan_snr), max_snr, 30 )
snr_hist = pyplot.figure()
snr_hist_sub = snr_hist.add_subplot(111)
snr_hist_sub.hist( preveto_snr, bins=snr_bins, log=True, histtype='step', edgecolor='r', label='before round 1' )
"""

rnd = 1

ref_trigs = 0

sigdrop = None
round_winner, efficiency, deadtime = [], [], []
run_summary = []

#
# Round loop
#
while True:
	subround_winner = {}
	subround_count = {}
	subround_coinc = {}
	print "Round %d" % rnd
	for snr_t in snr_thresh:

		# Mark triggers not meeting threshold criteria
		laldetchar.DetCharPruneTrigs( trig_seq, lalsegl, snr_t, ref_chan )
		if laldetchar.CountUnmarkedSnglBurst( trig_seq ) <= 0:
			if verbose:
				print "No triggers remaining at current threshold. Moving on."
			continue

		for wind in twind:
			if verbose:
				print "SNR threshold %f, time window %f" % (snr_t, wind)
			chancount = laldetchar.CreateGHashTable(laldetchar.LALGTYPE_INT)
			chancoinc = laldetchar.CreateGHashTable(laldetchar.LALGTYPE_INT)
			# TODO: Wrap this better
			laldetchar.DetCharScanTrigs( chancount, chancoinc, trig_seq, ref_chan, wind, coinc_type )

			if not gl.ghash_key_in_dict( chancount, ref_chan ):
				if verbose:
					print "No reference channel triggers found."
				continue

			t_ratio = wind / livetime
			rnd_sig, winner = laldetchar.DetCharVetoRound( chancount, chancoinc, ref_chan, t_ratio )
			if verbose:
				print "Subround winner %s (SNR: %f, wind: %f): %f" % (winner, snr_t, wind, rnd_sig)

			chancount = gl.ghash_to_dict( chancount )
			subround_count[(snr_t, wind)] = chancount
			chancoinc = gl.ghash_to_dict( chancoinc )
			subround_coinc[(snr_t, wind)] = chancoinc

			if chancoinc.has_key(winner):
				mu = t_ratio * chancount[winner] * chancount[ref_chan]
				subround_winner[(snr_t, wind)] = (winner, rnd_sig, mu)
			#print chancount, chancoinc

	if verbose:
		print "Determining subround winner"
	subround_winner = sorted(subround_winner.iteritems(), key=lambda tup: tup[1][1] )

	if len(subround_winner) == 0:
		if verbose:
			print "No triggers remain, exiting."
		exit()
	round_winner = subround_winner[-1]
	snr_win, wind_win = round_winner[0]
	winner, sig, mu = round_winner[1]

	# Create veto segments
	start = lal.LIGOTimeGPS(-wind_win/2.0)
	stop = lal.LIGOTimeGPS(wind_win/2.0)
	veto = lal.SegCreate( start, stop, rnd )

	# We need a list of the reference channel triggers for plotting later
	# TODO: Reenable for breakout
	#all_ref_trigs = filter(lambda sb: sb.channel == ref_chan, gl.seq_to_sbtable( trig_seq ) )

	print "Retrieving vetoed triggers from list."
	vetoed_trigs = laldetchar.DetCharRemoveTrigs( trig_seq, veto, winner, None, snr_win )

	vetoed_trigs = gl.seq_to_sbtable( vetoed_trigs )
	# Number of reference channel triggers vetoed
	ref_chan_trig_vetoed = len( filter(lambda sb: sb.channel == ref_chan, vetoed_trigs ) )

	print "Calculating veto segments."
	# Remove vetoed winner triggers from livetime
	vetosegs = segmentlist([])
	for trig in vetoed_trigs:
		if trig.channel != winner: continue
		#trig = hveto.lalburst_sb_to_glue_sb( trig )
		vs = segment( trig.get_peak() - wind_win/2.0, trig.get_peak() + wind_win/2.0 )
		vetosegs.append( vs )
	vetosegs.coalesce()

	# Remove vetoed time from livetime
	deadt = float(abs(livesegs))
	livesegs -= vetosegs
	deadt -= float(abs(livesegs))

	# Write out vetoes and vetoed triggers
	xmldoc = hveto.write_round_xml( vetosegs, vetoed_trigs, winner, ref_ifo, opt.__dict__ )
	if opt.write_coinc:
		if opt.verbose:
			print "Appending coincs to output."
		hveto.write_coinc_tables( vetoed_trigs, xmldoc, ref_chan, wind_win )
	fname = "%s-HVETO_ROUND_%d-%d-%d.xml.gz" % (ref_ifo, rnd, int(runseg[0]), int(abs(runseg)))
	utils.write_filename( xmldoc, fname, gz=True, verbose=verbose )

	# Other information
	if rnd == 1:
		# Lowest SNR and window, for the total reference channel count
		snr_l, wind_l = snr_thresh[0], twind[0]
		ref_trigs = subround_count[(snr_l, wind_l)][ref_chan]
		efficiency.append( float(ref_chan_trig_vetoed)/ref_trigs*100 )
		inc_eff = efficiency[-1]
		deadtime.append( deadt/runlive*100 )
	else:
		inc_eff = float(ref_chan_trig_vetoed)/ref_trigs*100
		efficiency.append( inc_eff )
		deadtime.append( deadt/runlive*100 )

	print """
Round statistics:
\treference channel / winner: %s / %s
\tsnr threshold / time window: %f / %f
\tmu / significance: %f / %f
\tN_ref / vetoed = %d / %d
\tN_aux = %d
\tefficiency %% (rnd/cum): %f / %f
\tdeadtime %% (rnd/cum) %f / %f""" % ( ref_chan, winner, snr_win, wind_win, mu, sig, ref_trigs, ref_chan_trig_vetoed, subround_count[(snr_l, wind_l)][winner], inc_eff, sum(efficiency), deadtime[-1], sum(deadtime) )
	run_summary.append( {"winner": winner, "twin": wind_win, "snr": snr_win, "significance": sig, "nveto": ref_chan_trig_vetoed, "efficiency": inc_eff, "deadtime": deadtime[-1]} )

	# We've dropped below the useful threshold, bail.
	if sig < sig_stop_thresh:
		print "Last round below significance threshold, we're done here."
		break

	# Plots and the like
	chansig = defaultdict(dict)

	# Plot the vetoed reference channel tringgers
	if vetoed_trigs is not None or len(vetoed_trigs) != 0:
		aux_trigs = filter(lambda sb: sb.channel == winner, vetoed_trigs)
		pname = "hveto_round_%d_summary.png" % rnd
		#breakout.breakout_plot( all_ref_trigs, aux_trigs, ref_chan, winner, livesegs, vetosegs, pname )

	# Calculate the significance of each of the channels against the reference
	# channel for each of the subrounds for the significance drop plot
	for subrnd, coinctbl in subround_coinc.iteritems():
		ref_cnt = subround_count[subrnd][ref_chan]
		snr_t_sub, twind_sub = subrnd
		for chan, coinc in coinctbl.iteritems():
			cnt = subround_count[subrnd][chan]
			mu = twind_sub*cnt*ref_cnt/livetime
			chansig[subrnd][chan] = laldetchar.DetCharHvetoSignificance( mu, coinc )

	# If this isn't the first round, calculate the significance drop from
	# last round
	if sigdrop is not None:
		pu.plot_sigdrop( sigdrop, chansig[(snr_win, wind_win)], rnd )
	sigdrop = chansig[(snr_win, wind_win)]

	livetime = float(abs(livesegs))
	print "Live time after round %d: %f" % (rnd, livetime)

	#break
	#if rnd_sig < sig_stop_thresh:
	if sig < sig_stop_thresh:
		break
	else:
		rnd +=	1

# Final plots
pyplot.figure()
efficiency = numpy.cumsum( efficiency )
deadtime = numpy.cumsum( deadtime )
pyplot.xlabel( "Deadtime (%)" )
pyplot.ylabel( "Efficiency (%)" )
pyplot.plot( deadtime, efficiency, "k-" )
pyplot.grid()
pyplot.savefig( "effdead.png" )

print "Run Summary " + "="*68
for i, round in enumerate(run_summary):
	if i == 0:
		print "\t".join( round.keys() )
	print "\t".join( map(str, round.values()) )
print "="*80

# Finish up the SNR histogram after the veto rounds
"""
postveto_snr = [sb.snr for sb in gl.seq_to_sbtable( trig_seq ) if sb.channel == ref_chan ]

snr_hist_sub.hist( postveto_snr, bins=snr_bins, log=True, histtype='step', edgecolor='b', label="after round %d" % rnd )
snr_hist_sub.set_title( "Summary of vetoed triggers in channel %s\nTime offset from GPS %d, number of rounds %d" % (ref_chan, runseg[0], rnd) )
snr_hist_sub.set_xlabel("SNR")
snr_hist_sub.set_ylabel("Number")
snr_hist_sub.legend()
snr_hist.savefig( "overall_snr_hist.png" )
"""
