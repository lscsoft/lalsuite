#!@PYTHON@
#
# Copyright (C) 2013  Leo Singer
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
from __future__ import division
"""
Synthesize triggers for simulated sources using a miniature matched-filter
detection pipeline. The input file (or stdin if the input file is omitted)
should be an optionally gzip-compressed LIGO-LW XML file of the form
produced by lalapps_inspinj. The output file (or stdout if omitted) will be an
optionally gzip-compressed LIGO-LW XML file containing single-detector triggers
and coincidences. A single template that has the same intrinsic parameters as
the injection is used.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Determine list of known detectors for command line arguments.
import lal
available_ifos = sorted(det.frDetector.prefix
    for det in lal.lalCachedDetectors)

# List of interpolation methods
available_interp_methods = [
    "catmull-rom", "lanczos", "nearest-neighbor", "quadratic-fit"]


# Command line interface.
from optparse import Option, OptionParser
from lalinference.bayestar import command

parser = OptionParser(
    formatter = command.NewlinePreservingHelpFormatter(),
    description = __doc__,
    usage="%prog [options] --template-bank TMPLTBANK.xml[.gz] [INPUT.xml[.gz]] [-o OUTPUT.xml[.gz]]",
    option_list = [
        Option("-o", "--output", metavar="OUTPUT.xml[.gz]", default="/dev/stdout",
            help="Name of output file [default: %default]"),
        Option("--detector", metavar='|'.join(available_ifos), action="append",
            help="Detectors to use.  May be specified multiple times.",
            choices=available_ifos),
        Option("--trigger-window", type=float, default=0.1, metavar="SECONDS",
            help="Search for a trigger across this many seconds before and after the time of the injection [default: %default]"),
        Option("--interp-method", metavar='|'.join(available_interp_methods),
               default="lanczos", choices=available_interp_methods,
            help="Trigger interpolation method [default: %default]"),
        Option("--interp-window", metavar="SAMPLES", type=int, default=2,
            help="Trigger interpolation window [default: %default]"),
        Option("--waveform",
            help="Waveform to use for injections"),
        Option("--snr-threshold", type=float, default=4.,
            help="Single-detector SNR threshold [default: %default]"),
        Option("--min-triggers", type=int, default=2,
            help="Emit coincidences only when at least this many triggers are found [default: %default]"),
        Option("-R", "--repeat-first-injection", type=int, default=None,
            help="Instead of performing each injection once, just perform the first injection this many times."),
        Option("--template-bank", metavar="TMPLTBANK.xml[.gz]",
            help="Name of template bank file (required)"),
        Option("--reference-psd", metavar="PSD.xml[.gz]",
            help="Name of PSD file (required)")
    ]
)
opts, args = parser.parse_args()
infilename = command.get_input_filename(parser, args)
command.check_required_arguments(parser, opts, 'waveform', 'template_bank', 'reference_psd')


# Python standard library imports.
import os
import signal
import sys

# LIGO-LW XML imports.
from glue.ligolw import ligolw
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw.utils import search_summary as ligolw_search_summary
from glue.ligolw import table as ligolw_table
from pylal import ligolw_thinca
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

# glue and LAL imports.
from glue import segments
import glue.text_progress_bar
import glue.lal
import lal, lalsimulation
import lal.series

# BAYESTAR imports.
from lalinference.bayestar import timing
from lalinference.bayestar import filter
from lalinference.bayestar import ligolw as ligolw_bayestar

# Other imports.
import numpy as np


progress = glue.text_progress_bar.ProgressBar()

template_approximant, template_amplitude_order, template_phase_order = \
    timing.get_approximant_and_orders_from_string(opts.waveform)


# FIXME: sample rate could be a command line option; template duration and data
# duration should be determined from chirp time
sample_rate = 4096 # sample rate in Hz
template_duration = 128 # template duration in seconds
template_length = sample_rate * template_duration # template length in samples
data_duration = 512 # data duration in seconds
data_length = sample_rate * data_duration # data length in samples


# Open output file.
out_xmldoc = ligolw.Document()
out_xmldoc.appendChild(ligolw.LIGO_LW())

# Write process metadata to output file.
process = ligolw_process.register_to_xmldoc(out_xmldoc, parser.get_prog_name(),
    opts.__dict__, ifos=opts.detector, comment="Little hope!")

# Add search summary to output file.
all_time = segments.segment([glue.lal.LIGOTimeGPS(0), glue.lal.LIGOTimeGPS(2e9)])
search_summary_table = lsctables.New(lsctables.SearchSummaryTable)
out_xmldoc.childNodes[0].appendChild(search_summary_table)
summary = ligolw_search_summary.append_search_summary(out_xmldoc, process,
    inseg=all_time, outseg=all_time)

# Read template bank file.
progress.update(-1, 'reading ' + opts.template_bank)
xmldoc = ligolw_utils.load_filename(opts.template_bank)

# Determine the low frequency cutoff from the template bank file.
template_bank_f_low = ligolw_bayestar.get_temlate_bank_f_low(xmldoc)

template_bank = ligolw_table.get_table(xmldoc,
    lsctables.SnglInspiralTable.tableName)

# Read injection file.
progress.update(-1, 'reading ' + infilename)
xmldoc = ligolw_utils.load_filename(infilename)

# Extract simulation table from injection file.
sim_inspiral_table = ligolw_table.get_table(xmldoc,
    lsctables.SimInspiralTable.tableName)

# Create a SnglInspiral table and initialize its row ID counter.
sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
out_xmldoc.childNodes[0].appendChild(sngl_inspiral_table)
sngl_inspiral_table.set_next_id(lsctables.SnglInspiralID(0))

# Create a time slide entry.  Needed for coinc_event rows.
try:
    time_slide_table = ligolw_table.get_table(out_xmldoc, lsctables.TimeSlideTable.tableName)
except ValueError:
    time_slide_table = lsctables.New(lsctables.TimeSlideTable)
    out_xmldoc.childNodes[0].appendChild(time_slide_table)
    time_slide_table.sync_next_id()
time_slide_id = time_slide_table.get_time_slide_id(
    dict((ifo, 0) for ifo in opts.detector), create_new=process)

# Create a CoincDef table and record a CoincDef row for
# sngl_inspiral <-> sngl_inspiral coincidences.
coinc_def_table = lsctables.New(lsctables.CoincDefTable)
out_xmldoc.childNodes[0].appendChild(coinc_def_table)
coinc_def = ligolw_thinca.InspiralCoincDef
coinc_def_id = coinc_def_table.get_next_id()
coinc_def.coinc_def_id = coinc_def_id
coinc_def_table.append(coinc_def)

# Create a CoincMap table.
coinc_map_table = lsctables.New(lsctables.CoincMapTable)
out_xmldoc.childNodes[0].appendChild(coinc_map_table)

# Create a CoincEvent table.
coinc_table = lsctables.New(lsctables.CoincTable)
out_xmldoc.childNodes[0].appendChild(coinc_table)

# Create a CoincInspiral table.
coinc_inspiral_table = lsctables.New(lsctables.CoincInspiralTable)
out_xmldoc.childNodes[0].appendChild(coinc_inspiral_table)

# Read PSD file.
progress.update(-1, 'reading ' + opts.reference_psd)
xmldoc = ligolw_utils.load_filename(opts.reference_psd)
psds = lal.series.read_psd_xmldoc(xmldoc)
psds = dict(
    (key, timing.InterpolatedPSD(filter.abscissa(psd), psd.data.data))
    for key, psd in psds.iteritems() if psd is not None)

# Detector noise PSD model
class psdfunc(object):

    def __init__(self, func):
        self.func = func

    def __call__(self, f):
        _f = np.atleast_1d(f)
        ret = np.empty(_f.shape, dtype=float)
        cond = (5 <= _f) & (_f <= 1800)
        ret[cond] = self.func(_f[cond])
        ret[~cond] = 0.
        if np.isscalar(f):
            ret = np.asscalar(ret)
        return ret

psdfuncs = [psdfunc(psds[ifo]) for ifo in opts.detector]

# Compute horizon distances for each template, for each detector
progress.update(-1, 'computing signal models')
horizons_bank = [
    [
        timing.SignalModel(sngl.mass1, sngl.mass2, S, template_bank_f_low, template_approximant, template_amplitude_order, template_phase_order).get_horizon_distance()
        for S in psdfuncs
    ] for sngl in template_bank]

# Generate templates for each unique set of intrinsic parameters
# FIXME: Get template_duration, template_approximant, template_amplitude_order, and template_phase_order from sngl_inspiral table.
progress.update(-1, 'computing template bank')
template_bank = [
    [
        filter.generate_template(sngl.mass1, sngl.mass2, S, template_bank_f_low, sample_rate, template_duration, template_approximant, template_amplitude_order, template_phase_order)
        for S in psdfuncs
    ] for sngl in template_bank]

# Generate PSDs for data coloring
progress.update(-1, 'computing PSDs')
def generate_psd(S):
    psd = lal.CreateREAL8FrequencySeries(None, lal.LIGOTimeGPS(0), 0, 1 / data_duration, filter.unitInverseHertz, data_length // 2 + 1)
    psd.data.data = S(filter.abscissa(psd))
    return psd
psds = [generate_psd(S) for ifo, S in zip(opts.detector, psdfuncs)]


if opts.repeat_first_injection:
    n_injections = opts.repeat_first_injection
else:
    n_injections = len(sim_inspiral_table)


def detect_sngls(ifos, data, horizons, templates):
    for ifo, x, horizon, zW in zip(ifos, data, horizons, templates):
        # Apply matched filter
        rho = filter.fftfilt(zW, x.data.data)[len(zW)-1:]

        # Find maximum index
        i0 = long(round(-(template_duration + float(x.epoch - end_time)) * sample_rate))
        di = long(round(sample_rate * opts.trigger_window))
        imax = np.argmax(filter.abs2(rho[i0 - di:i0 + di])) + i0 - di

        # If SNR < threshold, then the injection is not found. Skip it.
        if abs(rho[imax]) < opts.snr_threshold:
            continue

        # Interpolate time series
        imax, rhomax = filter.interpolate_max(imax, rho, opts.interp_window, method=opts.interp_method)
        tmax = x.epoch + (imax / sample_rate + template_duration)

        # Add SnglInspiral entry.
        sngl_inspiral = lsctables.SnglInspiral()
        for validcolumn in sngl_inspiral_table.validcolumns.iterkeys():
            setattr(sngl_inspiral, validcolumn, None)
        sngl_inspiral.process_id = process.process_id
        sngl_inspiral.ifo = ifo
        sngl_inspiral.mass1 = mass1
        sngl_inspiral.mass2 = mass2
        sngl_inspiral.mtotal = mass1 + mass2
        sngl_inspiral.mchirp = (mass1 * mass2)**0.6 * sngl_inspiral.mtotal**-0.2
        sngl_inspiral.end_time = tmax.gpsSeconds
        sngl_inspiral.end_time_ns = tmax.gpsNanoSeconds
        sngl_inspiral.snr = abs(rhomax)
        sngl_inspiral.coa_phase = -np.angle(rhomax) # minus sign to match gstlal_inspiral phase convention
        sngl_inspiral.eff_distance = horizon / sngl_inspiral.snr
        yield sngl_inspiral


def detect_net_snr_and_sngls(ifos, data, horizons, templates):
    sngls = list(detect_sngls(ifos, data, horizons, templates))
    return sum(sngl.snr for sngl in sngls), sngls


def inject(hplus, hcross, ifo, psd):
    # Generate colored noise
    x = filter.colored_noise(epoch, data_duration, sample_rate, psd)

    # Project injection for this detector.
    detector = lalsimulation.DetectorPrefixToLALDetector(ifo)
    lalsimulation.SimInjectDetectorStrainREAL8TimeSeries(
        x, hplus, hcross,
        ra, dec, psi,
        detector, None)

    # Done!
    return x


class keyboard_interrupt_handler(object):
    def __init__(self):
        self.interrupted = False
        signal.signal(signal.SIGINT, self)
    def __call__(self, signal, frame):
        self.interrupted = True
handler = keyboard_interrupt_handler()


for i_sim_inspiral in progress.iterate(range(n_injections), format='injection %d of ' + str(n_injections)):

    if handler.interrupted:
        print 'warning: interrupted, cleaning up'
        break

    if opts.repeat_first_injection:
        sim_inspiral = sim_inspiral_table[0]
    else:
        sim_inspiral = sim_inspiral_table[i_sim_inspiral]

    # Unpack some values from the row in the table.
    if not(opts.repeat_first_injection) or i_sim_inspiral == 0:
        mass1 = sim_inspiral.mass1
        mass2 = sim_inspiral.mass2
        spin1x = sim_inspiral.spin1x
        spin1y = sim_inspiral.spin1y
        spin1z = sim_inspiral.spin1z
        spin2x = sim_inspiral.spin2x
        spin2y = sim_inspiral.spin2y
        spin2z = sim_inspiral.spin2z
        f_low = sim_inspiral.f_lower
        DL = sim_inspiral.distance
        ra = sim_inspiral.longitude
        dec = sim_inspiral.latitude
        inc = sim_inspiral.inclination
        phi = sim_inspiral.coa_phase
        psi = sim_inspiral.polarization
        end_time = lal.LIGOTimeGPS(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
        epoch = end_time - 256 # GPS start time of data
        gmst = lal.GreenwichMeanSiderealTime(end_time)
        approximant, amplitude_order, phase_order = timing.get_approximant_and_orders_from_string(sim_inspiral.waveform)
        if approximant != lalsimulation.TaylorT4:
            raise ValueError("unrecognized approximant")

        # Generate injection
        hplus, hcross = lalsimulation.SimInspiralChooseTDWaveform(
            phi, 1 / sample_rate,
            mass1 * lal.LAL_MSUN_SI, mass2 * lal.LAL_MSUN_SI,
            spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
            f_low, f_low,
            DL * 1e6 * lal.LAL_PC_SI,
            inc, 0, 0,
            None, None,
            amplitude_order,
            phase_order,
            approximant)
        hplus.epoch += end_time
        hcross.epoch += end_time

    # Realize detector noise and add injection
    data = [inject(hplus, hcross, ifo, psd) for ifo, psd in zip(opts.detector, psds)]

    net_snr, sngl_inspirals = max(detect_net_snr_and_sngls(opts.detector, data, horizons, templates) for templates, horizons in zip(template_bank, horizons_bank))

    # If too few triggers were found, then skip this event.
    if len(sngl_inspirals) < opts.min_triggers:
        continue

    # Add Coinc table entry.
    coinc = lsctables.Coinc()
    coinc.coinc_event_id = coinc_table.get_next_id()
    coinc.process_id = process.process_id
    coinc.coinc_def_id = coinc_def_id
    coinc.time_slide_id = time_slide_id
    coinc.set_instruments(opts.detector)
    coinc.nevents = len(opts.detector)
    coinc.likelihood = None
    coinc_table.append(coinc)

    # Add CoincInspiral table entry.
    coinc_inspiral = lsctables.CoincInspiral()
    coinc_inspiral.coinc_event_id = coinc.coinc_event_id
    coinc_inspiral.snr = np.sqrt(np.sum(np.square([sngl_inspiral.snr
        for sngl_inspiral in sngl_inspirals])))
    coinc_inspiral.set_ifos([sngl_inspiral.ifo
        for sngl_inspiral in sngl_inspirals])
    coinc_inspiral.set_end(glue.lal.LIGOTimeGPS(float(np.mean([float(sngl_inspiral.get_end())
        for sngl_inspiral in sngl_inspirals]))))
    coinc_inspiral.false_alarm_rate = None
    coinc_inspiral.combined_far = None
    coinc_inspiral.mchirp = sngl_inspirals[0].mchirp
    coinc_inspiral.mass = sngl_inspirals[0].mtotal
    coinc_inspiral.minimum_duration = None
    coinc_inspiral_table.append(coinc_inspiral)

    # Record all sngl_inspiral records and associate them with coincidences.
    for sngl_inspiral in sngl_inspirals:
        # Give this sngl_inspiral record an id and add it to the table.
        sngl_inspiral.event_id = sngl_inspiral_table.get_next_id()
        sngl_inspiral_table.append(sngl_inspiral)

        # Add CoincMap entry.
        coinc_map = lsctables.CoincMap()
        coinc_map.coinc_event_id = coinc.coinc_event_id
        coinc_map.table_name = sngl_inspiral_table.tableName
        coinc_map.event_id = sngl_inspiral.event_id
        coinc_map_table.append(coinc_map)


signal.signal(signal.SIGINT, signal.SIG_DFL)


# Record process end time.
progress.update(-1, 'writing ' + opts.output)
ligolw_process.set_process_end_time(process)

# Write output file.
ligolw_utils.write_filename(out_xmldoc, opts.output,
    gz=(os.path.splitext(opts.output)[-1]==".gz"))


if handler.interrupted:
    sys.exit(1)
