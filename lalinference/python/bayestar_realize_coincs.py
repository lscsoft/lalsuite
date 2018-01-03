#
# Copyright (C) 2013-2016  Leo Singer
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
"""
Synthesize triggers for simulated sources by realizing Gaussian measurement
errors in SNR and time of arrival. The input file (or stdin if the input file
is omitted) should be an optionally gzip-compressed LIGO-LW XML file of the form
produced by lalapps_inspinj. The output file (or stdout if omitted) will be an
optionally gzip-compressed LIGO-LW XML file containing single-detector triggers
and coincidences.

The root-mean square measurement error depends on the SNR of the signal, so
there is a choice for how to generate perturbed time and phase measurements:

 - zero-noise: no measurement error at all
 - from-truth: use true, nominal SNR in each detector
 - from-measurement: first perturb SNR with measurement error, then use
   that perturbed SNR to compute covariance of time and phase errors
"""
from __future__ import division
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Determine list of known detectors for command line arguments.
import lal
available_ifos = sorted(det.frDetector.prefix
    for det in lal.CachedDetectors)

# Command line interface.
import argparse
from lalinference.bayestar import command

parser = command.ArgumentParser()
parser.add_argument(
    'input', metavar='IN.xml[.gz]', type=argparse.FileType('rb'),
    default='-', help='Name of input file [default: stdin]')
parser.add_argument(
    '-o', '--output', metavar='OUT.xml[.gz]', type=argparse.FileType('wb'),
    default='-', help='Name of output file [default: stdout]')
parser.add_argument(
    '--detector', metavar='|'.join(available_ifos), nargs='+',
    help='Detectors to use [required].', choices=available_ifos, required=True)
parser.add_argument('--waveform',
    help='Waveform to use for injections (overrides values in '
    'sim_inspiral table)')
parser.add_argument('--snr-threshold', type=float, default=4.,
    help='Single-detector SNR threshold [default: %(default)s]')
parser.add_argument('--net-snr-threshold', type=float, default=12.,
    help='Network SNR threshold [default: %(default)s]')
parser.add_argument('--keep-subthreshold', action='store_true',
    help='Keep sub-threshold triggers that do not contribute to network SNR '
    '[default: %(default)s]')
parser.add_argument('--min-triggers', type=int, default=2,
    help='Emit coincidences only when at least this many triggers '
    'are found [default: %(default)s]')
parser.add_argument(
    '--measurement-error',
    choices=('zero-noise', 'from-truth', 'from-measurement'),
    default='zero-noise',
    help='How to compute the measurement error [default: %(default)s]')
parser.add_argument(
    '--reference-psd', metavar='PSD.xml[.gz]', type=argparse.FileType('rb'),
    required=True, help='Name of PSD file [required]')
parser.add_argument('--f-low', type=float,
    help='Override low frequency cutoff found in sim_inspiral table')
opts = parser.parse_args()


# Python standard library imports.
import os

# LIGO-LW XML imports.
from glue.ligolw import ligolw
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw.utils import search_summary as ligolw_search_summary
from glue.ligolw import table as ligolw_table
from pylal import ligolw_thinca
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

# glue, LAL and pylal imports.
from glue import segments
import glue.lal
import lal
import lal.series
import lalsimulation
from glue.text_progress_bar import ProgressBar

# BAYESTAR imports.
from lalinference.bayestar import ligolw as ligolw_bayestar
from lalinference.bayestar import filter
from lalinference.bayestar import timing

# Other imports.
import numpy as np


progress = ProgressBar()

# Open output file.
progress.update(-1, 'setting up output document')
out_xmldoc = ligolw.Document()
out_xmldoc.appendChild(ligolw.LIGO_LW())

# Write process metadata to output file.
process = command.register_to_xmldoc(
    out_xmldoc, parser, opts, ifos=opts.detector,
    comment="Simulated coincidences")

# Add search summary to output file.
all_time = segments.segment(
    [glue.lal.LIGOTimeGPS(0), glue.lal.LIGOTimeGPS(2e9)])
search_summary_table = lsctables.New(lsctables.SearchSummaryTable)
out_xmldoc.childNodes[0].appendChild(search_summary_table)
summary = ligolw_search_summary.append_search_summary(out_xmldoc, process,
    inseg=all_time, outseg=all_time)

# Read PSDs.
progress.update(-1, 'reading ' + opts.reference_psd.name)
xmldoc, _ = ligolw_utils.load_fileobj(
    opts.reference_psd, contenthandler=lal.series.PSDContentHandler)
psds = lal.series.read_psd_xmldoc(xmldoc)
psds = dict(
    (key, timing.InterpolatedPSD(filter.abscissa(psd), psd.data.data))
    for key, psd in psds.items() if psd is not None)

# Read injection file.
progress.update(-1, 'reading ' + opts.input.name)
xmldoc, _ = ligolw_utils.load_fileobj(
    opts.input, contenthandler=ligolw_bayestar.LSCTablesContentHandler)

# Extract simulation table from injection file.
sim_inspiral_table = ligolw_table.get_table(xmldoc,
    lsctables.SimInspiralTable.tableName)

# Create a SnglInspiral table and initialize its row ID counter.
sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
out_xmldoc.childNodes[0].appendChild(sngl_inspiral_table)
sngl_inspiral_table.set_next_id(lsctables.SnglInspiralID(0))

# Create a time slide entry.  Needed for coinc_event rows.
time_slide_table = lsctables.New(lsctables.TimeSlideTable)
out_xmldoc.childNodes[0].appendChild(time_slide_table)
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

# Precompute values that are common to all simulations.
detectors = [lalsimulation.DetectorPrefixToLALDetector(ifo)
    for ifo in opts.detector]
responses = [det.response for det in detectors]
locations = [det.location for det in detectors]

for sim_inspiral in progress.iterate(sim_inspiral_table):

    # Unpack some values from the row in the table.
    m1 = sim_inspiral.mass1
    m2 = sim_inspiral.mass2
    f_low = sim_inspiral.f_lower if opts.f_low is None else opts.f_low
    DL = sim_inspiral.distance
    ra = sim_inspiral.longitude
    dec = sim_inspiral.latitude
    inc = sim_inspiral.inclination
    phi = sim_inspiral.coa_phase
    psi = sim_inspiral.polarization
    epoch = lal.LIGOTimeGPS(
        sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
    gmst = lal.GreenwichMeanSiderealTime(epoch)
    waveform = sim_inspiral.waveform if opts.waveform is None else opts.waveform

    # FIXME: Set tranverse spin components to 0
    sim_inspiral.spin1x = 0
    sim_inspiral.spin1y = 0
    sim_inspiral.spin2x = 0
    sim_inspiral.spin2y = 0

    # Pre-evaluate some trigonometric functions that we will need.
    u = np.cos(inc)
    u2 = np.square(u)

    # Signal models for each detector.
    H = filter.sngl_inspiral_psd(sim_inspiral, waveform, f_low)
    signal_models = [
        timing.SignalModel(filter.signal_psd_series(H, psds[ifo]))
        for ifo in opts.detector]

    # Get SNR=1 horizon distances for each detector.
    horizons = np.asarray([signal_model.get_horizon_distance()
        for signal_model in signal_models])

    # Get antenna factors for each detector.
    Fplus, Fcross = np.asarray([
        lal.ComputeDetAMResponse(response, ra, dec, psi, gmst)
        for response in responses]).T

    # Compute TOAs at each detector.
    toas = np.asarray([lal.TimeDelayFromEarthCenter(location, ra, dec,
        epoch) for location in locations])

    # Compute SNR in each detector.
    snrs = (0.5 * (1 + u2) * Fplus + 1j * u * Fcross) * horizons / DL

    abs_snrs = np.abs(snrs)
    arg_snrs = np.angle(snrs)

    if opts.measurement_error == 'zero-noise':
        pass
    elif opts.measurement_error == 'from-truth':
        # If user asked, apply noise to amplitudes /before/ adding noise to TOAs and phases.

        # Add noise to SNR estimates.
        abs_snrs += np.random.randn(len(abs_snrs))

        for i, signal_model in enumerate(signal_models):
            arg_snrs[i], toas[i]  = np.random.multivariate_normal(
                [arg_snrs[i], toas[i]], signal_model.get_cov(abs_snrs[i]))
    elif opts.measurement_error == 'from-measurement':
        # Otherwise, by defualt, apply noise to TOAs and phases first.

        for i, signal_model in enumerate(signal_models):
            arg_snrs[i], toas[i]  = np.random.multivariate_normal(
                [arg_snrs[i], toas[i]], signal_model.get_cov(abs_snrs[i]))

        # Add noise to SNR estimates.
        abs_snrs += np.random.randn(len(abs_snrs))
    else:
        raise RuntimeError("This code should not be reached.")

    sngl_inspirals = []
    net_snr = 0.0
    count_triggers = 0

    # Loop over individual detectors and create SnglInspiral entries.
    for ifo, abs_snr, arg_snr, toa, horizon in zip(
            opts.detector, abs_snrs, arg_snrs, toas, horizons):

        # If SNR < threshold, then the injection is not found. Skip it.
        if abs_snr >= opts.snr_threshold:
            count_triggers += 1
            net_snr += np.square(abs_snr)
        elif not opts.keep_subthreshold:
            continue

        # Create SnglInspiral entry.
        sngl_inspiral = lsctables.SnglInspiral()
        for validcolumn in sngl_inspiral_table.validcolumns.keys():
            setattr(sngl_inspiral, validcolumn, None)
        sngl_inspiral.process_id = process.process_id
        sngl_inspiral.ifo = ifo
        sngl_inspiral.mass1 = m1
        sngl_inspiral.mass2 = m2
        sngl_inspiral.end_time = (epoch + toa).gpsSeconds
        sngl_inspiral.end_time_ns = (epoch + toa).gpsNanoSeconds
        sngl_inspiral.snr = abs_snr
        sngl_inspiral.coa_phase = np.angle(np.exp(1j * arg_snr))
        sngl_inspiral.eff_distance = horizon / sngl_inspiral.snr
        sngl_inspirals.append(sngl_inspiral)

    net_snr = np.sqrt(net_snr)

    # If too few triggers were found, then skip this event.
    if count_triggers < opts.min_triggers:
        continue

    # If network SNR < threshold, then the injection is not found. Skip it.
    if net_snr < opts.net_snr_threshold:
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
    coinc_inspiral.ifos = lsctables.ifos_from_instrument_set(
        sngl_inspiral.ifo for sngl_inspiral in sngl_inspirals)
    coinc_inspiral.end = lal.LIGOTimeGPS(
        sum(sngl_inspiral.end.ns() for sngl_inspiral in sngl_inspirals)
        // len(sngl_inspirals) * 1e-9) # FIXME: should only be detected sngls
    coinc_inspiral.mass = sim_inspiral.mass1 + sim_inspiral.mass2
    coinc_inspiral.mchirp = sim_inspiral.mchirp
    coinc_inspiral.combined_far = 0.0 # Not provided
    coinc_inspiral.false_alarm_rate = 0.0 # Not provided
    coinc_inspiral.minimum_duration = None # Not provided
    coinc_inspiral.snr = net_snr
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


# Record process end time.
progress.update(-1, 'writing ' + opts.output.name)
ligolw_process.set_process_end_time(process)

# Write output file.
with ligolw_utils.SignalsTrap():
  ligolw_utils.write_fileobj(out_xmldoc, opts.output,
      gz=(os.path.splitext(opts.output.name)[-1]==".gz"))
