#
# Copyright (C) 2013-2017  Leo Singer
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
Produce GW sky maps for all coincidences in a LIGO-LW XML file.

The filename of the (optionally gzip-compressed) LIGO-LW XML input is an
optional argument; if omitted, input is read from stdin.

The distance prior is controlled by the --prior-distance-power argument.
If you set --prior-distance-power=k, then the distance prior is
proportional to r^k. The default is 2, uniform in volume.

If the --min-distance argument is omitted, it defaults to zero. If the
--max-distance argument is omitted, it defaults to the SNR=4 horizon
distance of the most sensitive detector.

A FITS file is created for each sky map, having a filename of the form

  "X.toa_phoa_snr.fits"
  "X.toa_snr_mcmc.fits"
  "X.toa_phoa_snr_mcmc.fits"

where X is the LIGO-LW row id of the coinc and "toa" or "toa_phoa_snr"
identifies whether the sky map accounts for times of arrival (TOA),
PHases on arrival (PHOA), and amplitudes on arrival (SNR).
"""


# Command line interface.
import argparse
from lalinference.bayestar import command

methods = '''
    toa_phoa_snr
    toa_phoa_snr_mcmc
    '''.split()
default_method = 'toa_phoa_snr'
command.skymap_parser.add_argument(
    '--method', choices=methods, default=[default_method], nargs='*',
    help='Sky localization methods [default: %(default)s]')
parser = command.ArgumentParser(
    parents=[
        command.waveform_parser, command.prior_parser, command.skymap_parser])
parser.add_argument(
    '--keep-going', '-k', default=False, action='store_true',
    help='Keep processing events if a sky map fails to converge [default: no]')
parser.add_argument(
    'input', metavar='INPUT.xml[.gz]', default='-', nargs='+',
    type=argparse.FileType('rb'),
    help='Input LIGO-LW XML file [default: stdin] or PyCBC HDF5 files. If '
    'PyCBC files, must be bank file, coinc file, and trigger files, in that '
    'order.')
parser.add_argument(
    '--pycbc-sample', default='foreground',
    help='sample population [PyCBC only; default: %(default)s]')
parser.add_argument(
    '--psd-files', nargs='*',
    help='pycbc-style merged HDF5 PSD files')
parser.add_argument(
    '--coinc-event-id', type=int, nargs='*',
    help='run on only these specified events')
parser.add_argument(
    '--output', '-o', default='.',
    help='output directory [default: current directory]')
parser.add_argument(
    '--condor-submit', action='store_true',
    help='submit to Condor instead of running locally')
opts = parser.parse_args()

#
# Logging
#

import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('BAYESTAR')

# LIGO-LW XML imports.
import lal.series
from glue.ligolw import utils as ligolw_utils

# BAYESTAR imports.
from lalinference.bayestar.decorator import memoized
from lalinference.io import fits
from lalinference.bayestar import ligolw as ligolw_bayestar
from lalinference.bayestar import filter
from lalinference.bayestar import timing
from lalinference.bayestar.sky_map import ligolw_sky_map

# Other imports.
import os
import sys
import numpy as np
import h5py

# Read coinc file.
log.info('%s:reading input files', ','.join(file.name for file in opts.input))
if os.path.splitext(opts.input[0].name)[1] in {'.hdf', '.hdf5', '.h5'}:
    hdf5files = [h5py.File(file.name, 'r') for file in opts.input]

    snr_dict = {}

    coinc_and_sngl_inspirals = ligolw_bayestar.coinc_and_sngl_inspirals_for_hdf(
        *hdf5files, sample=opts.pycbc_sample)
else:
    infile, = opts.input
    xmldoc, _ = ligolw_utils.load_fileobj(
        infile,
        contenthandler=ligolw_bayestar.LSCTablesAndSeriesContentHandler)

    snr_dict = ligolw_bayestar.snr_series_by_sngl_inspiral_id_for_xmldoc(xmldoc)

    coinc_and_sngl_inspirals = ligolw_bayestar.coinc_and_sngl_inspirals_for_xmldoc(xmldoc)

command.mkpath(opts.output)

if opts.condor_submit:
    if opts.coinc_event_id:
        raise ValueError('must not set --coinc-event-id with --condor-submit')
    cmd = ['condor_submit', 'accounting_group=ligo.dev.o3.cbc.pe.bayestar',
           'on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)',
           'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)',
           'request_memory = 1000 MB',
           'universe=vanilla', 'getenv=true', 'executable=' + sys.executable,
           'JobBatchName=BAYESTAR', 'environment="OMP_NUM_THREADS=1"',
           'error=' + os.path.join(opts.output, '$(CoincEventId).err'),
           'log=' + os.path.join(opts.output, '$(CoincEventId).log'),
           'arguments="-B ' + ' '.join(arg for arg in sys.argv
               if arg != '--condor-submit') + ' --coinc-event-id $(CoincEventId)"',
           '-append', 'queue CoincEventId in ' + ' '.join(
               str(coinc_event_id) for coinc_event_id in coinc_and_sngl_inspirals),
           '/dev/null']
    os.execvp('condor_submit', cmd)

if opts.psd_files:  # read pycbc psds here
    import lal
    from glue.segments import segment, segmentlist
    import h5py

    class psd_segment(segment):

        def __new__(cls, psd, *args):
            return segment.__new__(cls, *args)

        def __init__(self, psd, *args):
            self.psd = psd

    psdseglistdict = {}
    for psd_file in opts.psd_files:
        (ifo, group), = h5py.File(psd_file, 'r').items()
        psd = [group['psds'][str(i)] for i in range(len(group['psds']))]
        psdseglistdict[ifo] = segmentlist(
            psd_segment(*segargs) for segargs in zip(
                psd, group['start_time'], group['end_time']))

    def reference_psd_for_sngl(sngl):
        psd = psdseglistdict[sngl.ifo]
        try:
            psd = psd[psd.find(sngl.end)].psd
        except ValueError:
            raise ValueError(
                'No PSD found for detector {0} at GPS time {1}'.format(
                    sngl.ifo, sngl.end))

        flow = psd.file.attrs['low_frequency_cutoff']
        df = psd.attrs['delta_f']
        kmin = int(flow / df)

        fseries = lal.CreateREAL8FrequencySeries(
            'psd', 0, kmin * df, df,
            lal.StrainUnit**2 / lal.HertzUnit, len(psd.value) - kmin)
        fseries.data.data = psd.value[kmin:] / np.square(
            ligolw_bayestar.PYCBC_DYN_RANGE_FAC)

        return timing.InterpolatedPSD(
            filter.abscissa(fseries), fseries.data.data,
            f_high_truncate=opts.f_high_truncate)

    def reference_psds_for_sngls(sngl_inspirals):
        return [reference_psd_for_sngl(sngl) for sngl in sngl_inspirals]
else:
    reference_psd_filenames_by_process_id = ligolw_bayestar.psd_filenames_by_process_id_for_xmldoc(xmldoc)

    @memoized
    def reference_psds_for_filename(filename):
        xmldoc = ligolw_utils.load_filename(
            filename, contenthandler=lal.series.PSDContentHandler)
        psds = lal.series.read_psd_xmldoc(xmldoc, root_name=None)
        return {
            key: timing.InterpolatedPSD(filter.abscissa(psd), psd.data.data,
                f_high_truncate=opts.f_high_truncate)
            for key, psd in psds.items() if psd is not None}

    def reference_psd_for_ifo_and_filename(ifo, filename):
        return reference_psds_for_filename(filename)[ifo]

    def reference_psds_for_sngls(sngl_inspirals):
        return tuple(
            reference_psd_for_ifo_and_filename(sngl_inspiral.ifo,
            reference_psd_filenames_by_process_id[sngl_inspiral.process_id])
            for sngl_inspiral in sngl_inspirals)

count_sky_maps_failed = 0

# Loop over all coinc_event <-> sim_inspiral coincs.
if opts.coinc_event_id:
    coinc_and_sngl_inspirals = {
        coinc_event_id: coinc_and_sngl_inspirals[coinc_event_id]
        for coinc_event_id in opts.coinc_event_id}

for int_coinc_event_id, sngl_inspirals in coinc_and_sngl_inspirals.items():
    coinc_event_id = 'coinc_event:coinc_event_id:{}'.format(int_coinc_event_id)

    instruments = {sngl_inspiral.ifo for sngl_inspiral in sngl_inspirals}

    # Look up PSDs
    log.info('%s:reading PSDs', coinc_event_id)
    psds = reference_psds_for_sngls(sngl_inspirals)

    # Look up SNR time series
    try:
        snrs = [snr_dict[sngl.event_id] for sngl in sngl_inspirals]
    except KeyError:
        snrs = None

    # Loop over sky localization methods
    for method in opts.method:
        log.info("%s:method '%s':computing sky map", coinc_event_id, method)
        if opts.chain_dump:
            chain_dump = '%s.chain.npy' % int_coinc_event_id
        else:
            chain_dump = None
        try:
            sky_map = ligolw_sky_map(
                sngl_inspirals, opts.waveform, opts.f_low, opts.min_distance,
                opts.max_distance, opts.prior_distance_power, opts.cosmology,
                psds=psds, method=method, nside=opts.nside,
                chain_dump=chain_dump, phase_convention=opts.phase_convention,
                snr_series=snrs, enable_snr_series=opts.enable_snr_series)
            sky_map.meta['objid'] = coinc_event_id
        except (ArithmeticError, ValueError):
            log.exception(
                "%s:method '%s':sky localization failed",
                coinc_event_id, method)
            count_sky_maps_failed += 1
            if not opts.keep_going:
                raise
        else:
            log.info(
                "%s:method '%s':saving sky map",
                coinc_event_id, method)
            filename = '%d.%s.fits' % (int_coinc_event_id, method)
            fits.write_sky_map(
                os.path.join(opts.output, filename), sky_map, nest=True)


if count_sky_maps_failed > 0:
    raise RuntimeError("{0} sky map{1} did not converge".format(
        count_sky_maps_failed, 's' if count_sky_maps_failed > 1 else ''))
