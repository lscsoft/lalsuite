#
# Copyright (C) 2013-2015  Leo Singer
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
r"""
Listen for new events from LVAlert and perform sky localization.

`bayestar_localize_lvalert` supports two modes of operation. You can
explicitly specify the GraceDb ID on the command line, as in:

    $ bayestar_localize_lvalert T90713

Or, `bayetar_localize_lvalert` can accept an LVAlert XML packet read from
stdin. This is handy for quickly setting up automatic processing of events in
response to LVAlert notifications. To do this, first put your LVAlert login
credentials in the file `~/.netrc` in your home directory:

    machine lvalert.cgca.uwm.edu login albert.einstein password ligorocks

replacing `albert.einstein` and `ligorocks` with your actual username and
password. Then, subscribe to events of types that you are interested in
(probably `cbc_lowmass` or `test_lowmass`):

    $ lvalert_admin --subscribe --node cbc_lowmass
    $ lvalert_admin --subscribe --node test_lowmass

Create a configuration file that will tell `lvalert_listen` what to do in
response to those event types. For example, you might create a file called
`lvalert_listen.ini` with the following contents:

    [test_lowmass]
    executable = bayestar_localize_lvalert
    [cbc_lowmass]
    executable = bayestar_localize_lvalert

Finally, start `lvalert_listen`:

    $ lvalert_listen --config-file lvalert_listen.ini
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


#
# Command line interface
#

from lalinference.bayestar import command
import logging
import ligo.lvalert.utils
import sys
import os

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('BAYESTAR')

methods = '''
    toa_phoa_snr
    toa_snr_mcmc
    toa_phoa_snr_mcmc
    toa_snr_mcmc_kde
    toa_phoa_snr_mcmc_kde
    '''.split()
default_method = 'toa_phoa_snr'
command.skymap_parser.add_argument(
    '--method', choices=methods, default=default_method,
    help='Sky localization method [default: %(default)s]')
parser = command.ArgumentParser(
    parents=[command.waveform_parser, command.prior_parser, command.skymap_parser])
parser.add_argument('-N', '--dry-run', default=False, action='store_true',
    help='Dry run; do not update GraceDB entry [default: %(default)s]')
parser.add_argument('-o', '--output', metavar='FILE.fits[.gz]',
    default='bayestar.fits.gz',
    help='Name for uploaded file [default: %(default)s]')
parser.add_argument('graceid', metavar='G123456', nargs='?',
    help='Run on this GraceDB ID [default: listen to lvalert]')
opts = parser.parse_args()

if opts.graceid is None:
    # No GraceDB ID; read LVAlert data from stdin
    lvadata = ligo.lvalert.utils.get_LVAdata_from_stdin(sys.stdin, as_dict=True)
    log.info("received lvalert event of type='%s' for uid='%s' and file='%s'",
        lvadata['alert_type'], lvadata['uid'], lvadata['file'])
    if lvadata['alert_type'] == 'update' and lvadata['file'] == 'psd.xml.gz':
        graceid = lvadata['uid']
    else:
        log.info('ignoring')
        raise SystemExit
else:
    # One command line argument; manual start from GraceDB id
    graceid = opts.graceid


#
# Hook logging into GraceDb
#

import ligo.gracedb.logging
import ligo.gracedb.rest

# Fire up a GraceDb client
# FIXME: Mimic the behavior of the GraceDb command line client, where the
# environment variable GRACEDB_SERVICE_URL overrides the default service URL.
# It would be nice to get this behavior into the gracedb package itself.
gracedb = ligo.gracedb.rest.GraceDb(
    os.environ.get('GRACEDB_SERVICE_URL',
    ligo.gracedb.rest.DEFAULT_SERVICE_URL))

# Send log messages to GraceDb too
if not opts.dry_run:
    handler = ligo.gracedb.logging.GraceDbLogHandler(gracedb, graceid)
    handler.setLevel(logging.INFO)
    logging.root.addHandler(handler)


#
# Perform sky localization
#

import os
import shutil
import tempfile
from lalinference.bayestar import distance
from lalinference.bayestar.ligolw_sky_map import gracedb_sky_map
from lalinference import fits

# A little bit of Cylon humor
log.info('by your command...')

try:
    # download coinc.xml
    coinc_file = gracedb.files(graceid, "coinc.xml")

    # download psd.xml.gz
    psd_file = gracedb.files(graceid, "psd.xml.gz")

    if opts.chain_dump:
        chain_dump = opts.output.replace('.fits.gz', '').replace('.fits', '') + '.chain.npy'
    else:
        chain_dump = None

    # perform sky localization
    log.info("starting sky localization")
    sky_map, epoch, elapsed_time, instruments = gracedb_sky_map(
        coinc_file, psd_file, opts.waveform, opts.f_low,
        opts.min_distance, opts.max_distance, opts.prior_distance_power,
        phase_convention=opts.phase_convention, nside=opts.nside,
        f_high_truncate=opts.f_high_truncate,
        method=opts.method, chain_dump=chain_dump)
    prob, distmu, distsigma, _ = sky_map
    distmean, diststd = distance.parameters_to_marginal_moments(
        prob, distmu, distsigma)
    log.info("sky localization complete")

    # upload FITS file
    fitsdir = tempfile.mkdtemp()
    try:
        fitspath = os.path.join(fitsdir, opts.output)
        fits.write_sky_map(fitspath, sky_map, gps_time=float(epoch),
            creator=parser.prog, objid=str(graceid),
            url='https://gracedb.ligo.org/events/{0}'.format(graceid),
            runtime=elapsed_time, instruments=instruments,
            distmean=distmean, diststd=diststd,
            origin='LIGO/Virgo', nest=True)
        if not opts.dry_run:
            gracedb.writeLog(graceid, "INFO:BAYESTAR:uploaded sky map",
                filename=fitspath, tagname=("sky_loc", "lvem"))
        else:
            command.rename(fitspath, os.path.join('.', opts.output))
    finally:
        shutil.rmtree(fitsdir)
except:
    # Produce log message for any otherwise uncaught exception
    log.exception("sky localization failed")
    # Then re-raise the exception
    raise
