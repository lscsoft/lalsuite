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
Match sky maps with injections in an inspinjfind-style sqlite database and print
summary values for each sky map:

 * event ID
 * false alarm rate
 * searched area
 * searched posterior probability
 * angle between true sky location and maximum a posteriori estimate
 * runtime in seconds
 * (optional) areas of and numbers of modes within specified probability contours

The filenames of the sky maps may be provided as positional command line
arguments, and may also be provided as globs (such as '*.fits.gz').
"""
from __future__ import division
from __future__ import print_function
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Command line interface.
import argparse
from lalinference.bayestar import command
if __name__ == '__main__':
    parser = command.ArgumentParser()
    parser.add_argument(
        '-o', '--output', metavar='OUT.dat',
        type=argparse.FileType('w'), default='-',
        help='Name of output file [default: stdout]')
    parser.add_argument(
        '-j', '--jobs', type=int, default=1, const=None, nargs='?',
        help='Number of threads [default: %(default)s]')
    parser.add_argument(
        '-p', '--contour', default=[], nargs='+', type=float, metavar='PERCENT',
        help='Report the area of the smallest contour and the number of modes '
        'containing this much probability.')
    parser.add_argument(
        '-a', '--area', default=[], nargs='+', type=float, metavar='DEG2',
        help='Report the largest probability contained within any region '
        'of this area in square degrees. Can be repeated multiple times.')
    parser.add_argument(
        '--modes', default=False, action='store_true',
        help='Compute number of disjoint modes [default: %(default)s]')
    parser.add_argument(
        'db', type=command.SQLiteType('r'), metavar='DB.sqlite',
        help='Input SQLite database from search pipeline')
    parser.add_argument(
        'fitsfileglobs', metavar='GLOB.fits[.gz]', nargs='+',
        help='Input FITS filenames and/or globs')
    opts = parser.parse_args()


# Imports.
import sqlite3
from lalinference.io import fits
from lalinference.bayestar import postprocess


def startup(dbfilename, opts_contour, opts_modes, opts_area):
    global db, contours, modes, areas
    db = sqlite3.connect(dbfilename)
    contours = opts_contour
    modes = opts_modes
    areas = opts_area


def process(fitsfilename):
    sky_map, metadata = fits.read_sky_map(fitsfilename, nest=None)

    coinc_event_id = metadata['objid']
    try:
        runtime = metadata['runtime']
    except KeyError:
        runtime = float('nan')

    row = db.execute("""
        SELECT DISTINCT sim.simulation_id AS simulation_id, sim.longitude AS ra, sim.latitude AS dec,
        ci.combined_far AS far, ci.snr AS snr
        FROM coinc_event_map AS cem1 INNER JOIN coinc_event_map AS cem2
        ON (cem1.coinc_event_id = cem2.coinc_event_id)
        INNER JOIN sim_inspiral AS sim ON (cem1.event_id = sim.simulation_id)
        INNER JOIN coinc_inspiral AS ci ON (cem2.event_id = ci.coinc_event_id)
        WHERE cem1.table_name = 'sim_inspiral'
        AND cem2.table_name = 'coinc_event' AND cem2.event_id = ?""",
        (coinc_event_id,)).fetchone()
    if row is None:
        raise ValueError(
            "No database record found for event '{0}' in '{1}'".format(
            coinc_event_id, command.sqlite_get_filename(db)))
    simulation_id, true_ra, true_dec, far, snr = row
    searched_area, searched_prob, offset, searched_modes, contour_areas, area_probs, contour_modes = postprocess.find_injection(
        sky_map, true_ra, true_dec, contours=[0.01 * p for p in contours],
        areas=areas, modes=modes, nest=metadata['nest'])

    if snr is None:
        snr = float('nan')
    if far is None:
        far = float('nan')
    distmean = metadata.get('distmean', float('nan'))
    diststd = metadata.get('diststd', float('nan'))

    ret = [coinc_event_id, simulation_id, far, snr, searched_area, searched_prob, offset, runtime, distmean, diststd] + contour_areas + area_probs
    if modes:
        ret += [searched_modes] + contour_modes
    return ret


if __name__ == '__main__':
    from glue.text_progress_bar import ProgressBar
    progress = ProgressBar()

    db = opts.db
    contours = opts.contour
    modes = opts.modes
    areas = opts.area

    progress.update(-1, 'spawning workers')
    if opts.jobs == 1:
        from six.moves import map
    else:
        try:
            from emcee.interruptible_pool import InterruptiblePool as Pool
        except ImportError:
            from multiprocessing import Pool
        map = Pool(
            opts.jobs, startup,
            (command.sqlite_get_filename(db), contours, modes, areas)
            ).imap_unordered

    progress.update(-1, 'obtaining filenames of sky maps')
    fitsfilenames = tuple(command.chainglob(opts.fitsfileglobs))

    colnames = (
        ['coinc_event_id', 'simulation_id', 'far', 'snr', 'searched_area',
        'searched_prob', 'offset', 'runtime', 'distmean', 'diststd'] +
        ["area({0:g})".format(p) for p in contours] +
        ["prob({0:g})".format(a) for a in areas])
    if modes:
        colnames += ['searched_modes'] + ["modes({0:g})".format(p) for p in contours]
    print(*colnames, sep="\t", file=opts.output)

    count_records = 0
    progress.max = len(fitsfilenames)
    for record in map(process, fitsfilenames):
        count_records += 1
        progress.update(count_records, record[0])
        print(*record, sep="\t", file=opts.output)
