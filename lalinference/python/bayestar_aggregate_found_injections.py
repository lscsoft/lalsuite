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
from __future__ import print_function
__doc__ = """
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
__author__ = "Leo Singer <leo.singer@ligo.org>"


from lalinference.bayestar import command

if __name__ == '__main__':
    # Command line interface.
    from optparse import Option, OptionParser

    parser = OptionParser(
        formatter=command.NewlinePreservingHelpFormatter(),
        description=__doc__,
        usage="%prog [-o OUTPUT] DATABASE.sqlite FILE1.fits[.gz] FILE2.fits[.gz] ...",
        option_list=[
            Option("-o", "--output", default="/dev/stdout",
                help="Name of output file [default: %default]"),
            Option("-j", "--jobs", default=1, type=int,
                help="Number of threads [default: %default]"),
            Option("-p", "--contour", default=[], action="append",
                type=float, metavar="PERCENT",
                help="Report the area of the smallest contour and the number of modes"
                + "containing this much probability. Can be repeated mulitple times"),
            Option("--modes", default=False, action="store_true",
                help="Compute number of disjoint modes [default: %default]")
        ]
    )
    opts, args = parser.parse_args()

    try:
        dbfilename = args[0]
        fitsfileglobs = args[1:]
    except IndexError:
        parser.error("not enough command line arguments")
    if opts.jobs < 1:
        parser.error("invalid value for -j, --jobs: must be >= 1")

    outfile = open(opts.output, "w")


# Imports.
import sqlite3
from lalinference import fits
from lalinference.bayestar import postprocess


def startup(dbfilename, opts_contour, opts_modes):
    global db, contours, modes
    db = command.sqlite3_connect_nocreate(dbfilename)
    contours = opts_contour
    modes = opts_modes


def process(fitsfilename):
    sky_map, metadata = fits.read_sky_map(fitsfilename)

    coinc_event_id = metadata['objid']
    try:
        runtime = metadata['runtime']
    except KeyError:
        runtime = float('nan')

    simulation_id, true_ra, true_dec, far, snr = db.execute("""
        SELECT DISTINCT sim.simulation_id AS simulation_id, sim.longitude AS ra, sim.latitude AS dec,
        ci.combined_far AS far, ci.snr AS snr
        FROM coinc_event_map AS cem1 INNER JOIN coinc_event_map AS cem2
        ON (cem1.coinc_event_id = cem2.coinc_event_id)
        INNER JOIN sim_inspiral AS sim ON (cem1.event_id = sim.simulation_id)
        INNER JOIN coinc_inspiral AS ci ON (cem2.event_id = ci.coinc_event_id)
        WHERE cem1.table_name = 'sim_inspiral'
        AND cem2.table_name = 'coinc_event' AND cem2.event_id = ?""",
        (coinc_event_id,)).fetchone()
    searched_area, searched_prob, offset, searched_modes, contour_areas, contour_modes = postprocess.find_injection(
        sky_map, true_ra, true_dec, contours=[0.01 * p for p in contours], modes=modes)

    if snr is None:
        snr = float('nan')
    if far is None:
        far = float('nan')

    ret = [coinc_event_id, simulation_id, far, snr, searched_area, searched_prob, offset, runtime] + contour_areas
    if modes:
        ret += [searched_modes] + contour_modes
    return ret


if __name__ == '__main__':
    from glue.text_progress_bar import ProgressBar
    progress = ProgressBar()

    progress.update(-1, 'spawning {0} workers'.format(opts.jobs))
    startupargs = (dbfilename, opts.contour, opts.modes)
    if opts.jobs == 1:
        from itertools import imap
    else:
        import multiprocessing
        imap = multiprocessing.Pool(opts.jobs, startup, startupargs).imap_unordered
    startup(*startupargs)

    progress.update(-1, 'obtaining filenames of sky maps')
    fitsfilenames = tuple(command.chainglob(fitsfileglobs))

    colnames = ['coinc_event_id', 'simulation_id', 'far', 'snr', 'searched_area',
        'searched_prob', 'offset', 'runtime'] + ["area({0:g})".format(p)
        for p in contours]
    if modes:
        colnames += ['searched_modes'] + ["modes({0:g})".format(p) for p in contours]
    print(*colnames, sep="\t", file=outfile)

    count_records = 0
    progress.max = len(fitsfilenames)
    for record in imap(process, fitsfilenames):
        count_records += 1
        progress.update(count_records, record[0])
        print(*record, sep="\t", file=outfile)
