# -*- coding: utf-8 -*-
#
# Copyright (C) 2018 Matthew Pitkin <matthew.pitkin@ligo.org>.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with with program; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA  02110-1301  USA

from __future__ import print_function, division

import os
import sys
import argparse

try:
    # function to get positions and velocities
    from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
    from astropy.time import Time, TimeDelta
    from astropy import constants as const
except ImportError:
    raise Exception("Could not import functions from astropy")

try:
    from lalapps import git_version
    __author__ = git_version.builder # person who built the code and therefore created the ephemeris!
    __date__ = git_version.date
    __id__ = git_version.id          # git tag of build
    __branch__ = git_version.branch  # git branch of build
except ImportError:
    __author__ = os.environ['USER']
    import datetime
    __date__ = datetime.datatime.now()
    __id__ = 'Unknown'
    __branch__ = 'Unknown'

HEADER = """\
# Build information for {exec}
# Author: {author}
# LALApps Commit ID: {gitid}
# LALApps Commit Date: {gitdate}
# LALApps Commit Branch: {gitbranch}
#
# Ephemeris creation command:-
#       {command}
#
# The JPL {ephem} ephemeris {ephemurl} has been used.
#
# This file consists of a header line containing:
#       GPS time of the start, interval between entries (secs), no. of entries
# Each entry consists of:
#       GPS time                Pos. x (lt sec)         Pos. y (lt sec)
#       Pos. z (lt sec)         Vel. x (lt sec/sec)     Vel. y (lt sec/sec)
#       Vel. z (lt sec/sec)     Acc. x (lt sec/sec^2)   Acc. y (lt sec/sec^2)
#       Acc. z (lt sec/sec^2)
"""

# set locations of JPL ephemeris files for downloading
EPH_URLS = {'DE436': 'ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de436.bsp',
            'DE435': 'ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de435.bsp',
            'DE432S': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp',
            'DE430': 'http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp',
            'DE421': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp',
            'DE414': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de414.bsp',
            'DE410': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de410.bsp',
            'DE406S': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de406s.bsp',
            'DE405': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp',
            'DE200': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de200.bsp'}

# code description
DESCRIPTION = """Create an ephemeris file for a given solar system body containing positions and \
velocities at a requested set of times.
"""

# set of allowed solar system bodies
BODIES = ['sun',
          'mercury',
          'venus',
          'earth-moon-barycenter',
          'earth-moon-barycentre', # allow UK spelling
          'earth',
          'moon',
          'mars',
          'jupiter',
          'saturn',
          'uranus',
          'neptune',
          'pluto']


if __name__=='__main__':
    parser = argparse.ArgumentParser( )
    parser.add_argument("-e", "--ephemeris", dest="ephemeris", required=True, help="Set the ephemeris to use (e.g. DE405)")
    parser.add_argument("-o", "--output-file", dest="output", default=None, help="Set the output file (defaults to stdout)")
    parser.add_argument("-g", "--gps-start", dest="gpsstart", type=float, default=None, help="Set the GPS time at which to start generating the ephemeris")
    parser.add_argument("-y", "--year-start", dest="yearstart", type=float, default=None, help="Set the (decimal) year at which to start generating the ephemeris")
    parser.add_argument("-n", "--num-years", dest="nyears", type=float, required=True, help="Set the number of years over which to generate the ephemeris")
    parser.add_argument("-i", "--interval", dest="interval", type=float, required=True, help="Set the time step (in hours) between successive output points")
    parser.add_argument("-t", "--target", dest="target", required=True, help="Set the solar system body to generate the ephemeris for")

    args = parser.parse_args()

    # check ephemeris is in our current list
    if args.ephemeris.upper() not in EPH_URLS.keys():
        print("Ephemeris '{}' is not allowed, use one of: {}".format(args.ephemeris, EPH_URLS.keys()))
        sys.exit(1)
    else:
        ephemfile = EPH_URLS[args.ephemeris.upper()]

    # check that the body is in our current list
    if args.target.lower() not in BODIES:
        print("Target body '{}' is not in the allowed list: {}".format(args.target, BODIES))
        sys.exit(1)
    else:
        body = args.target.lower() if args.target.lower() != 'earth-moon-barycentre' else 'earth-moon-barycenter'

    # set the ephemeris file
    solar_system_ephemeris.set(ephemfile)

    # set the start time
    if args.gpsstart is not None and args.yearstart is not None:
        print("Specify either '--gps-start' or '--year-start', but not both")
        sys.exit(1)
    elif args.gpsstart is not None:
        try:
            starttime = Time(args.gpsstart, format='gps', scale='utc')
        except ValueError:
            Exception("Could not parse start GPS time: {}".format(args.gpsstart))
    else:
        try:
            starttime = Time(args.yearstart, format='decimalyear', scale='utc')
        except ValueError:
            Exception("Could not parse start year: {}".format(args.yearstart))

    # set the time step
    try:
        dt = TimeDelta(args.interval*3600., format='sec')
    except ValueError:
        Exception("Could not parse time interval: {}".format(args.interval))

    # set the end time
    try:
        endtime = Time(starttime.decimalyear + args.nyears, format='decimalyear', scale='utc') + dt
    except ValueError:
        Exception("Could not parse total timespan")

    pos = []
    vel = []
    acc = []

    # get positions, velocities and accelerations
    curtime = starttime
    while curtime <= endtime:
        tpos, tvel = get_body_barycentric_posvel(body, curtime)

        # convert positions to light seconds
        pos.append(tpos.xyz.to('m')/const.c)

        # convert velocities to light seconds per second
        vel.append(tvel.xyz.to('m/s')/const.c)

        # calculate accelerations (using velocities +/- dt/2 around the current time)
        ctminus = curtime - dt/2.
        _, mvel = get_body_barycentric_posvel(body, ctminus)
        ctplus = curtime + dt/2.
        _, pvel = get_body_barycentric_posvel(body, ctplus)
        acc.append(((pvel.xyz.to('m/s')-mvel.xyz.to('m/s'))/const.c)/dt.to('s'))

        curtime += dt

    # set output header
    headdic = {}
    headdic['exec'] = os.path.basename(sys.argv[0])
    headdic['author'] = __author__
    headdic['gitid'] = __id__
    headdic['gitbranch'] = __branch__
    headdic['gitdate'] = __date__
    headdic['command'] = ' '.join([os.path.basename(sys.argv[0])]+sys.argv[1:])
    headdic['ephem'] = args.ephemeris.upper()
    headdic['ephemurl'] = ephemfile

    outfile = args.output
    if outfile is None:
        fp = sys.stdout
    else:
        # gzip if extension ends with '.gz'
        if outfile[-3:] == '.gz':
            try:
                import gzip
                fp = gzip.open(outfile, 'wb')
            except IOError:
                Exception("Problem opening gzip output file '{}'".format(outfile))
        else:
            try:
                fp = open(outfile, 'w')
            except IOError:
                Exception("Problem opening output file '{}'".format(outfile))

    # write out header
    fp.write(HEADER.format(**headdic))

    # write out start time (GPS), time interval (secs), number of entries
    fp.write('{}\t{}\t{}\n'.format(int(starttime.gps), dt, len(pos)))

    curtime = starttime
    for tpos, tvel, tacc in zip(pos, vel, acc):
        fp.write('{0:.6f}\t'.format(curtime.gps))
        fp.write('{0:.16e}\t{1:.16e}\t{2:.16e}\t'.format(tpos[0].value, tpos[1].value, tpos[2].value))
        fp.write('{0:.16e}\t{1:.16e}\t{2:.16e}\t'.format(tvel[0].value, tvel[1].value, tvel[2].value))
        fp.write('{0:.16e}\t{1:.16e}\t{2:.16e}\n'.format(tacc[0].value, tacc[1].value, tacc[2].value))
        curtime += dt

    if outfile is not None:
        fp.close()
