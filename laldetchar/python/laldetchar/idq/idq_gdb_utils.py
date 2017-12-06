# Copyright (C) 2014 Reed Essick
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## addtogroup pkg_py_laldetchar_idq
## Synopsis
# ~~~
# from laldetchar.idq import idq_gdb_utils
# ~~~
# \author Reed Essick (<reed.essick@ligo.org>)

"""Module with utility functions used for generating iDQ input to GraceDB.
"""

from laldetchar import git_version

__author__ = 'Reed Essick <reed.essick@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

## addtogroup pkg_py_laldetchar_idq_auxmvc
# @{

# Utility functions used for generating iDQ input to GraceDB.

import numpy as np
import re as re
from laldetchar.idq import event
from laldetchar.idq import idq


def combine_ts(filenames):
    """ 
....combine multiple files into a single time-series. Assumes filenames have the standard LIGO naming covention: *-start-dur.suffix
....Also assumes that filenames are sorted into chronological order

....returns lists of arrays, with each array consisting of only contiguous data
........return timeseries, times
...."""

    t = np.array([])  # the array storing continuous data
    ts = np.array([])
    times = []  # a list that will contain stretches of continuous data
    timeseries = []

    matchfile = re.compile('.*-([0-9]*)-([0-9]*).*$')

    end = False
    for filename in filenames:
        m = matchfile.match(filename)
        (_start, _dur) = (int(m.group(1)), int(m.group(2)))

        # ## check to see if we have continuous data

        if not end or end == _start:  # beginning of data
            end = _start + _dur

            _file = event.gzopen(filename)
            _ts = np.load(_file)
            _file.close()

            ts = np.concatenate((ts, _ts))
            t = np.concatenate((t, np.arange(_start, _start + _dur, 1.0
                               * _dur / len(_ts))))
        else:

            # gap in the data!

            times.append(t)  # put old continuous data into lists
            timeseries.append(ts)

            _file = event.gzopen(filename)  # start new continuous data
            ts = np.load(_file)
            _file.close()
            t = np.arange(_start, _start + _dur, 1.0 * _dur / len(ts))
            end = _start + _dur

    times.append(t)
    timeseries.append(ts)

    return (times, timeseries)


def stats_ts(ts):
    """ 
....compute basic statistics about ts 

....return min(ts), max(ts), mean(ts), stdv(ts)
...."""

    return (np.min(ts), np.max(ts), np.mean(ts), np.var(ts) ** 0.5)



def execute_gdb_timeseries(
    gps_start,
    gps_end,
    gps,
    gracedb_id,
    ifo,
    classifier,
    cp,
    input_dir,
    exec_prog,
    usertag='',
    gch_xml=[],
    cln_xml=[]):
    """ Function that sets up and runs idq-gdb-timeseries script as one of the tasks of idq-gdb-processor."""
    # form the command line
    cmd_line = [exec_prog, '-s', gps_start, '-e', gps_end, '--gps', gps,\
        '-g', gracedb_id, '--ifo', ifo, '-c', classifier, '-i', input_dir,\
        '-t', usertag]
	
    # add extra options from config file
    if cp.has_option("general","gdb_url"):
        cmd_line += ["--gdb-url", cp.get("general","gdb_url")]
    for gch in gch_xml:
        cmd_line += ["--gch-xml", gch]
    for cln in cln_xml:
        cmd_line += ["--cln-xml", cln]

    for (option,value) in cp.items('gdb-time-series'):
        cmd_line.extend([option, value])
	print cmd_line
    exit_status = idq.submit_command(cmd_line, 'gdb_timeseries', verbose=True)
	
    return exit_status
	
	
	
def execute_gdb_glitch_tables(
    gps_start,
    gps_end,
    gracedb_id,
    ifo,
    classifier,
    cp,
    input_dir,
    exec_prog,
    usertag=''):
    """ Function that sets up and runs idq-gdb-timeseries script as one of the tasks of idq-gdb-processor."""
    # form the command line
    cmd_line = [exec_prog, '-s', gps_start, '-e', gps_end, '-g', gracedb_id,\
        '--ifo', ifo, '-c', classifier, '-i', input_dir, '-t', usertag]
	
    # add extra options from config file
    if cp.has_option("general","gdb_url"):
        cmd_line += ["--gdb-url", cp.get("general","gdb_url")]
    for (option,value) in cp.items('gdb-glitch-tables'):
        cmd_line.extend([option, value])
	print cmd_line
    exit_status = idq.submit_command(cmd_line, 'gdb_glitch_tables', verbose=True)
	
    return exit_status	
##@}

