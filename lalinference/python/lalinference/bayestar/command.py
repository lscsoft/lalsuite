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
"""
Functions that support the command line interface.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


import argparse
from optparse import IndentedHelpFormatter
import glob
import inspect
import itertools
import os
import sys


class NewlinePreservingHelpFormatter(IndentedHelpFormatter):
    """A help formatter for optparse that preserves paragraphs and bulleted
    lists whose lines start with a whitespace character."""

    def _format_text(self, text):
        __doc__ = IndentedHelpFormatter._format_text
        return "\n\n".join(
            t if len(t) == 0 or t[0].isspace()
            else IndentedHelpFormatter._format_text(self, t)
            for t in text.strip().split("\n\n")
        )


def check_required_arguments(parser, opts, *keys):
    """Raise an error if any of the specified command-line arguments are missing."""
    for key in keys:
        if getattr(opts, key) is None:
            parser.error("Missing required argument: --" + key.replace("_", "-"))


def get_input_filename(parser, args):
    """Determine name of input: either the sole positional command line argument,
    or /dev/stdin."""
    if len(args) == 0:
        infilename = '/dev/stdin'
    elif len(args) == 1:
        infilename = args[0]
    else:
        parser.error("Too many command line arguments.")
    return infilename


def chainglob(patterns):
    """Generate a list of all files matching a list of globs."""
    return itertools.chain.from_iterable(glob.iglob(s) for s in patterns)


waveform_parser = argparse.ArgumentParser(add_help=False)
group = waveform_parser.add_argument_group(
    'waveform options', 'Options that affect template waveform generation')
group.add_argument('--f-low', type=float, metavar='Hz', default=10,
    help='Low frequency cutoff [default: %(default)s]')
group.add_argument('--waveform', default='o1-uberbank',
    help='Template waveform approximant (e.g., TaylorF2threePointFivePN) '
    '[default: O1 uberbank mass-dependent waveform]')
del group


prior_parser = argparse.ArgumentParser(add_help=False)
group = prior_parser.add_argument_group(
    'prior options', 'Options that affect the BAYESTAR prior')
group.add_argument('--min-distance', type=float, metavar='Mpc',
    help='Minimum distance of prior in megaparsecs '
    '[default: infer from effective distance]')
group.add_argument('--max-distance', type=float, metavar='Mpc',
    help='Maximum distance of prior in megaparsecs '
    '[default: infer from effective distance]')
group.add_argument('--prior-distance-power', type=int, metavar='-1|2',
    default=2, help='Distance prior '
    '[-1 for uniform in log, 2 for uniform in volume, default: %(default)s]')
del group


class ArgumentParser(argparse.ArgumentParser):
    def __init__(self,
                 prog=None,
                 usage=None,
                 description=None,
                 epilog=None,
                 parents=[],
                 formatter_class=argparse.RawDescriptionHelpFormatter,
                 prefix_chars='-',
                 fromfile_prefix_chars=None,
                 argument_default=None,
                 conflict_handler='error',
                 add_help=True):
        if prog is None:
            prog = os.path.basename(sys.argv[0]).replace('.py', '')
        if description is None:
            parent_frame = inspect.currentframe().f_back
            description = parent_frame.f_locals.get('__doc__', None)
        super(ArgumentParser, self).__init__(
                 prog=prog,
                 usage=usage,
                 description=description,
                 epilog=epilog,
                 parents=parents,
                 formatter_class=argparse.RawDescriptionHelpFormatter,
                 prefix_chars=prefix_chars,
                 fromfile_prefix_chars=fromfile_prefix_chars,
                 argument_default=argument_default,
                 conflict_handler=conflict_handler,
                 add_help=add_help)


class SQLiteType(argparse.FileType):
    """Open an SQLite database, or fail if it does not exist.
    FIXME: use SQLite URI when we drop support for Python < 3.4.
    See: https://docs.python.org/3.4/whatsnew/3.4.html#sqlite3"""

    def __init__(self, mode='r'):
        super(SQLiteType, self).__init__(mode + 'b')

    def __call__(self, string):
        if string == '-':
            raise argparse.ArgumentTypeError(
                'Cannot open stdin/stdout as an SQLite database')
        with super(SQLiteType, self).__call__(string):
            import sqlite3
            return sqlite3.connect(string)


def sqlite_get_filename(connection):
    """Get the name of the file associated with an SQLite connection"""
    result = connection.execute('pragma database_list').fetchall()
    try:
        (_, _, filename), = result
    except ValueError:
        raise RuntimeError('Expected exactly one attached database')
    return filename
