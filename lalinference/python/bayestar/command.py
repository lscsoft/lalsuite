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
"""
Functions that support the command line interface.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


from optparse import IndentedHelpFormatter
import glob
import itertools


class NewlinePreservingHelpFormatter(IndentedHelpFormatter):
    """A help formatter for optparse that preserves paragraphs and bulleted
    lists whose lines start with a whitespace character."""

    def _format_text(self, text):
        __doc__ = IndentedHelpFormatter._format_text
        return "\n\n".join(
            t if len(t) == 0 or t[0].isspace()
            else IndentedHelpFormatter._format_text(self, t)
            for t in text.split("\n\n")
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


def sqlite3_connect_nocreate(dbfilename):
    """Open an SQLite database, or fail if it does not exist.
    FIXME: use SQLite URI when we drop support for Python < 3.4.
    See: https://docs.python.org/3.4/whatsnew/3.4.html#sqlite3"""
    import sqlite3
    with open(dbfilename, 'rb') as testfile:
        pass
    return sqlite3.connect(dbfilename)
