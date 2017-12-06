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
Collection of Python decorators.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


def memoized(func):
    """Memoize a function or class by caching its return values for any given
    arguments."""
    from functools import wraps
    from collections import Hashable

    cache = {}

    # FIXME: In Python 2.7, use inspect.getcallargs to bind function arguments.
    # This will allow us to handle default arguments better.

    @wraps(func)
    def memo(*args, **kwargs):
        # Create a key out of the arguments.
        key = (args, frozenset(kwargs.iteritems()))

        if isinstance(args, Hashable): # The key is immutable.
            try:
                # Look up the return value for these arguments.
                ret = cache[key]
            except KeyError:
                # Not found; invoke function and store return value.
                ret = cache[key] = func(*args, **kwargs)
        else: # The key is mutable. We can't cache it.
            ret = func(*args, **kwargs)

        # Done!
        return ret

    # Return wrapped function.
    return memo
