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
Collection of Python decorators.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"
__all__ = ('memoized', 'with_numpy_random_seed', 'as_dict')


from functools import wraps
from collections import Hashable
from astropy.utils.misc import NumpyRNGContext


def memoized(func):
    """Memoize a function or class by caching its return values for any given
    arguments."""
    cache = {}

    # FIXME: In Python 3.4, use inspect.getcallargs to bind function arguments.
    # This will allow us to handle default arguments better.
    # (Though inspect.getcallargs was added in Python 2.7, it won't work with
    # built-in functions or functions from C extensions until Python 3.4.
    # See https://bugs.python.org/issue17481.)

    @wraps(func)
    def memo(*args, **kwargs):
        # Create a key out of the arguments.
        key = (args, frozenset(kwargs.items()))

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


def with_numpy_random_seed(func, seed=0):
    """Decorate a function so that it is called with a pre-defined random seed.
    The random seed is restored when the function returns."""

    @wraps(func)
    def wrapped_func(*args, **kwargs):
        with NumpyRNGContext(seed):
            ret = func(*args, **kwargs)
        return ret

    return wrapped_func


def as_dict(func):
    """Decorate a generator function to turn its output into a dict."""

    @wraps(func)
    def wrapped_func(*args, **kwargs):
        return dict(func(*args, **kwargs))

    return wrapped_func
