# Copyright (C) 2010  Peter Couvares
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


import os
import errno


"""
Check for pid existence efficiently and reliably.  (Using the null
signal is faster and more unix-portable than looking in /proc).

Inspired by Larz Wirzenius
<http://stackoverflow.com/questions/1005972>
"""
def pid_exists(pid):
    """ Returns true if the given pid exists, false otherwise. """
    try:
        # signal 0 is harmless and can be safely used to probe pid existence
        # faster and more unix-portable than looking in /proc
        os.kill(pid, 0)
    except OSError as e:
        # "permission denied" proves existence; otherwise, no such pid
        return e.errno == errno.EPERM
    else:
        return True


"""
Performs the equivalent of "mkdir -p", creating any intermediate
directories needed to create the leaf directory -- but unlike
os.makedirs(), produces no error if the path already exists.

Inspired by Christos Georgiou
<http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python>
"""
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else: raise
    return path
