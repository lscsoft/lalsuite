# Copyright (C) 2017  Leo Singer
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
from . import base, detector_disabled, ligolw, gracedb, hdf, magic
from .base import *
from .detector_disabled import *
from .ligolw import *
from .gracedb import *
from .hdf import *
from .magic import *
from .magic import open
__all__ = (base.__all__ + detector_disabled.__all__ + ligolw.__all__ +
           gracedb.__all__ + hdf.__all__ + magic.__all__)
