#
# Copyright (C) 2012-2016  Leo Singer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Plotting classes and methods
"""
from . import allsky, healpix, marker, poly, pp
from .allsky import *
from .healpix import *
from .marker import *
from .poly import *
from .pp import *
__all__ = (allsky.__all__ + healpix.__all__ + marker.__all__ + poly.__all__
           + pp.__all__)
