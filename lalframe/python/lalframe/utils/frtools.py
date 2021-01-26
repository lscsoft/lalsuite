# -*- coding: utf-8 -*-
# Copyright (C) 2013,2019 Duncan Macleod
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

"""Native python versions of LIGOTools and FrXXX frame utilities
"""

import lalframe

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"


def get_channels(framefile):
    """Return a list of all channels written into the given GWF-format
    framefile

    @param framefile
        `string` path to input frame file

    @return A list of channel names

    Example:

    \code
    >>> get_channels('H-H1_LDAS_C02_L2-9668/H-H1_LDAS_C02_L2-966802176-128.gwf')
    ['H1:LSC-DATA_QUALITY_VECTOR', 'H1:IFO-SV_STATE_VECTOR', 'H1:LDAS-STRAIN']
    \endcode
    """
    return sorted(iter_channels(framefile))


def iter_channels(framefile):
    """Yield channels from the table-of-contents of the given GWF file

    @param framefile
        `str` path to input frame file
    @return An iterator of channel names
    """
    frfile = lalframe.FrameUFrFileOpen(str(framefile), "r")
    frtoc = lalframe.FrameUFrTOCRead(frfile)
    nsim = lalframe.FrameUFrTOCQuerySimN(frtoc)
    for i in range(lalframe.FrameUFrTOCQueryAdcN(frtoc)):
        yield lalframe.FrameUFrTOCQueryAdcName(frtoc, i)
    for i in range(lalframe.FrameUFrTOCQueryProcN(frtoc)):
        yield lalframe.FrameUFrTOCQueryProcName(frtoc, i)
    for i in range(lalframe.FrameUFrTOCQuerySimN(frtoc)):
        yield lalframe.FrameUFrTOCQuerySimName(frtoc, i)
