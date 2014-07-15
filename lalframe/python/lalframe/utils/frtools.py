# Copyright (C) 2013 Duncan Macleod
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
# Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

"""Native python versions of LIGOTools and FrXXX frame utilities utilities
"""

import lalframe

from lalframe import git_version
__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date


def get_channels(framefile):
    """Return a list of all channels written into the given GWF-format
    framefile

    @param framefile
        `string` path to input frame file

    @returns a list of all channels present in the table of contents

    Example:

    \code
    >>> get_channels('H-H1_LDAS_C02_L2-9668/H-H1_LDAS_C02_L2-966802176-128.gwf')
    ['H1:LSC-DATA_QUALITY_VECTOR', 'H1:IFO-SV_STATE_VECTOR', 'H1:LDAS-STRAIN']
    \endcode
    """
    frfile = lalframe.FrameUFrFileOpen(framefile, "r")
    frtoc = lalframe.FrameUFrTOCRead(frfile)
    nadc = lalframe.FrameUFrTOCQueryAdcN(frtoc)
    nproc = lalframe.FrameUFrTOCQueryProcN(frtoc)
    nsim = lalframe.FrameUFrTOCQuerySimN(frtoc)
    adcchannels = [lalframe.FrameUFrTOCQueryAdcName(frtoc, i) for
                   i in range(nadc)]
    procchannels = [lalframe.FrameUFrTOCQuerySimName(frtoc, i) for
                    i in range(nsim)]
    simchannels = [lalframe.FrameUFrTOCQueryProcName(frtoc, i) for
                   i in range(nproc)]
    return sorted(adcchannels + procchannels + simchannels)
