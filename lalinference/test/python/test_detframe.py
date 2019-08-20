# Test program for functions mapping between Equatorial and Detector-based coordinate systems
# Copyright (C) 2016 John Veitch <john.veitch@ligo.org>
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

import lal
import lalinference as li
import numpy as np
import sys

Ntest=1000

LHO=lal.CachedDetectors[lal.LALDetectorIndexLHODIFF]
LLO=lal.CachedDetectors[lal.LALDetectorIndexLLODIFF]

def test_invertable(ra,dec,time,tolerance=[1e-6,1e-5,1e-5]):
    forward = li.EquatorialToDetFrame(LHO,LLO,time,ra,dec)
    back = li.DetFrameToEquatorial(LHO,LLO,*forward)
    delta=np.array([abs(time-back[0]),abs(ra-back[1]),abs(dec-back[2])])
    if any(delta > tolerance):
        return False
    else:
        return True


ras=np.random.uniform(low=0,high=lal.TWOPI,size=Ntest)
decs=np.random.uniform(low=-lal.PI/2.0,high=lal.PI/2.0,size=Ntest)
times=np.random.uniform(low=0,high=lal.DAYSID_SI,size=Ntest)

# Test mapping back and forward
if all(map(lambda x: test_invertable(*x), zip(ras,decs,times) ) ):
    sys.exit(0)
else:
    sys.exit(1)

