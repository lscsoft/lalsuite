# Copyright (C) 2020 Rodrigo Tenorio, David Keitel
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
import sys
import numpy as np
from numpy.testing import assert_allclose
import lal
import lalpulsar
import pytest

# load Earth and Sun ephemerides
earth_ephem_file = os.path.join(
    os.environ["LAL_TEST_PKGDATADIR"], "earth00-19-DE405.dat.gz"
)
sun_ephem_file = os.path.join(
    os.environ["LAL_TEST_PKGDATADIR"], "sun00-19-DE405.dat.gz"
)
ephemerides = lalpulsar.InitBarycenter(earth_ephem_file, sun_ephem_file)

# parameters to make fake data
IFOs = ["H1","L1"]
numDet = len(IFOs)
sqrtSX = ["1e-23"]
f_min = 20
f_max = 21
T0 = 900000000
Tspan = 86400
Tsft = 1800
Toverlap = 0
dataParams = lalpulsar.CWMFDataParams()
dataParams.fMin = f_min
dataParams.Band = f_max - f_min
lalpulsar.ParseMultiLALDetector(dataParams.multiIFO, IFOs)
lalpulsar.ParseMultiNoiseFloor(dataParams.multiNoiseFloor, sqrtSX, numDet)
multiTS = lalpulsar.MakeMultiTimestamps(T0, Tspan, Tsft, Toverlap, numDet)
dataParams.multiTimestamps = multiTS
dataParams.randSeed = 42
dataParams.inputMultiTS = None
injectionSources = lalpulsar.PulsarParamsVector() # empty struct, don't actually inject a signal

# make fake data directly in memory
multiSFTs = 0
print("Generating SFTs with CWMakeFakeMultiData...")
multiSFTs, _ = lalpulsar.CWMakeFakeMultiData(
    multiSFTs,
    None,
    injectionSources=injectionSources,
    dataParams=dataParams,
    edat=ephemerides,
)
assert multiSFTs.length == numDet
numSFTs = lal.CreateUINT4Vector(numDet)
numSFTs.data = [sft.length for sft in multiSFTs.data]
f0 = multiSFTs.data[0].data[0].f0
deltaF = multiSFTs.data[0].data[0].deltaF
numBins = multiSFTs.data[0].data[0].data.length
print("numIFOs:", numDet)
print("numSFTs per IFO:", numSFTs.data)
print("numBins per SFT:", numBins)

# make an additional hard copy of the SFTs
# to test PSD computation with/without normalizing in place
print("Copying SFTs...")
multiSFTs_copy = lalpulsar.CreateEmptyMultiSFTVector(numSFTs)
for X,sftvec in enumerate(multiSFTs.data):
    for k,sft in enumerate(sftvec.data):
        lalpulsar.CopySFT(multiSFTs_copy.data[X].data[k],sft)

def test_ComputePSDandNormSFTPower():

    rtol = 1e-15 # float comparison tolerance

    print("Calling ComputePSDandNormSFTPower() with normalizeSFTsInPlace=False...")
    psd, multiPSDVector, normSFT = lalpulsar.ComputePSDandNormSFTPower(
        inputSFTs=multiSFTs,
        returnMultiPSDVector=True,
        returnNormSFT=True,
        blocksRngMed=101,
        PSDmthopSFTs=1,  # arithmean
        PSDmthopIFOs=0,  # arithsum
        nSFTmthopIFOs=1,  # arithmean
        nSFTmthopSFTs=0,  # arithsum
        normalizeByTotalNumSFTs=False,
        FreqMin=f_min,
        FreqBand=f_max - f_min,
        normalizeSFTsInPlace=False,
    )
    print("PSD from ComputePSDandNormSFTPower() with normalizeSFTsInPlace=False:", psd.data)
    assert psd.length == numBins
    print("normSFT:", normSFT.data)
    assert normSFT.length == numBins
    assert multiPSDVector.length == numDet
    summed_average_IFO_PSD_vector = np.zeros_like(psd.data)
    for X in range(numDet):
        assert multiPSDVector.data[X].length == numSFTs.data[X]
        assert multiPSDVector.data[X].data[0].data.length == psd.length
        assert multiPSDVector.data[X].data[0].f0 == f0
        assert multiPSDVector.data[X].data[0].deltaF == deltaF
        average_IFO_PSD_vector = np.array(
            [
                np.mean(
                    [multiPSDVector.data[X].data[n].data.data[k] for n in range(numSFTs.data[X])]
                )
                for k in range(numBins)
            ]
        )
        print("multiPSDVector[{:d}] averaged: {}".format(X, average_IFO_PSD_vector))
        summed_average_IFO_PSD_vector += average_IFO_PSD_vector
    print("sum of per-IFO averaged multiPSDVector entries:", summed_average_IFO_PSD_vector)
    assert_allclose(
        summed_average_IFO_PSD_vector,
        psd.data,
        rtol=rtol,
    )
    print("Relative agreement to within {:e}.".format(rtol))

    print("Calling ComputePSDandNormSFTPower() again, with normalizeSFTsInPlace=True...")
    psd2, multiPSDVector, normSFT2 = lalpulsar.ComputePSDandNormSFTPower(
        inputSFTs=multiSFTs,
        returnMultiPSDVector=True,
        returnNormSFT=True,
        blocksRngMed=101,
        PSDmthopSFTs=1,  # arithmean
        PSDmthopIFOs=0,  # arithsum
        nSFTmthopIFOs=1,  # arithmean
        nSFTmthopSFTs=0,  # arithsum
        normalizeByTotalNumSFTs=False,
        FreqMin=f_min,
        FreqBand=f_max - f_min,
        normalizeSFTsInPlace=True,
    )
    print("PSD from ComputePSDandNormSFTPower() with normalizeSFTsInPlace=True:", psd2.data)
    assert psd2.length == numBins
    assert_allclose(
        psd2.data,
        psd.data,
        rtol=rtol,
    )
    print("Relative PSD agreement to within {:e}.".format(rtol))
    print("normSFT:", normSFT2.data)
    assert normSFT2.length == numBins
    assert_allclose(
        normSFT2.data,
        normSFT.data,
        rtol=rtol,
    )
    print("Relative normSFT agreement to within {:e}.".format(rtol))

    print("Calling ComputePSDfromSFTs() on original (now normalized) SFTs...")
    psd3 = lalpulsar.ComputePSDfromSFTs(
        inputSFTs=multiSFTs,
        blocksRngMed=101,
        PSDmthopSFTs=1,  # arithmean
        PSDmthopIFOs=0,  # arithsum
        normalizeByTotalNumSFTs=False,
        FreqMin=f_min,
        FreqBand=f_max - f_min,
    )
    print("PSD from ComputePSDfromSFTs() on normalized SFTs:", psd3.data)
    assert psd3.length == numBins
    # these will be different, so not comparing the output

    print("Calling ComputePSDfromSFTs() on copied SFTs...")
    psd4 = lalpulsar.ComputePSDfromSFTs(
        inputSFTs=multiSFTs_copy,
        blocksRngMed=101,
        PSDmthopSFTs=1,  # arithmean
        PSDmthopIFOs=0,  # arithsum
        normalizeByTotalNumSFTs=False,
        FreqMin=f_min,
        FreqBand=f_max - f_min,
    )
    print("PSD from ComputePSDfromSFTs() on copied SFTs:", psd4.data)
    assert psd4.length == numBins
    assert_allclose(
        psd4.data,
        psd.data,
        rtol=rtol,
    )
    print("Relative agreement to within {:e}.".format(rtol))

if __name__ == "__main__":
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-computePSD.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
