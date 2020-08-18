# Copyright (C) 2020 Pep Covas, David Keitel, Rodrigo Tenorio
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

import sys
import pytest
import lalpulsar
import numpy as np
import os
try:
    from pathlib import Path
except ImportError as exc:  # probably macports
    import warnings
    warnings.warn(str(exc))
    sys.exit(77)


def test_ReadSFDB():

    # reference parameters of the test data file
    testfilename = "V1_mfdv5_20191114_043831.SFDB09"
    Tsft = 1024
    startTime = 1257741529
    SFToverlap = 512
    IFOs = ["V1"]
    fmin = 49
    fmax = 51
    signal_freq = 50

    # find reference SFDB file in a portable way
    TEST_PATH = Path(os.getenv(
        "LAL_TEST_SRCDIR",
        Path(__file__).parent,
    )).absolute()
    sfdb_file = str(TEST_PATH / testfilename)

    # read SFTs from SFDB
    multi_sfts_from_sfdb = lalpulsar.ReadSFDB(
        f_min=fmin,
        f_max=fmax,
        file_pattern=sfdb_file,
        timeStampsStarting=None,
        timeStampsFinishing=None,
    )
    print("Got SFTs for {:d} detectors from SFDBs.".format(multi_sfts_from_sfdb.length))
    assert multi_sfts_from_sfdb.length == len(IFOs)
    sfts_from_sfdb = multi_sfts_from_sfdb.data[0]
    nSFTs = sfts_from_sfdb.length
    print("Got {:d} SFTs for {:s} from SFDBs.".format(nSFTs, sfts_from_sfdb.data[0].name))
    assert nSFTs == 3
    for k,sft in enumerate(sfts_from_sfdb.data):
        print("SFT {:d}/{:d}: epoch={:d}".format(k,nSFTs,sft.epoch.gpsSeconds))
        f0 = sft.f0
        deltaF = sft.deltaF
        nBins = sft.data.length
        absval = np.abs(sft.data.data)
        maxabs = np.max(absval)
        maxbin = np.argmax(absval)
        maxfreq = f0 + deltaF * maxbin

        # simple sanity checks of returned properties
        print("f0 SFDB: {}".format(f0))  # Starting frequency
        print("df SFDB: {}".format(deltaF))  # Frequency spacing between bins
        print("nBins SFDB: {}".format(nBins))  # Number of frequency bins in one SFT
        print("max abs value {:.4e} at freq[{:d}]={:.6f} (expected signal frequency: {:.6f})".format(maxabs,maxbin,maxfreq,signal_freq))
        if k==0:
            maxbin0 = maxbin
            maxabs0 = maxabs
        assert sft.name == IFOs[0]
        assert f0 == fmin
        assert abs(deltaF-1./Tsft) < 1e-16
        assert nBins == 2*Tsft+1
        assert maxbin == maxbin0
        assert abs(maxfreq-signal_freq) < 4*deltaF


if __name__ == "__main__":
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-ReadSFDB.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
