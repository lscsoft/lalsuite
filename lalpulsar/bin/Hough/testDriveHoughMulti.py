# Copyright (C) 2022 Rodrigo Tenorio, David Keitel
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

import subprocess

# run this test with "make check TESTS=testDriveHoughMulti"
# everything will be automatically put into testDriveHoughMulti.testdir
testdir = "."

mfd_cml = " ".join(
    [
        "lalpulsar_Makefakedata_v5",
        "--outSingleSFT=TRUE",
        f"--outSFTdir={testdir}",
        '--outLabel="simulated_signal"',
        '--IFOs="H1","L1"',
        '--sqrtSX="1e-22"',
        "--startTime=1000000000",
        f"--duration={10 * 1800}",
        "--fmin=100",
        "--Band=1",
        "--Tsft=1800",
    ]
)

dhm_cml = " ".join(
    [
        "lalpulsar_DriveHoughMulti",
        "--sunEphemeris",
        "sun00-40-DE405.dat.gz",
        "--earthEphemeris",
        "earth00-40-DE405.dat.gz",
        "--f0",
        "100.5",
        "--freqBand",
        "0.001",
        "--skyRegion",
        '"allsky"',
        "--dAlpha",
        "0.5",
        "--dDelta",
        "0.5",
        "--peakThreshold",
        "1.6",
        "--weighAM",
        "--weighNoise",
        "--sftData",
        f'"{testdir}/*-10*.sft"',
        "--numCand",
        "100",
    ]
)

print(f"Running MFD to generate SFTs: {mfd_cml}")
subprocess.check_call(mfd_cml, shell=True)
print(f"Running test: {dhm_cml}")
subprocess.check_call(dhm_cml, shell=True)
