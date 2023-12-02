import sys
import subprocess as sp

import numpy as np

import pytest

# from LALConstants.h
LAL_C_SI = 299792458e0
LAL_REARTH_SI = 6378136.6
LAL_AU_SI = 149597870700e0
MAX_DT_DETS_SSB = (LAL_AU_SI + LAL_REARTH_SI) / LAL_C_SI

# compute times for a variety of telescopes and sky positions
GPS_SSB = 1100000000.0
TELESCOPES = (
    "H1",
    "L1",
    "V1",
    "G1",
    "T1",
    "GBT",
    "PKS",
    "JBO",
    "AO",
    "EFF",
    "NRT",
    "HOB",
    "HART",
    "VLA",
    "WSRT",
)


@pytest.mark.parametrize("telescope", TELESCOPES)
@pytest.mark.parametrize("ra", np.arange(0, 24, 8))
def test_ssb(telescope, ra):
    cmd = [
        "lalpulsar_ssbtodetector",
        "--gps",
        str(GPS_SSB),
        "--ra",
        f"{ra}:21:34.76",
        "--dec",
        "-1:53:12.36",
        "-t",
        str(telescope),
    ]
    out = sp.check_output(cmd)
    gps_det = float(out)
    assert abs(gps_det - GPS_SSB) < MAX_DT_DETS_SSB


if __name__ == "__main__":
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-ssbtodetector.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
