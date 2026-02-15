# Copyright (C) 2022 Karl Wette
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

from lalpulsar.public_sft_directory import public_sft_directory


SFT_FILENAME_TO_DIRECTORY = (
    (
        "H-1_H1_1800SFT_O4RUN+R1+Choft+WTKEY5-1257800000-1800.sft",
        "H1_1800SFT_O4RUN+R1+Choft+WTKEY5_BROADBAND-1257",
    ),
    (
        "H-1_H1_1800SFT_O4RUN+R1+Choft+WTKEY5-1257901800-1800.sft",
        "H1_1800SFT_O4RUN+R1+Choft+WTKEY5_BROADBAND-1257",
    ),
    (
        "H-1_H1_1800SFT_O4RUN+R1+Choft+WTKEY5-1258003600-1800.sft",
        "H1_1800SFT_O4RUN+R1+Choft+WTKEY5_BROADBAND-1258",
    ),
    (
        "H-1_H1_1800SFT_O4RUN+R1+Choft+WTKEY5-1258105400-1800.sft",
        "H1_1800SFT_O4RUN+R1+Choft+WTKEY5_BROADBAND-1258",
    ),
    (
        "H-1_H1_1800SFT_O4RUN+R1+Choft+WTKEY5-1258207200-1800.sft",
        "H1_1800SFT_O4RUN+R1+Choft+WTKEY5_BROADBAND-1258",
    ),
    (
        "H-5_H1_1800SFT_O5SIM+R2+Choft+WHANN_NBF0010Hz0W0008Hz0-1257800000-9000.sft",
        "H1_1800SFT_O5SIM+R2+Choft+WHANN_NARROWBAND-NBF0010Hz0W0008Hz0",
    ),
    (
        "H-5_H1_1800SFT_O5SIM+R2+Choft+WHANN_NBF0018Hz0W0008Hz0-1257900000-9000.sft",
        "H1_1800SFT_O5SIM+R2+Choft+WHANN_NARROWBAND-NBF0018Hz0W0008Hz0",
    ),
    (
        "H-5_H1_1800SFT_O5SIM+R2+Choft+WHANN_NBF0026Hz0W0008Hz0-1258000000-9000.sft",
        "H1_1800SFT_O5SIM+R2+Choft+WHANN_NARROWBAND-NBF0026Hz0W0008Hz0",
    ),
    (
        "H-5_H1_1800SFT_O5SIM+R2+Choft+WHANN_NBF0034Hz0W0008Hz0-1258100000-9000.sft",
        "H1_1800SFT_O5SIM+R2+Choft+WHANN_NARROWBAND-NBF0034Hz0W0008Hz0",
    ),
    (
        "H-5_H1_1800SFT_O5SIM+R2+Choft+WHANN_NBF0042Hz0W0008Hz0-1258200000-9000.sft",
        "H1_1800SFT_O5SIM+R2+Choft+WHANN_NARROWBAND-NBF0042Hz0W0008Hz0",
    ),
)


@pytest.mark.parametrize("filename,directory", SFT_FILENAME_TO_DIRECTORY)
def test_public_sft_directory(filename, directory):
    assert public_sft_directory(filename) == directory


if __name__ == "__main__":
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-public_sft_directory.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
