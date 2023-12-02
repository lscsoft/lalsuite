# Copyright (C) 2021 Rodrigo Tenorio
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

"""
Test code for lineFileParser.py
"""

import sys
import tempfile

try:
    from pathlib import Path
except ImportError as exc:  # probably macports
    import warnings

    warnings.warn(str(exc))
    sys.exit(77)

import numpy as np
from numpy.testing import assert_allclose

from lalpulsar.lineFileParser import LineFileParser

import pytest


@pytest.fixture
def identified_lines():
    return {
        "left_side": np.array(
            [
                15.95,
                0.95,
                1.95,
                2.95,
                3.95,
                4.95,
                9.95,
                19.9,
                29.85,
                39.8,
                49.75,
                109.0,
                209.0,
                309.0,
            ]
        ),
        "right_side": np.array(
            [
                16.05,
                1.05,
                2.05,
                3.05,
                4.05,
                5.05,
                10.05,
                20.1,
                30.15,
                40.2,
                50.25,
                111.0,
                211.0,
                311.0,
            ]
        ),
    }


@pytest.fixture
def unidentified_lines():
    return {
        "left_side": np.array([10.0, 10.2, 10.4]),
        "right_side": np.array([10.2, 10.4, 10.6]),
    }


@pytest.yield_fixture
def identified_lines_file():
    filename = "IdentifiedLines.csv"
    line_data = (
        "# git repo: <url-to-repo>\n# git tag:  <tag>\n# git hash: <hash>\n"
        "Frequency or frequency spacing [Hz],Type (0:line; 1:comb; 2:comb with scaling width),"
        "Frequency offset [Hz],First visible harmonic, Last visible harmonic, "
        "Left width [Hz], Right width [Hz], Comments\n"
        "16.0,0,0,1,1,0.05,0.05,Calibration line\n"
        "1.0,1,0,1,5,0.05,0.05,Comb\n"
        "10.0,2,0,1,5,0.05,0.05,Scaling Comb\n"
        "100.0,1,10,1,3,1.0,1.0,Comb with offset\n"
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = Path(tmpdir) / filename
        filepath.write_text(line_data)
        yield filepath


@pytest.yield_fixture
def unidentified_lines_file():
    filename = "UnidentifiedLines.csv"
    line_data = (
        "# git repo: <url-to-repo>\n# git tag:  <tag>\n# git hash: <hash>\n"
        "Frequency [Hz], Comments\n"
        "10.1,Unidentified line\n"
        "10.3,Unidentified line\n"
        "10.5,Unidentified line\n"
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = Path(tmpdir) / filename
        filepath.write_text(line_data)
        yield filepath


@pytest.fixture
def identified_lines_header():
    return {
        "frequency": 0,
        "left_wing": 5,
        "right_wing": 6,
        "first_harmonic": 3,
        "last_harmonic": 4,
        "comb_type": 1,
        "offset": 2,
    }


@pytest.fixture
def unidentified_lines_header():
    return {"frequency": 0}


def test_identified_lines_parsing(
    identified_lines, identified_lines_file, identified_lines_header
):
    """
    Test that identified line files are properly read and expanded
    by comparing left, right limits to the expected values.
    """

    for columns in [identified_lines_header, None]:
        line_parser = LineFileParser()
        line_parser.parse_identified_lines_csv(
            lines_file=identified_lines_file,
            columns=columns,
            genfromtxt_kwargs={"delimiter": ",", "skip_header": 4},
        )

        assert_allclose(identified_lines["left_side"], line_parser.lines_left_side)
        assert_allclose(identified_lines["right_side"], line_parser.lines_right_side)


def test_unidentified_lines_parsing(
    unidentified_lines, unidentified_lines_file, unidentified_lines_header
):
    """
    Test that unidentified line files are properly read and expanded
    by comparing left, right limits to the expected values.
    """

    for columns in [unidentified_lines_header, None]:
        line_parser = LineFileParser()
        line_parser.parse_unidentified_lines_csv(
            lines_file=unidentified_lines_file,
            columns=columns,
            extra_wing_Hz=0.1,
            genfromtxt_kwargs={"delimiter": ",", "skip_header": 4},
        )

        assert_allclose(unidentified_lines["left_side"], line_parser.lines_left_side)
        assert_allclose(unidentified_lines["right_side"], line_parser.lines_right_side)


def test_concatenation(
    identified_lines,
    identified_lines_file,
    identified_lines_header,
    unidentified_lines,
    unidentified_lines_file,
    unidentified_lines_header,
):
    """
    Test that files are properly read and concatenated.
    """

    line_parser = LineFileParser()
    line_parser.parse_identified_lines_csv(
        lines_file=identified_lines_file,
        columns=identified_lines_header,
        genfromtxt_kwargs={"delimiter": ",", "skip_header": 4},
    )
    line_parser.parse_unidentified_lines_csv(
        lines_file=unidentified_lines_file,
        columns=unidentified_lines_header,
        extra_wing_Hz=0.1,
        genfromtxt_kwargs={"delimiter": ",", "skip_header": 4},
    )

    assert_allclose(
        np.hstack([identified_lines["left_side"], unidentified_lines["left_side"]]),
        line_parser.lines_left_side,
    )
    assert_allclose(
        np.hstack([identified_lines["right_side"], unidentified_lines["right_side"]]),
        line_parser.lines_right_side,
    )


if __name__ == "__main__":
    args = sys.argv[1:] or ["-v", "-rs"]
    sys.exit(pytest.main(args=[__file__] + args))
