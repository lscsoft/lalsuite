# Copyright (C) 2019 Duncan Macleod
# Copyright (C) 2025 Leo Singer
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

"""Tests for lal.series"""

import io
import sys

import igwn_ligolw.ligolw
import igwn_ligolw.utils
import numpy as np
import pytest

import lal.series


@pytest.mark.parametrize(
    "scalar_type", ["REAL4", "REAL8", "COMPLEX8", "COMPLEX16"]
)
@pytest.mark.parametrize("domain", ["Frequency", "Time"])
@pytest.mark.parametrize(
    "encoding,assert_array_equal",
    [
        ["Text", np.testing.assert_array_almost_equal],
        ["base64", np.testing.assert_array_equal]
    ]
)
def test_build_series(scalar_type, domain, encoding, assert_array_equal):
    class_name = f"{scalar_type}{domain}Series"
    build_series = getattr(lal.series, f"build_{class_name}")
    parse_series = getattr(lal.series, f"parse_{class_name}")
    create_series = getattr(lal, f"Create{class_name}")

    series = create_series(
        "name", lal.LIGOTimeGPS(3141, 59), 42.0, 2.5, lal.KiloGramUnit, 1234)
    series.data.data = np.random.uniform(size=series.data.data.shape)

    fileobj = io.BytesIO()
    xmldoc = igwn_ligolw.ligolw.Document()
    xmldoc.appendChild(
        build_series(series, comment="This is a comment", encoding=encoding))
    igwn_ligolw.utils.write_fileobj(xmldoc, fileobj)
    print(fileobj.getvalue())
    xml_str = fileobj.getvalue().decode()
    if encoding == igwn_ligolw.ligolw.Encoding.enc(igwn_ligolw.ligolw.Stream.Encoding.default):
        assert 'Encoding=' not in xml_str
    else:
        assert f'Encoding="{encoding}' in xml_str
    assert "This is a comment" in xml_str

    fileobj.seek(0)
    xmldoc = igwn_ligolw.utils.load_fileobj(fileobj)
    new_series = parse_series(xmldoc)
    assert new_series.name == series.name
    assert new_series.f0 == series.f0
    assert_array_equal(new_series.data.data, series.data.data)


if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-test-series.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
