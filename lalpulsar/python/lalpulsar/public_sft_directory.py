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

## \defgroup lalpulsar_py_public_sft_directory Public SFT Directory Function
## \ingroup lalpulsar_python
"""
Implements the public SFT directory naming convention detailed in the SFT spec (T040164)
"""

import os

__author__ = "Karl Wette <karl.wette@ligo.org>"
__version__ = "Consistent with SFT spec LIGO-T040164-v2"
__date__ = "2022"


def public_sft_directory(filename):
    # split off the filename, if any
    _, filename = os.path.split(filename)

    # break into S, D, G, T tokens
    lvl1 = filename.split("-")
    if len(lvl1) == 4:
        S, D, G, T_extn = lvl1
    else:
        raise ValueError(f'"{filename}" does not contain 4 tokens separated by "-"')

    # break D into 4 or 5 tokens
    lvl2 = D.split("_")
    if len(lvl2) == 5:
        D1, D2, D3, D4, D5 = lvl2
        if not D5.startswith("NBF"):
            raise ValueError(
                f'"{filename}" token "{D5}" does not match the narrow-band SFT filename spec'
            )
    elif len(lvl2) == 4:
        D1, D2, D3, D4 = lvl2
        D5 = None
    else:
        raise ValueError(
            f'"{filename}" token "{D}" does not contain 4 or 5 tokens separated by "_"'
        )
    if not "+" in D4:
        raise ValueError(
            f'"{filename}" token "{D4}" does not match the public SFT filename spec'
        )

    # build common SFT directory part
    SFT_base_directory = f"{D2}_{D3}_{D4}"

    if D5:
        # build narrow-band SFT directory
        SFT_base_directory += "_NARROWBAND"
        SFT_subdirectory = D5

    else:
        # build broad-band SFT directory
        SFT_base_directory += "_BROADBAND"
        G_million = int(G) // 1000000
        SFT_subdirectory = f"GPS{G_million}M"

    return os.path.join(SFT_base_directory, SFT_subdirectory)


def public_sft_directory_readme_md():
    # return README.md describing the SFT file/directory naming scheme
    path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(path, "public_sft_directory_README.md")) as f:
        return f.read()
