##python
# Copyright (C) 2025 Karl Wette
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
# Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA  02110-1301  USA

## \file
## \ingroup lalpulsar_bin_SFTTools
"""Output the public directory for the given SFT name(s), following the
convention detailed in the SFT spec (T040164)."""


import argparse
import sys
import os

from lalpulsar.public_sft_directory import public_sft_directory
from lalpulsar import git_version

__author__ = "Karl Wette <karl.wette@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date


def parse_command_line():
    # parse command line
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-F",
        "--format",
        dest="fmt",
        type=str,
        default="{dir}",
        help="format output using this Python {}-style formatter, where: {dir}=public SFT directory, {name}=SFT name, {orig}=original SFT path",
    )
    parser.add_argument(
        "SFT_paths",
        type=str,
        nargs="*",
        help="SFT paths; if not given, read from standard input",
    )
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_command_line()

    # iterate either over command lines or standard input
    SFT_path_source = args.SFT_paths
    if not SFT_path_source:
        SFT_path_source = sys.stdin

    for SFT_path_str in SFT_path_source:

        # original SFT path
        fmt_vals = {"orig": SFT_path_str.rstrip()}

        # get original SFT path and SFT name
        _, fmt_vals["name"] = os.path.split(fmt_vals["orig"])

        # get public SFT directory
        fmt_vals["dir"] = public_sft_directory(fmt_vals["name"])

        # print output
        print(args.fmt.format(**fmt_vals))
