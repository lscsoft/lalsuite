##python
# Copyright (C) 2024 Evan Goetz
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
"""Move SFTs between directories."""

import argparse
import os
import sys
import shutil
from contextlib import contextmanager
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

from lal import LALERRORBIT, LALWARNINGBIT, LALINFOBIT, LALTRACEBIT
from lal import GetDebugLevel, ClobberDebugLevel

from lalpulsar import git_version
from lalpulsar import ValidateSFTFile, SFTErrorMessage

__author__ = "Evan Goetz <evan.goetz@ligo.org>"
__credits__ = "Karl Wette <karl.wette@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date


@contextmanager
def silence_xlal_error_messages():
    saveDebugLevel = GetDebugLevel()
    silentDebugLevel = saveDebugLevel & ~(
        LALERRORBIT | LALWARNINGBIT | LALINFOBIT | LALTRACEBIT
    )
    ClobberDebugLevel(silentDebugLevel)
    try:
        yield None
    finally:
        ClobberDebugLevel(saveDebugLevel)


def parse_command_line():
    # parse command line
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-p", "--processes", type=int, default=1, help="number of moving processes"
    )
    parser.add_argument(
        "-n",
        "--no-validate",
        dest="validate",
        action="store_false",
        help="do not validate destination SFTs",
    )
    parser.add_argument(
        "-s", "--source-directory", type=str, help="SFT source directory"
    )
    parser.add_argument(
        "-c",
        "--channels",
        type=str,
        nargs="+",
        help="Channel names (must be the same number as " "--dest-directory arguments)",
    )
    parser.add_argument(
        "-d",
        "--dest-directory",
        type=str,
        nargs="+",
        help="SFT destination directory (must be the same "
        "number as --channels arguments)",
    )
    args = parser.parse_args()

    # check arguments
    if args.processes <= 0:
        parser.error("--processes must be strictly positive")
    if not os.path.isdir(args.source_directory):
        parser.error("source_directory is not a directory")
    if len(args.channels) != len(args.dest_directory):
        parser.error(
            "Number of channel arguments must equal number of " "directory arguments"
        )

    return args


def find_SFT_files(source_directory, channel, dest_directory):
    # find SFT files for a specific channel and return the source-destiation tuple
    src_dest_paths = []

    # channel format in filename
    chan = channel.split(":")[1].replace("-", "").replace("_", "")

    # find source SFT files
    for src_root, _, src_files in os.walk(source_directory):
        for src_file in src_files:
            if src_file.endswith(".sft") and chan in src_file:
                src_path = os.path.join(src_root, src_file)
                _, src_name = os.path.split(src_path)

                # build SFT destination path
                dest_path = os.path.join(dest_directory, src_name)

                # add to outputs
                src_dest_paths.append((src_path, dest_path))

    return src_dest_paths


def make_dest_dirs(dest_dirs):
    # make destination SFT directories
    print(f"{__file__}: making {len(dest_dirs)} directories ...", flush=True)
    for dest_dir in dest_dirs:
        os.makedirs(dest_dir, exist_ok=True)
    print(f"{__file__}: making {len(dest_dirs)} directories ... done\n", flush=True)


def move_SFT_file(src_path, dest_path, validate):
    # move SFT with a temporary extension
    tmp_dest_path = dest_path + "_TO_BE_VALIDATED"
    shutil.move(src_path, tmp_dest_path)

    # validate SFT if requested
    if validate:
        with silence_xlal_error_messages() as _:
            validate_errorcode = ValidateSFTFile(tmp_dest_path)
        if validate_errorcode != 0:
            validate_errorstr = SFTErrorMessage(validate_errorcode)
            return (tmp_dest_path, validate_errorstr)

    # move destination SFT to final location
    os.rename(tmp_dest_path, dest_path)

    return None


def move_all_SFT_files(src_dest_paths, validate, processes):
    validate_errors = []

    # create executor
    print(f"{__file__}: copying {len(src_dest_paths)} SFTs ...", flush=True)
    with ProcessPoolExecutor(max_workers=args.processes) as executor:
        # submit tasks
        pool = [
            executor.submit(move_SFT_file, src_path, dest_path, validate)
            for src_path, dest_path in src_dest_paths
        ]

        # collect tasks
        for task in tqdm(as_completed(pool), total=len(pool)):
            validate_error = task.result()
            if validate_error is not None:
                validate_errors.append(validate_error)

    print("")

    # show any validation errors
    if validate_errors:
        print(
            f"{__file__}: failed to validate {len(validate_errors)} SFTs after copying:",
            flush=True,
        )
        for tmp_dest_path, validate_errorstr in validate_errors:
            print(f"  {tmp_dest_path}\n    {validate_errorstr}", flush=True)
        sys.exit(1)

    print(f"{__file__}: copying {len(src_dest_paths)} SFTs ... done\n", flush=True)


if __name__ == "__main__":
    args = parse_command_line()

    dest_dirs = set()
    src_dest_paths = []
    for idx, c in enumerate(args.channels):
        src_dest = find_SFT_files(args.source_directory, c, args.dest_directory[idx])

        dest_dirs.add(args.dest_directory[idx])
        src_dest_paths.extend(src_dest)

    make_dest_dirs(dest_dirs)

    move_all_SFT_files(src_dest_paths, args.validate, args.processes)

    print(f"{__file__}: DONE", flush=True)
