# Copyright (C) 2022 Karl Wette
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

"""Copy SFTs between directories. The destination directory is organised
following the convention detailed in the SFT spec (T040164)."""

import argparse
import os
import sys
import time
import shutil
from contextlib import contextmanager
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

from lal import LALERRORBIT, LALWARNINGBIT, LALINFOBIT, LALTRACEBIT
from lal import GetDebugLevel, ClobberDebugLevel

from lalpulsar import git_version
from lalpulsar import ValidateSFTFile, SFTErrorMessage
from lalpulsar.public_sft_directory import public_sft_directory
from lalpulsar.public_sft_directory import public_sft_directory_readme_md

__author__ = "Karl Wette <karl.wette@ligo.org>"
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
        "-p", "--processes", type=int, default=1, help="number of copying processes"
    )
    parser.add_argument(
        "-n",
        "--no-validate",
        dest="validate",
        action="store_false",
        help="do not validate destination SFTs",
    )
    parser.add_argument(
        "--readme-md",
        dest="readme_md",
        action="store_true",
        help="write README.md in the destination directory",
    )
    parser.add_argument("source_directory", type=str, help="SFT source directory")
    parser.add_argument("dest_directory", type=str, help="SFT destination directory")
    args = parser.parse_args()

    # check arguments
    if args.processes <= 0:
        parser.error("--processes must be strictly positive")
    if not os.path.isdir(args.source_directory):
        parser.error("source_directory is not a directory")
    if not os.path.isdir(args.dest_directory):
        parser.error("dest_directory is not a directory")

    return args


def find_SFT_files(source_directory, dest_directory):
    dest_dirs = set()
    src_dest_paths = []

    # find source SFT files
    t0 = time.time()
    num_SFTs = 0
    print_progress = 100
    print_progress_step = 100
    print_progress_max = 1000
    for src_root, _, src_files in os.walk(source_directory):
        for src_file in src_files:
            if src_file.endswith(".sft"):
                src_path = os.path.join(src_root, src_file)
                _, src_name = os.path.split(src_path)

                # build SFT destination directory
                dest_dir = os.path.join(dest_directory, public_sft_directory(src_name))
                dest_path = os.path.join(dest_dir, src_name)

                # add to outputs
                dest_dirs.add(dest_dir)
                src_dest_paths.append((src_path, dest_path))

                # print progress
                num_SFTs += 1
                if num_SFTs % print_progress == 0:
                    dt = time.time() - t0
                    print(
                        f"{__file__}: found {num_SFTs} SFTs in {dt:0.1f} seconds",
                        flush=True,
                    )
                    print_progress += print_progress_step
                    if print_progress == print_progress_max:
                        print_progress_step *= 10
                        print_progress_max *= 10

    print(f"{__file__}: found {num_SFTs} SFTs\n", flush=True)

    return dest_dirs, src_dest_paths


def make_dest_dirs(dest_dirs):
    # make destination SFT directories
    print(f"{__file__}: making {len(dest_dirs)} directories ...", flush=True)
    for dest_dir in dest_dirs:
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
    print(f"{__file__}: making {len(dest_dirs)} directories ... done\n", flush=True)


def copy_SFT_file(src_path, dest_path, validate):
    # copy SFT with a temporary extension
    tmp_dest_path = dest_path + "_TO_BE_VALIDATED"
    shutil.copyfile(src_path, tmp_dest_path)

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


def copy_all_SFT_files(src_dest_paths, validate, processes):
    validate_errors = []

    # create executor
    print(f"{__file__}: copying {len(src_dest_paths)} SFTs ...", flush=True)
    with ProcessPoolExecutor(max_workers=args.processes) as executor:
        # submit tasks
        pool = [
            executor.submit(copy_SFT_file, src_path, dest_path, validate)
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


def write_readme_md(dest_directory):
    # write README.md
    with open(os.path.join(dest_directory, "README.md"), "w") as f:
        f.write(public_sft_directory_readme_md())


if __name__ == "__main__":
    args = parse_command_line()

    dest_dirs, src_dest_paths = find_SFT_files(
        args.source_directory, args.dest_directory
    )

    make_dest_dirs(dest_dirs)

    copy_all_SFT_files(src_dest_paths, args.validate, args.processes)

    if args.readme_md:
        write_readme_md(args.dest_directory)

    print(f"{__file__}: DONE", flush=True)
