##python
# Copyright (C) 2013, 2014, 2020--2024 Evan Goetz
# Copyright (C) 2011, 2021, 2022 Karl Wette
# Copyright (C) 2005, 2007 Gregory Mendell
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
"""Creates DAGs to run jobs that generates SFTs"""

import math
import argparse
import os
import re
from pathlib import Path
from urllib.parse import urlparse

from gwdatafind import find_urls
from gwdatafind.utils import filename_metadata, file_segment

from igwn_segments import segment, segmentlist

from lalpulsar import (
    git_version,
    SFTFilenameSpec,
    FillSFTFilenameSpecStrings,
    BuildSFTFilenameFromSpec,
)

__author__ = "Evan Goetz <evan.goetz@ligo.org>, Greg Mendell"
__version__ = git_version.id
__date__ = git_version.date

cache_re = re.compile(r"^([A-Z])(\s+)(\w+)(\s+)(\d+)(\s+)(\d+)(\s+)(.+gwf)")


# REVISIONS:
# 12/02/05 gam; generate datafind.sub and MakeSFTs.sub as well as dag file in
#               PWD, with log files based subLogPath and dag filename.
# 12/28/05 gam; Add option --make-gps-dirs, -D <num>, to make directory based
#               on this many GPS digits.
# 12/28/05 gam; Add option --misc-desc, -X <string> giving misc. part of the
#               SFT description field in the filename.
# 12/28/05 gam; Add options --start-freq -F and --band -B options to enter
#               these.
# 12/28/05 gam; Add in --window-type, -w options; 0 = no window, 1 = default =
#               Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 =
#               Hann window.
# 12/28/05 gam; Add option --overlap-fraction -P (for use with windows; e.g.,
#               use -P 0.5 with -w 3 Hann windows; default is 0.0)
# 12/28/05 gam; Add --sft-version, -v option to select output SFT version (1 =
#               default is version 1 SFTs; 2 = version 2 SFTs.
# 12/28/05 gam; Add --comment-field, -c option, for comment for version 2 SFTs.
# 12/28/05 gam; Remove sample rate option
# 01/09/06 gam; Add -Z option; write SFT to .*.tmp file, then move to final
#               file name.
# 01/14/07 gam; Add -u option to specify frame struct and type; add -i option
#               to specify IFO name.
# 07/24/07 gam; Add in -q option to read in list of nodes on which to output
#               SFTs, -Q option to give node path, and -R option for number of
#               jobs per node.
# 04/XX/13 eag; Add -y option to synchronize the start times of SFTs.
# 07/24/14 eag; Change default to version 2 SFTs
# 12/2020  eag; Update script to conform to modern python3 and pep8
# 10/2020  kww; Pass args directly to writeToDag(), use Python f-strings
# 10/2022  kww; Deprecate options that have been removed from MakeSFTs
# 10/2022  kww; Parse window type as a string, parameter separated by colon
# 10/2022  kww; Merge -O and -o log path options to free up -O option
# 10/2022  kww; Implement public SFT file naming convention
# 11/2022  kww; -R command line option now used for --observing-revision
#               instead of --output-jobs-per-node, which now uses -r
# 11/2022  kww; --datafind-path and --makesfts-path accept executable names
# 03/2023  eag; Allow user to pass a frame cache file --cache-file
# 04/2023  kww; Improve documentation of --window-type argument
# 05/2023  eag; Add the --gaps flag to gw_data_find
# 08/2023  eag; Allow for multiple channel names to be provided
# 09/2023  eag; Modify use of environment variables
# 01/2024  eag; Allow skipping of channels if not in frames or too low data
#               rate
# 10/2024  eag; Modify workflow for version 3 SFTs and HTCondor file transfer
#               workflow
# 12/2024  eag; Modify workflow to use lalpulsar_MoveSFTs script instead of
#               remapping files
# 03/2025  eag; Rewrite to generate cache files running lalpulsar_MakeSFTDAG
#               enabling frame files on /home or OSDF


def sft_name_from_vars(
    obs,
    gpsstart,
    Tsft,
    channel=None,
    kind=None,
    rev=None,
    window="unknown",
    par=None,
    miscstr=None,
):
    """Create SFT file name from specification"""

    spec = SFTFilenameSpec()

    FillSFTFilenameSpecStrings(
        spec=spec,
        path=None,
        extn=None,
        detector=channel[:2],
        window_type=window,
        privMisc=miscstr,
        pubObsKind=kind,
        pubChannel=channel,
    )
    spec.pubObsRun = obs or 0
    spec.pubRevision = rev or 0
    spec.window_param = par or 0
    spec.numSFTs = 1  # MakeSFTDAG will only ever generate 1 SFT per file
    spec.SFTtimebase = Tsft
    spec.gpsStart = gpsstart
    spec.SFTspan = Tsft  # MakeSFTDAG will only ever generate 1 SFT per file

    return BuildSFTFilenameFromSpec(spec)


def get_urls(args):
    """Get frame file URL list from gwdatafind or cache file"""

    if not args.cache_file:
        urls = find_urls(
            site,
            args.input_data_type,
            segList[0][0],
            segList[-1][-1],
            match=args.datafind_match,
            urltype=args.datafind_urltype,
        )
    else:
        urls = []
        with open(args.cache_file, "r") as f:
            for line in f:
                m = cache_re.match(line)
                if m:
                    framefile = m.group(9)
                    urls.append(framefile)

    # sort the urls by gps time since find_urls() may not return a sorted list
    sorted_urls = sorted(urls, key=lambda x: file_segment(x)[0])

    return sorted_urls


def make_cache(
    urls,
    job_seg,
):
    """Make a frame list and cache list from a list of URLs"""

    cache = []  # list of lines for the cache file
    frames = []  # list of frame filenames used in the job
    for idx, url in enumerate(urls):
        obs, desc, dataseg = filename_metadata(url)
        dataseg = segment(dataseg)
        if dataseg.disjoint(job_seg) < 0:
            continue
        if dataseg.disjoint(job_seg) > 0:
            break
        if dataseg.intersects(job_seg):
            framefileurl = urlparse(url)
            framefilepath = Path(framefileurl.path)

            # list in cache file if files not visible on execute node
            # otherwise use file url
            if "/home" in str(framefilepath.parent) or "osdf" in framefileurl.scheme:
                newcache = (
                    f"{obs}\t{desc}\t{dataseg[0]}\t{abs(dataseg)}\t{framefilepath.name}"
                )
            else:
                newcache = f"{obs}\t{desc}\t{dataseg[0]}\t{abs(dataseg)}\t{url}"
            cache.append(newcache)

            if "/home" in str(framefilepath.parent):
                frames.append(framefilepath)
            else:
                frames.append(url)

    return frames, cache


def writeToDag(dagFID, nodeCount, startTimeThisNode, endTimeThisNode, urls, args):
    """Write one job to DAG file"""

    MakeSFTs = f"MakeSFTs_{nodeCount}"
    tagStringOut = f"{args.tag_string}_{nodeCount}"

    job_segment = segment(
        startTimeThisNode - args.extra_datafind_time,
        endTimeThisNode + args.extra_datafind_time,
    )

    frames, cache = make_cache(urls, job_segment)

    obs, desc, dataseg = filename_metadata(urls[0])
    cacheFile = args.cache_path / f"{obs}-{job_segment[0]}-{job_segment[1]}.cache"
    with open(cacheFile, "w") as f:
        for l in cache:
            f.write(f"{l}\n")

    argList = []
    argList.append(f"-O {args.observing_run}")
    if args.observing_run > 0:
        argList.append(f"-K {args.observing_kind}")
        argList.append(f"-R {args.observing_revision}")
    elif args.misc_desc:
        argList.append(f"-X {args.misc_desc}")
    argList.append(f"-f {args.filter_knee_freq}")
    argList.append(f"-t {args.time_baseline}")
    # To work with the condor file transfer protocol, we save everything to the
    # scratch directory because the files will have unique names. Since
    # lalpulsar_MakeSFTs wants to have a path to save the file to, we provide .
    argList.append(f"-p {','.join(['.' for p in args.output_sft_path])}")
    # To work with the condor file transfer protocol, the cache file is saved
    # to the scratch directory on transfer so we just need the name, not the
    # full path
    argList.append(f"-C {cacheFile.name}")
    argList.append(f"-s {startTimeThisNode}")
    argList.append(f"-e {endTimeThisNode}")
    argList.append(f"-N {','.join(args.channel_name)}")
    argList.append(f"-F {args.start_freq}")
    argList.append(f"-B {args.band}")
    if args.comment_field:
        argList.append(f"-c {args.comment_field}")
    if ":" in args.window_type:
        window_type, window_param = args.window_type.split(":")
        window_param = float(window_param)
        argList.append(f"-w {window_type} -r {window_param}")
    else:
        window_type = args.window_type
        window_param = None
        argList.append(f"-w {window_type}")
    if args.overlap_fraction:
        argList.append(f"-P {args.overlap_fraction}")
    if args.allow_skipping:
        argList.append("--allow-skipping TRUE")
    argStr = " ".join(argList)

    # The files are going to go to specific directories, so we need to map
    # the files to their output directories
    outputfiles = []
    remap = []
    sft_start = startTimeThisNode
    sft_end = sft_start + args.time_baseline
    # loop over start times
    while sft_end <= endTimeThisNode:
        # loop over channels
        for idx, c in enumerate(args.channel_name):
            filename = sft_name_from_vars(
                args.observing_run,
                sft_start,
                args.time_baseline,
                c,
                kind=args.observing_kind,
                rev=args.observing_revision,
                window=window_type,
                par=window_param,
                miscstr=args.misc_desc,
            )
            outputfiles.append(filename)
            remap.append(f"{filename}={args.output_sft_path[idx]/filename}")

        # update start and end times
        if args.overlap_fraction:
            sft_start += int(round((1 - args.overlap_fraction) * args.time_baseline))
        else:
            sft_start += args.time_baseline
        sft_end = sft_start + args.time_baseline

    # MakeSFT job
    dagFID.write(f"JOB {MakeSFTs} {Path(dagFID.name).parent / 'MakeSFTs.sub'}\n")
    dagFID.write(f"RETRY {MakeSFTs} 1\n")
    dagFID.write(f'VARS {MakeSFTs} argList="{argStr}" cachefile="{cacheFile}" ')
    if args.transfer_frame_files:
        framefiles = ",".join([str(fr) for fr in frames])
        dagFID.write(f'framefiles="{framefiles}" ')
    dagFID.write(f'tagstring="{tagStringOut}"\n')


#
# MAIN CODE START HERE
#

parser = argparse.ArgumentParser(
    description="This script creates MakeSFTs.sub, MoveSFTs.sub, and a dag \
                 file that generates SFTs based on the options given.",
    fromfile_prefix_chars="@",
)

dag_group = parser.add_argument_group(
    "DAG organization", "Options for workflow control"
)
datafind_group = parser.add_argument_group(
    "Datafind", "Options for locating frame files"
)
makesfts_group = parser.add_argument_group(
    "SFT creation", "Options for SFT creation and output"
)
deprecated_group = parser.add_argument_group("DEPRECATED")

dag_group.add_argument(
    "-f",
    "--dag-file",
    required=True,
    type=Path,
    help="filename for .dag file (should end in .dag)",
)
dag_group.add_argument(
    "-G",
    "--tag-string",
    required=True,
    type=str,
    help="tag string used in names of various files unique to \
          jobs that will run under the DAG",
)
dag_group.add_argument(
    "-a",
    "--analysis-start-time",
    type=int,
    help="GPS start time of data from which to generate \
          SFTs (optional and unused if a segment file is given)",
)
dag_group.add_argument(
    "-b",
    "--analysis-end-time",
    type=int,
    help="GPS end time of data from which to generate SFTs \
          (optional and unused if a segment file is given)",
)
dag_group.add_argument(
    "-L",
    "--max-length-all-jobs",
    type=int,
    help="maximum total amount of data to process, in seconds \
          (optional and unused if a segment file is given)",
)
dag_group.add_argument(
    "-g",
    "--segment-file",
    type=Path,
    help="alternative file with segments to use, rather than \
          the input times",
)
dag_group.add_argument(
    "-l",
    "--min-seg-length",
    type=int,
    default=0,
    help="minimum length segments to process in seconds (used \
          only if a segment file is given)",
)
dag_group.add_argument(
    "-y",
    "--synchronize-start",
    action="store_true",
    help="synchronize the start times of the SFTs so that the \
          start times are synchronized when there are gaps in the \
          data",
)
dag_group.add_argument(
    "-o",
    "--log-path",
    type=Path,
    default="logs",
    help="path to log, output, and error files (default \
          is $PWD/logs; this directory is created if it does not \
          exist and usually should be under a local file system)",
)
dag_group.add_argument(
    "-m",
    "--max-num-per-node",
    type=int,
    default=1,
    help="maximum number of SFTs to generate on one node",
)
dag_group.add_argument(
    "-J",
    "--makesfts-path",
    type=Path,
    help="string specifying the lalpulsar_MakeSFTs executable, \
          or a path to it; if not set, will use \
          MAKESFTS_PATH env variable or system default (in that \
          order)",
)
dag_group.add_argument(
    "--movesfts-path",
    type=Path,
    help="string specifying the lalpulsar_MoveSFTs executable, \
          or a path to it; if not set, will use \
          MOVESFTS_PATH env variable or system default (in that \
          order)",
)
dag_group.add_argument(
    "-Y",
    "--request-memory",
    type=int,
    default=4096,
    help="memory allocation in MB to request from condor for \
          lalpulsar_MakeSFTs step",
)
dag_group.add_argument(
    "-s",
    "--request-disk",
    type=int,
    default=4096,
    help="disk space allocation in MB to request from condor \
          for lalpulsar_MakeSFTs step",
)
dag_group.add_argument(
    "-A",
    "--accounting-group",
    required=True,
    type=str,
    help="Condor tag for the production of SFTs",
)
dag_group.add_argument(
    "-U",
    "--accounting-group-user",
    required=True,
    type=str,
    help="albert.einstein username (do not add @LIGO.ORG)",
)
dag_group.add_argument(
    "-t",
    "--transfer-frame-files",
    action="store_true",
    help="Transfer frame files via HTCondor file transfer system. \
          This should be specified if frames are not visible to the \
          compute node file system. Ex. this should be specified if \
          frames are on /home or running on the open science grid. \
          Usually frame files are visible on CIT, LHO, LLO clusters \
          so that this does not need to be specified in that case.",
)

datafind_group.add_argument(
    "-d",
    "--input-data-type",
    required=True,
    type=str,
    help="input data type for use with the gw_data_find --type \
          option",
)
datafind_group.add_argument(
    "-x",
    "--extra-datafind-time",
    type=int,
    default=0,
    help="extra time to subtract/add from/to start/end time \
          arguments of gw_data_find",
)
datafind_group.add_argument(
    "-M",
    "--datafind-match",
    type=str,
    help="string to use with the gw_data_find --match option",
)
datafind_group.add_argument(
    "--datafind-urltype",
    type=str,
    default="file",
    choices=["file", "osdf"],
    help="String for the gw_data_find --urltype option. \
          Use 'file' if creating SFTs on a local LDG cluster. \
          Use 'osdf' if creating SFTs on the open science grid",
)
datafind_group.add_argument(
    "-e",
    "--cache-file",
    type=Path,
    help="path and filename to frame cache file to use instead \
          of gw_data_find",
)

makesfts_group.add_argument(
    "-O",
    "--observing-run",
    required=True,
    type=int,
    help="For public SFTs, observing run data the SFTs are generated from, or \
          (in the case of mock data challenge data) the observing \
          run on which the data is most closely based",
)
makesfts_group.add_argument(
    "-K",
    "--observing-kind",
    type=str,
    choices=["RUN", "AUX", "SIM", "DEV"],
    help='For public SFTs, one of: "RUN" for production SFTs of h(t) channels; \
          "AUX" for SFTs of non-h(t) channels; \
          "SIM" for mock data challenge or other simulated data; or \
          "DEV" for development/testing purposes',
)
makesfts_group.add_argument(
    "-R",
    "--observing-revision",
    type=int,
    help="For public SFTs: revision number starts at 1, and should be incremented once \
          SFTs have been widely distributed across clusters, advertised \
          as being ready for use, etc.  For example, if mistakes are found \
          in the initial SFT production run after they have been published, \
          regenerated SFTs should have a revision number of at least 2",
)
makesfts_group.add_argument(
    "-X",
    "--misc-desc",
    type=str,
    help="For private SFTs, miscellaneous part of the SFT \
          description field in the filename",
)
makesfts_group.add_argument(
    "-k",
    "--filter-knee-freq",
    required=True,
    type=float,
    help="high pass filter knee frequency used on time domain \
          data before generating SFTs",
)
makesfts_group.add_argument(
    "-T",
    "--time-baseline",
    required=True,
    type=int,
    help="time baseline of SFTs  (e.g., 60 or 1800 seconds)",
)
makesfts_group.add_argument(
    "-p",
    "--output-sft-path",
    nargs="+",
    type=Path,
    help="Path where to save the SFT files. Can specify multiple options, \
          If specifying multiple options then it is required to specify the \
          same number of output-sft-path options as the number of channels. \
          The first listed channel will have the SFTs go into the first \
          listed output-sft-path. Otherwise specify only one output path. \
          If one path is specified and more than 1 channels are specified \
          then --observing-run must be >= 1 and --observing-kind and \
          --observing-revision must be set",
)
makesfts_group.add_argument(
    "-C",
    "--cache-path",
    type=Path,
    default="cache",
    help="path to cache files that will be produced by \
          gw_data_find (default is $PWD/cache; this directory is \
          created if it does not exist and must agree with that \
          given in .sub files)",
)
makesfts_group.add_argument(
    "-N",
    "--channel-name",
    nargs="+",
    type=str,
    help="Name of input time-domain channel to read from frames. \
          Can specify multiple options. The number of channels must be \
          equal to the number of output-sft-path options given. The \
          first listed channel will have the SFTs go to the first listed \
          output-sft-path. Can only specify one channel when generating \
          private SFTs (--observing-run=0)",
)
makesfts_group.add_argument(
    "-c", "--comment-field", type=str, help="comment for SFT header"
)
makesfts_group.add_argument(
    "-F", "--start-freq", type=int, default=10, help="start frequency of the SFTs"
)
makesfts_group.add_argument(
    "-B", "--band", type=int, default=1990, help="frequency band of the SFTs"
)
makesfts_group.add_argument(
    "-w",
    "--window-type",
    type=str,
    default="tukey:0.001",
    help='type of windowing of time-domain to do \
          before generating SFTs, e.g. "rectangular", \
          "hann", "tukey:<beta in [0,1], required>"; \
          (default is "tukey:0.001", standard choice for LVK production SFTs)',
)
makesfts_group.add_argument(
    "-P",
    "--overlap-fraction",
    type=float,
    default=0,
    help="overlap fraction (for use with windows; e.g., use \
          --overlap-fraction 0.5 with --window-type hann windows)",
)
makesfts_group.add_argument(
    "--allow-skipping",
    action="store_true",
    help="allow channels to be skipped if not in frames or too low sampling \
          frequency",
)
makesfts_group.add_argument(
    "--no-validate",
    dest="validate",
    action="store_false",
    help="do not validate created SFTs",
)

##### DEPRECATED OPTIONS #####
class DeprecateAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        parser.error(
            f"Argument {self.option_strings} has been deprecated in lalpulsar_MakeSFTs"
        )


deprecated_group.add_argument(
    "-u",
    "--frame-struct-type",
    nargs=0,
    action=DeprecateAction,
    help="DEPRECATED. No longer required; \
          the frame channel type is determined automatically",
)
deprecated_group.add_argument(
    "-H",
    "--use-hoft",
    nargs=0,
    action=DeprecateAction,
    help="DEPRECATED. No longer required; \
          the frame channel type is determined automatically",
)
deprecated_group.add_argument(
    "-i",
    "--ifo",
    nargs=0,
    action=DeprecateAction,
    help="DEPRECATED. No longer required; \
          the detector prefix is deduced from the channel name",
)
deprecated_group.add_argument(
    "-D",
    "--make-gps-dirs",
    nargs=0,
    action=DeprecateAction,
    help="DEPRECATED. No longer supported",
)
deprecated_group.add_argument(
    "-Z",
    "--make-tmp-file",
    nargs=0,
    action=DeprecateAction,
    help="DEPRECATED. Default behaviour",
)
deprecated_group.add_argument(
    "-v",
    "--sft-version",
    nargs=0,
    action=DeprecateAction,
    help="DEPRECATED. No longer supported",
)
deprecated_group.add_argument(
    "-S",
    "--use-single",
    nargs=0,
    action=DeprecateAction,
    help="DEPRECATED. No longer supported",
)
deprecated_group.add_argument(
    "-q",
    "--list-of-nodes",
    type=str,
    action=DeprecateAction,
    help="DEPCRECATED. No longer supported",
)
deprecated_group.add_argument(
    "-Q",
    "--node-path",
    type=Path,
    action=DeprecateAction,
    help="DEPCRECATED. No longer supported",
)
deprecated_group.add_argument(
    "-r",
    "--output-jobs-per-node",
    type=int,
    default=0,
    action=DeprecateAction,
    help="DEPRECATED. No longer supported",
)

args = parser.parse_args()

# Some basic argument value checking
if args.observing_run < 0:
    raise parser.error("--observing-run must be >= 0")

if args.observing_run > 0 and not args.observing_kind:
    raise parser.error("--observing-run requires --observing-kind")

if args.observing_run > 0 and not args.observing_revision:
    raise parser.error("--observing-run requires --observing-revision")

if args.observing_revision and args.observing_revision <= 0:
    raise parser.error("--observing-revision must be > 0")

if args.observing_run > 0 and args.misc_desc:
    raise parser.error(
        f"--observing-run={args.observing_run} incompatible with --misc-desc"
    )

if args.misc_desc and not re.compile(r"^[A-Za-z0-9]+$").match(args.misc_desc):
    raise parser.error("--misc-desc may only contain A-Z, a-z, 0-9 characters")

if args.extra_datafind_time < 0:
    raise parser.error("--extra-datafind-time must be >= 0")

if args.filter_knee_freq < 0:
    raise parser.error("--filter-knee-freq must be >= 0")

if args.time_baseline <= 0:
    raise parser.error("--time-baseline must be > 0")

if args.overlap_fraction < 0.0 or args.overlap_fraction >= 1.0:
    raise parser.error("--overlap-fraction must be in the range [0,1)")

if args.start_freq < 0.0 or args.start_freq >= 7192.0:
    raise parser.error("--start-freq must be in the range [0,7192)")

if args.band <= 0 or args.band >= 8192.0:
    raise parser.error("--band must be in the range (0,8192)")

if args.start_freq + args.band >= 8192.0:
    raise parser.error("--start-freq + --band must be < 8192")

if args.max_num_per_node <= 0:
    raise parser.error("--max-num-per-node must be > 0")

if (
    len(args.channel_name) != len(args.output_sft_path)
    and len(args.output_sft_path) != 1
):
    raise parser.error(
        "--channel-name and --output-sft-path must be the "
        "same length or --output-sft-path must be length of 1"
    )

if len(args.channel_name) > 1 and args.observing_run == 0:
    raise parser.error(
        "When creating SFTs from multiple channels, public SFT naming "
        "convention must be used: --observing-run > 0 and set "
        "--observing-kind and --observing-revision"
    )

if args.datafind_urltype == "osdf" and not args.transfer_frame_files:
    raise parser.error(
        "--transfer-frame-files must be specified when --datafind-urltype=osdf"
    )

# Set executables for lalpulsar_MakeSFTs, and lalpulsar_MoveSFTs
makeSFTsExe = "lalpulsar_MakeSFTs"
if args.makesfts_path:
    if args.makesfts_path.is_file():
        makeSFTsExe = args.makesfts_path
    else:
        makeSFTsExe = args.makesfts_path / makeSFTsExe
elif "MAKESFTS_PATH" in os.environ:
    makeSFTsExe = Path("$ENV(MAKESFTS_PATH)") / makeSFTsExe
else:
    makeSFTsExe = Path("@LALSUITE_BINDIR@") / makeSFTsExe

moveSFTsExe = "lalpulsar_MoveSFTs"
if args.movesfts_path:
    if args.movesfts_path.is_file():
        moveSFTsExe = args.movesfts_path
    else:
        moveSFTsExe = args.movesfts_path / moveSFTsExe
elif "MOVESFTS_PATH" in os.environ:
    moveSFTsExe = Path("$ENV(MOVESFTS_PATH)") / moveSFTsExe
else:
    moveSFTsExe = Path("@LALSUITE_BINDIR@") / moveSFTsExe

# make directories to store the cache files, job logs, and SFTs
args.log_path.mkdir(exist_ok=True)
args.cache_path.mkdir(exist_ok=True)
for p in args.output_sft_path:
    p.mkdir(exist_ok=True)

# Check if segment file was given, else set up one segment from the command line
segList = segmentlist()
adjustSegExtraTime = False
if args.segment_file is not None:
    if args.min_seg_length < 0:
        raise parser.error("--min-seg-length must be >= 0")

    # the next flag causes extra time that cannot be processes to be trimmed
    # from the start and end of a segment
    adjustSegExtraTime = True

    with open(args.segment_file) as fp_segfile:
        for idx, line in enumerate(fp_segfile):
            splitLine = line.split()
            oneSeg = segment(int(splitLine[0]), int(splitLine[1]))
            if abs(oneSeg) >= args.min_seg_length:
                segList.append(oneSeg)

    if len(segList) < 1:
        raise ValueError(f"No segments found in segment file: {args.segment_file}")
else:
    if args.analysis_start_time is None:
        raise parser.error(
            "--analysis-start-time must be specified if no segment file is " "given"
        )

    if args.analysis_end_time is None:
        raise parser.error(
            "--analysis-start-time must be specified if no segment file is " "given"
        )

    if args.max_length_all_jobs is None:
        raise parser.error(
            "--max-length-all-jobs must be specified if no segment file is " "given"
        )

    # Make sure not to exceed maximum allow analysis
    if args.analysis_end_time > (args.analysis_start_time + args.max_length_all_jobs):
        args.analysis_end_time = args.analysis_start_time + args.max_length_all_jobs

    oneSeg = segment(args.analysis_start_time, args.analysis_end_time)
    segList.append(oneSeg)
# END if (args.segment_file != None)
segList.coalesce()

# Get the IFO site, which is the first letter of the channel name.
site = args.channel_name[0][0]

# Get the frame file URL list
urls = get_urls(args)

# Basic check that the frame file url list are traditionally visible on EPs
if not args.transfer_frame_files:
    for f in urls:
        if "/home" in f:
            raise parser.error(
                "--transfer-frame-files must be specified when frame files are in /home"
            )

# data segments created from the list of frame URLs
dataSegs = segmentlist()
for url in urls:
    dataSegs.append(file_segment(url))
dataSegs.coalesce()

# intersection of segList with dataSegs
segList &= dataSegs
segList.coalesce()  # just in case

# initialize count of nodes
nodeCount = 0

# Create .sub files
path_to_dag_file = args.dag_file.parent
dag_filename = args.dag_file.name
makesfts_sub = path_to_dag_file / "MakeSFTs.sub"
movesfts_sub = path_to_dag_file / "MoveSFTs.sub"

# create MakeSFTs.sub
with open(makesfts_sub, "w") as MakeSFTsFID:
    MakeSFTsLogFile = f"{args.log_path}/MakeSFTs_{dag_filename}.log"
    MakeSFTsFID.write("universe = vanilla\n")
    MakeSFTsFID.write(f"executable = {makeSFTsExe}\n")
    MakeSFTsFID.write("arguments = $(argList)\n")
    MakeSFTsFID.write(f"accounting_group = {args.accounting_group}\n")
    MakeSFTsFID.write(f"accounting_group_user = {args.accounting_group_user}\n")
    MakeSFTsFID.write(f"log = {MakeSFTsLogFile}\n")
    MakeSFTsFID.write(f"error = {args.log_path}/MakeSFTs_$(tagstring).err\n")
    MakeSFTsFID.write(f"output = {args.log_path}/MakeSFTs_$(tagstring).out\n")
    MakeSFTsFID.write("notification = never\n")
    MakeSFTsFID.write(f"request_memory = {args.request_memory}MB\n")
    MakeSFTsFID.write(f"request_disk = {args.request_disk}MB\n")
    MakeSFTsFID.write("RequestCpus = 1\n")
    MakeSFTsFID.write("should_transfer_files = yes\n")
    if args.transfer_frame_files:
        MakeSFTsFID.write("transfer_input_files = $(cachefile),$(framefiles)\n")
    else:
        MakeSFTsFID.write("transfer_input_files = $(cachefile)\n")
    if "MAKESFTS_PATH" in os.environ and not args.makesfts_path:
        MakeSFTsFID.write("getenv = MAKESFTS_PATH\n")
    if args.datafind_urltype == "osdf":
        MakeSFTsFID.write("use_oauth_services = scitokens\n")
        MakeSFTsFID.write(
            "environment = BEARER_TOKEN_FILE=$$(CondorScratchDir)/.condor_creds/scitokens.use\n"
        )
    MakeSFTsFID.write("queue 1\n")

# create MoveSFTs.sub
with open(movesfts_sub, "w") as MoveSFTsFID:
    MoveSFTsLogFile = f"{args.log_path}/MoveSFTs_{dag_filename}.log"
    MoveSFTsFID.write("universe = local\n")
    MoveSFTsFID.write(f"executable = {moveSFTsExe}\n")
    MoveSFTsFID.write("arguments = ")
    if not args.validate:
        MoveSFTsFID.write("$(opts) ")
    MoveSFTsFID.write("-s $(sourcedirectory) -c $(channels) -d $(destdirectory)\n")
    MoveSFTsFID.write(f"accounting_group = {args.accounting_group}\n")
    MoveSFTsFID.write(f"accounting_group_user = {args.accounting_group_user}\n")
    MoveSFTsFID.write(f"log = {MoveSFTsLogFile}\n")
    MoveSFTsFID.write(f"error = {args.log_path}/MoveSFTs.err\n")
    MoveSFTsFID.write(f"output = {args.log_path}/MoveSFTs.out\n")
    MoveSFTsFID.write("notification = never\n")
    MoveSFTsFID.write(f"request_memory = 1GB\n")
    MoveSFTsFID.write(f"request_disk = 10MB\n")
    MoveSFTsFID.write("RequestCpus = 1\n")
    if "MOVESFTS_PATH" in os.environ and not args.movesfts_path:
        MoveSFTsFID.write("getenv = MOVESFTS_PATH\n")
    MoveSFTsFID.write("queue 1\n")

# create the DAG file with the jobs to run
with open(args.dag_file, "w") as dagFID:
    startTimeAllNodes = None
    firstSFTstartTime = 0  # need this for the synchronized start option

    # Loop over the segment list to generate the SFTs for each segment
    for seg in segList:
        # Each segment in the segList runs on one or more nodes;
        # initialize the number SFTs produced by the current node:
        numThisNode = 0
        numThisSeg = 0

        # Case 1: a segment file was given but the SFTs do not need their
        # start times to be synchronized
        if adjustSegExtraTime and not args.synchronize_start:
            segStartTime = seg[0]
            segEndTime = seg[1]

            # First we figure out how much extra time is in the segment so that
            # SFTs are fit within the segment:
            # |..<SFT><SFT><SFT>..|
            # where the .. represent the extra time in the segment
            # The amount of extra time in a segment is given as the remainder
            # of (total segment time) / (SFT time baseline)
            segExtraTime = (segEndTime - segStartTime) % args.time_baseline

            # If there is overlap of SFTs requested, then we compute the extra
            # time as:
            # the remainder of (end - start - Tsft) / (non-overlap time)
            # provided there was at least one SFT that is in the segment
            if args.overlap_fraction != 0.0:
                if (segEndTime - segStartTime) > args.time_baseline:
                    segExtraTime = (
                        segEndTime - segStartTime - args.time_baseline
                    ) % int((1.0 - args.overlap_fraction) * args.time_baseline)

            # We'll add half the extra time to the start of the SFTs to be
            # created in this segment and half at the end
            segExtraStart = int(segExtraTime / 2)
            segExtraEnd = segExtraTime - segExtraStart
            args.analysis_start_time = segStartTime + segExtraStart

            # This shift may have pushed past the end time of the segment. In
            # that case, just fix the start time to the end time of the segment
            if args.analysis_start_time > segEndTime:
                args.analysis_start_time = segEndTime

            # shifting the end time by the other portion of the extra time
            # amount ...
            args.analysis_end_time = segEndTime - segExtraEnd

            # Again, this shift could have pushed the end time beyond the start
            # of the segment, so just fix the end time to the segment start
            if args.analysis_end_time < segStartTime:
                args.analysis_end_time = segStartTime

        # Case 2: SFTs need a synchronized start. This is a special case for
        # methods like TwoSpect, where signal periodicity spacing must be
        # maintained
        elif args.synchronize_start:
            segStartTime = seg[0]
            segEndTime = seg[1]

            # If we haven't set the first SFT start time, then set it equal to
            # the start time of the first segment
            if firstSFTstartTime == 0:
                firstSFTstartTime = segStartTime

            # This is a tricky bit of math to set the start time based on when
            # the first SFT start time of all the segments
            args.analysis_start_time = (
                int(
                    round(
                        math.ceil(
                            (segStartTime - firstSFTstartTime)
                            / ((1.0 - args.overlap_fraction) * args.time_baseline)
                        )
                        * (1.0 - args.overlap_fraction)
                        * args.time_baseline
                    )
                )
                + firstSFTstartTime
            )

            # This shift may have pushed past the end time of the segment. In
            # that case, just fix the start time to the end time of the segment
            if args.analysis_start_time > segEndTime:
                args.analysis_start_time = segEndTime

            # This is a tricky bit of math to set the end time based on when
            # the first SFT start time of all the segments
            args.analysis_end_time = (
                int(
                    round(
                        math.floor(
                            (segEndTime - args.analysis_start_time - args.time_baseline)
                            / ((1.0 - args.overlap_fraction) * args.time_baseline)
                        )
                        * (1.0 - args.overlap_fraction)
                        * args.time_baseline
                    )
                )
                + args.time_baseline
                + args.analysis_start_time
            )

            # Again, this shift could have pushed the end time beyond the start
            # of the segment, so just fix the end time to the segment start
            if args.analysis_end_time < segStartTime:
                args.analysis_end_time = segStartTime

        # If no segment file given and no synchronized starts, just set the
        # start time and end time to the segment start and end
        else:
            args.analysis_start_time = seg[0]
            args.analysis_end_time = seg[1]

        # Loop through the analysis time; make sure no more than
        # args.max_num_per_node SFTs are produced by any one node
        startTimeThisNode = args.analysis_start_time
        endTimeThisNode = args.analysis_start_time
        endTimeAllNodes = args.analysis_start_time
        while endTimeAllNodes < args.analysis_end_time:
            # increment endTimeAllNodes by the args.time_baseline until we get
            # past the args.analysis_end_time
            if args.overlap_fraction != 0.0:
                # handle overlap
                if numThisSeg == 0:
                    endTimeAllNodes = endTimeAllNodes + args.time_baseline
                else:
                    endTimeAllNodes = endTimeAllNodes + int(
                        (1.0 - args.overlap_fraction) * args.time_baseline
                    )
            else:
                # default case, no overlap
                endTimeAllNodes = endTimeAllNodes + args.time_baseline
            if endTimeAllNodes <= args.analysis_end_time:
                # increment the number of SFTs output from this node, and
                # update the end time this node.
                numThisNode = numThisNode + 1
                numThisSeg = numThisSeg + 1
                endTimeThisNode = endTimeAllNodes
                if numThisNode < args.max_num_per_node:
                    continue
                else:
                    # write jobs to dag for this node
                    nodeCount = nodeCount + 1

                    if nodeCount == 1:
                        startTimeAllNodes = startTimeThisNode
                    writeToDag(
                        dagFID,
                        nodeCount,
                        startTimeThisNode,
                        endTimeThisNode,
                        urls,
                        args,
                    )
                    # Update for next node
                    numThisNode = 0
                    if args.overlap_fraction != 0.0:
                        # handle overlap
                        startTimeThisNode = endTimeThisNode - int(
                            (args.overlap_fraction) * args.time_baseline
                        )
                    else:
                        # default case, no overlap
                        startTimeThisNode = endTimeThisNode
        else:
            # we are at or past the args.analysis_end_time; output job for last
            # node if needed.
            if numThisNode > 0:
                # write jobs to dag for this node
                nodeCount = nodeCount + 1

                if nodeCount == 1:
                    startTimeAllNodes = startTimeThisNode
                writeToDag(
                    dagFID, nodeCount, startTimeThisNode, endTimeThisNode, urls, args
                )
        # END while (endTimeAllNodes < args.analysis_end_time)
    # END for seg in segList

    # Write the move SFTs job to the DAG
    # Move SFTs
    dagFID.write(f"JOB MoveSFTs {Path(dagFID.name).parent / 'MoveSFTs.sub'}\n")
    dagFID.write(f"RETRY MoveSFTs 1\n")
    dagFID.write(f"VARS MoveSFTs ")
    if not args.validate:
        dagFID.write('opts="--no-validate" ')
    dagFID.write(
        f'sourcedirectory="." '
        f"channels=\"{' '.join(args.channel_name)}\" "
        f"destdirectory=\"{' '.join([str(p) for p in args.output_sft_path])}\"\n"
    )
    dagFID.write(
        f"PARENT {' '.join([f'MakeSFTs_{n}' for n in range(1, nodeCount+1)])} CHILD MoveSFTs\n"
    )

# Close the DAG file

# Update actual end time of the last job and print out the times all jobs will run on:
endTimeAllNodes = endTimeThisNode

if startTimeAllNodes is None:
    raise Exception("The startTimeAllNodes == none; the DAG file contains no jobs!")

if endTimeAllNodes <= startTimeAllNodes:
    raise Exception(
        "The endTimeAllNodes <= startTimeAllNodes; the DAG file contains no jobs!"
    )

print(startTimeAllNodes, endTimeAllNodes)
