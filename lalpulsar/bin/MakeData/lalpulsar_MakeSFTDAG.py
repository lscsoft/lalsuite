import math
import argparse
import os


# MakeSFTDAG.py - Creates DAGs to run jobs that generates SFTs
# can act as a dag generator for use with onasys.
# Loosely based on lalapps_strain_pipe


__author__ = 'Greg Mendell<gmendell@ligo-wa.caltech.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

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

#
# FUNCTION THAT WRITE ONE JOB TO DAG FILE
#
def writeToDag(dagFID, nodeCount, startTimeThisNode, endTimeThisNode, site, args):
    LSCdataFind = f'LSCdataFind_{nodeCount}'
    MakeSFTs = f'MakeSFTs_{nodeCount}'
    startTimeDatafind = startTimeThisNode - args.extra_datafind_time
    endTimeDatafind = endTimeThisNode + args.extra_datafind_time
    tagStringOut = f'{args.tag_string}_{nodeCount}'
    cacheFile = f'{args.cache_path}/{site}-{startTimeDatafind}-{endTimeDatafind}.cache'

    argList = []
    argList.append(f'-f {args.filter_knee_freq}')
    argList.append(f'-t {args.time_baseline}')
    argList.append(f'-p {args.output_sft_path}')
    argList.append(f'-C {cacheFile}')
    argList.append(f'-s {startTimeThisNode}')
    argList.append(f'-e {endTimeThisNode}')
    argList.append(f'-N {args.channel_name}')
    argList.append(f'-F {args.start_freq}')
    argList.append(f'-B {args.band}')
    if args.comment_field:
        argList.append(f'-c {args.comment_field}')
    if args.window_type:
        argList.append(f'-w {args.window_type}')
    if args.overlap_fraction:
        argList.append(f'-P {args.overlap_fraction}')
    argStr = ' '.join(argList)

    dagFID.write(f'JOB {LSCdataFind} datafind.sub\n')
    dagFID.write(f'RETRY {LSCdataFind} 10\n')
    dagFID.write(f'VARS {LSCdataFind} gpsstarttime="{startTimeDatafind}" gpsendtime="{endTimeDatafind}" observatory="{site}" inputdatatype="{args.input_data_type}" tagstring="{tagStringOut}"\n')
    dagFID.write(f'JOB {MakeSFTs} MakeSFTs.sub\n')
    dagFID.write(f'RETRY {MakeSFTs} 5\n')
    dagFID.write(f'VARS {MakeSFTs} argList="{argStr}" tagstring="{tagStringOut}"\n')
    dagFID.write(f'PARENT {LSCdataFind} CHILD {MakeSFTs}\n')

#
# MAIN CODE START HERE
#

parser = argparse.ArgumentParser(
    description='This script creates datafind.sub, MakeSFTs.sub, and a dag \
    file that generates SFTs based on the options given.',
    fromfile_prefix_chars='@')
parser.add_argument('-a', '--analysis-start-time', type=int,
                    help='GPS start time of data from which to generate \
                    SFTs (optional and unused if a segment file is given)')
parser.add_argument('-b', '--analysis-end-time', type=int,
                    help='GPS end time of data from which to generate SFTs \
                    (optional and unused if a segment file is given)')
parser.add_argument('-f', '--dag-file', required=True, type=str,
                    help='filename for .dag file (should end in .dag)')
parser.add_argument('-G', '--tag-string', required=True, type=str,
                    help='tag string used in names of various files unique to \
                    jobs that will run under the DAG')
parser.add_argument('-d', '--input-data-type', required=True, type=str,
                    help='input data type for use with the LSCdataFind --type \
                    option')
parser.add_argument('-x', '--extra-datafind-time', type=int, default=0,
                    help='extra time to subtract/add from/to start/end time \
                    arguments of LSCdataFind')
parser.add_argument('-M', '--datafind-match', type=str,
                    help='string to use with the LSCdataFind --match option')
parser.add_argument('-y', '--synchronize-start', action='store_true',
                    help='synchronize the start times of the SFTs so that the \
                    start times are synchronized when there are gaps in the \
                    data')
parser.add_argument('-k', '--filter-knee-freq', required=True, type=int,
                    help='high pass filter knee frequency used on time domain \
                    data before generating SFTs')
parser.add_argument('-T', '--time-baseline', required=True, type=int,
                    help='time baseline of SFTs  (e.g., 60 or 1800 seconds)')
parser.add_argument('-p', '--output-sft-path', required=True, type=str,
                    help='path to output SFTs')
parser.add_argument('-C', '--cache-path', type=str, default='cache',
                    help='path to cache files that will be produced by \
                    LSCdataFind (default is $PWD/cache; this directory is \
                    created if it does not exist and must agree with that \
                    given in .sub files)')
parser.add_argument('-O', '--log-path', type=str, default='logs',
                    help='path to some log, output, and error files (default \
                    is $PWD/logs; this directory is created if it does not \
                    exist and should agree with that given in .sub files)')
parser.add_argument('-o', '--sub-log-path', type=str, default='logs',
                    help='path to log file to give in datafind.sub and \
                    MakeSFTs.sub (default is $PWD/logs; this directory must \
                    exist and usually should be under a local file system)')
parser.add_argument('-N', '--channel-name', required=True, type=str,
                    help='name of input time-domain channel to read from \
                    frames')
parser.add_argument('-c', '--comment-field', type=str,
                    help='comment for version 2 SFT header')
parser.add_argument('-F', '--start-freq', type=int, default=10,
                    help='start frequency of the SFTs')
parser.add_argument('-B', '--band', type=int, default=1990,
                    help='frequency band of the SFTs')
parser.add_argument('-w', '--window-type', type=int, choices=[0, 1, 2, 3],
                    default=1, help='type of windowing of time-domain to do \
                    before generating SFTs (0 = None, 1 = Tukey given by \
                    Matlab, 2 = Tukey given in \
                    lalapps/src/pulsar/make_sfts.c, 3 = Hann given by Matlab')
parser.add_argument('-P', '--overlap-fraction', type=float, default=0,
                    help='overlap fraction (for use with windows; e.g., use \
                    -P 0.5 with -w 3 Hann windows)')
parser.add_argument('-m', '--max-num-per-node', type=int, default=1,
                    help='maximum number of SFTs to generate on one node')
parser.add_argument('-L', '--max-length-all-jobs', type=int,
                    help='maximum total amount of data to process, in seconds \
                    (optional and unused if a segment file is given)')
parser.add_argument('-g', '--segment-file', type=str,
                    help='alternative file with segments to use, rather than \
                    the input times')
parser.add_argument('-l', '--min-seg-length', type=int, default=0,
                    help='minimum length segments to process in seconds (used \
                    only if a segment file is given)')
parser.add_argument('-q', '--list-of-nodes', type=str,
                    help='file with list of nodes on which to output SFTs')
parser.add_argument('-Q', '--node-path', type=str,
                    help='path to nodes to output SFTs; the node name is \
                    appended to this path, followed by path given by the -p \
                    option; for example, if -q point to file with the list \
                    node1 node2 ... and the -Q /data/ -p /frames/S5/sfts/LHO \
                    options are given, the first output file will go into \
                    /data/node1/frames/S5/sfts/LHO; the next node in the list \
                    is used in constructing the path when the number of jobs \
                    given by the -R option reached, and so on')
parser.add_argument('-R', '--output-jobs-per-node', type=int, default=0,
                    help='number of jobs to output per node in the list of \
                    nodes given with the -q option')
parser.add_argument('-j', '--datafind-path', type=str,
                    help='string specifying a path to look for the \
                    gw_data_find executable; if not set, will use \
                    LSC_DATAFIND_PATH env variable or system default (in \
                    that order)')
parser.add_argument('-J', '--makesfts-path', type=str,
                    help='string specifying a path to look for the \
                    lalpulsar_MakeSFTs executable; if not set, will use \
                    MAKESFTS_PATH env variable or system default (in that \
                    order)')
parser.add_argument('-Y', '--request-memory', type=int, default=2048,
                    help='memory allocation in MB to request from condor for \
                    lalpulsar_MakeSFTs step')
parser.add_argument('-s', '--request-disk', type=int, default=1024,
                    help='disk space allocation in MB to request from condor \
                    for lalpulsar_MakeSFTs step')
parser.add_argument('-A', '--accounting-group', required=True, type=str,
                    help='Condor tag for the production of SFTs')
parser.add_argument('-U', '--accounting-group-user', required=True, type=str,
                    help='albert.einstein username (do not add @LIGO.ORG)')


##### DEPRECATED OPTIONS #####
class DeprecateAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        parser.error(f'Argument {self.option_strings} has been deprecated in lalpulsar_MakeSFTs')

parser.add_argument('-u', '--frame-struct-type', nargs=0, action=DeprecateAction,
                    help='DEPRECATED. No longer required; \
                    the frame channel type is determined automatically')
parser.add_argument('-H', '--use-hot', nargs=0, action=DeprecateAction,
                    help='DEPRECATED. No longer required; \
                    the frame channel type is determined automatically')
parser.add_argument('-i', '--ifo', nargs=0, action=DeprecateAction,
                    help='DEPRECATED. No longer required; \
                    the detector prefix is deduced from the channel name')
parser.add_argument('-D', '--make-gps-dirs', nargs=0, action=DeprecateAction,
                    help='DEPRECATED. No longer supported')
parser.add_argument('-Z', '--make-tmp-file', nargs=0, action=DeprecateAction,
                    help='DEPRECATED. Default behaviour')
parser.add_argument('-X', '--misc-desc', nargs=0, action=DeprecateAction,
                    help='DEPRECATED. No longer supported \
                    in the public SFT filename specification')
parser.add_argument('-v', '--sft-version', nargs=0, action=DeprecateAction,
                    help='DEPRECATED. No longer supported')
parser.add_argument('-S', '--use-single', nargs=0, action=DeprecateAction,
                    help='DEPRECATED. No longer supported')


args = parser.parse_args()


# Some basic argument value checking
if args.extra_datafind_time < 0:
    raise argparse.error('--extra-datafind-time must be >= 0')

if args.filter_knee_freq < 0:
    raise argparse.error('--filter-knee-freq must be >= 0')

if args.time_baseline <= 0:
    raise argparse.error('--time-baseline must be > 0')

if args.overlap_fraction < 0.0 or args.overlap_fraction >= 1.0:
    raise argparse.error('--overlap-fraction must be in the range [0,1)')

if args.start_freq < 0.0 or args.start_freq >= 7192.0:
    raise argparse.error('--start-freq must be in the range [0,7192)')

if args.band <= 0 or args.band >= 8192.0:
    raise argparse.error('--band must be in the range (0,8192)')

if args.start_freq + args.band >= 8192.0:
    raise argparse.error('--start-freq + --band must be < 8192')

if args.max_num_per_node <= 0:
    raise argparse.error('--max-num-per-node must be > 0')

# Set the data find executable and lalpulsar_MakeSFTs executable
dataFindExe = 'gw_data_find'
if args.datafind_path:
    dataFindExe = os.path.join(args.datafind_path, dataFindExe)
elif 'LSC_DATAFIND_PATH' in os.environ:
    dataFindExe = os.path.join('$ENV(LSC_DATAFIND_PATH)', dataFindExe)
else:
    dataFindExe = os.path.join('/usr/bin', dataFindExe)

makeSFTsExe = 'lalpulsar_MakeSFTs'
if args.makesfts_path:
    makeSFTsExe = os.path.join(args.makesfts_path, makeSFTsExe)
elif 'MAKESFTS_PATH' in os.environ:
    makeSFTsExe = os.path.join('$ENV(MAKESFTS_PATH)', makeSFTsExe)
else:
    makeSFTsExe = os.path.join('@LALSUITE_BINDIR@', makeSFTsExe)

# try and make a directory to store the cache files and job logs
try: os.mkdir(args.log_path)
except: pass
try: os.mkdir(args.cache_path)
except: pass

# Check if list of nodes is given, on which to output SFTs.
nodeList = []
useNodeList = False
savedOutputSFTPath = None
if args.list_of_nodes is not None:

    if args.node_path is None:
        raise argparse.error(
            'Node file list given, but no node path specified')

    if args.output_jobs_per_node < 1:
        raise argparse.error(
            'Node file list given, but invalid output jobs per node specified')

    with open(args.list_of_nodes) as fp_nodelist:
        for idx, line in enumerate(fp_nodelist):
            splitLine = line.split()
            nodeList.append(splitLine[0])
        if len(nodeList) < 1:
            raise ValueError('No nodes found in node list file: {}'.format(
                args.list_of_nodes))

    # Set flag to use list of nodes in constructing output files
    useNodeList = True
    savedOutputSFTPath = args.output_sft_path
# END if (args.list_of_nodes != None)

# Check if segment file was given, else set up one segment from the command line
segList = []
adjustSegExtraTime = False
if args.segment_file is not None:

    if args.min_seg_length < 0:
        raise argparse.error('--min-seg-length must be >= 0')

    # the next flag causes extra time that cannot be processes to be trimmed
    # from the start and end of a segment
    adjustSegExtraTime = True

    with open(args.segment_file) as fp_segfile:
        for idx, line in enumerate(fp_segfile):
            splitLine = line.split()
            oneSeg = []
            oneSeg.append(int(splitLine[0]))
            oneSeg.append(int(splitLine[1]))
            if (oneSeg[1] - oneSeg[0]) >= args.min_seg_length:
                segList.append(oneSeg)

    if len(segList) < 1:
        raise ValueError('No segments found in segment file: {}'.format(
            args.segment_file))
else:
    if args.analysis_start_time is None:
        raise argparse.error(
            '--analysis-start-time must be specified if no segment file is \
            given')

    if args.analysis_end_time is None:
        raise argparse.error(
            '--analysis-start-time must be specified if no segment file is \
            given')

    if args.max_length_all_jobs is None:
        raise argparse.error(
            '--max-length-all-jobs must be specified if no segment file is \
            given')

    # Make sure not to exceed maximum allow analysis
    if args.analysis_end_time > (args.analysis_start_time +
                                 args.max_length_all_jobs):
        args.analysis_end_time = args.analysis_start_time + \
          args.max_length_all_jobs

    oneSeg = [];
    oneSeg.append(args.analysis_start_time);
    oneSeg.append(args.analysis_end_time);
    segList.append(oneSeg);
# END if (args.segment_file != None)

# Get the IFO site, which is the first letter of the channel name.
site = args.channel_name[0]

# initialize count of nodes
nodeCount = 0

# Create .sub files
path_to_dag_file = os.path.dirname(args.dag_file)
dag_filename = os.path.basename(args.dag_file)
datafind_sub = os.path.join(path_to_dag_file, 'datafind.sub')
makesfts_sub = os.path.join(path_to_dag_file, 'MakeSFTs.sub')

# create datafind.sub
with open(datafind_sub, 'w') as datafindFID:
    datafindLogFile = '{}/datafind_{}.log'.format(args.sub_log_path,
                                                  dag_filename)
    datafindFID.write('universe = vanilla\n')
    datafindFID.write('executable = {}\n'.format(dataFindExe))
    if not args.datafind_match:
        dataFindMatchString = ''
    else:
        dataFindMatchString = '--match {}'.format(args.datafind_match)
    datafindFID.write('arguments = -r $ENV(LIGO_DATAFIND_SERVER) ')
    datafindFID.write('--observatory $(observatory) --url-type file ')
    datafindFID.write('--gps-start-time $(gpsstarttime) ')
    datafindFID.write('--gps-end-time $(gpsendtime) --lal-cache ')
    datafindFID.write('--type $(inputdatatype) {}\n'.format(
        dataFindMatchString))
    datafindFID.write('getenv = True\n')
    datafindFID.write('request_disk = 5MB\n')
    datafindFID.write('accounting_group = {}\n'.format(args.accounting_group))
    datafindFID.write('accounting_group_user = {}\n'.format(
        args.accounting_group_user))
    datafindFID.write('log = {}\n'.format(datafindLogFile))
    datafindFID.write('error = {}/datafind_$(tagstring).err\n'.format(
        args.log_path))
    datafindFID.write('output = {}/'.format(args.cache_path))
    datafindFID.write('$(observatory)-$(gpsstarttime)-$(gpsendtime).cache\n')
    datafindFID.write('notification = never\n')
    datafindFID.write('queue 1\n')

# create MakeSFTs.sub
with open(makesfts_sub, 'w') as MakeSFTsFID:
    MakeSFTsLogFile = '{}/MakeSFTs_{}.log'.format(args.sub_log_path,
                                                  dag_filename)
    MakeSFTsFID.write('universe = vanilla\n')
    MakeSFTsFID.write('executable = {}\n'.format(makeSFTsExe))
    MakeSFTsFID.write('arguments = $(argList)\n')
    MakeSFTsFID.write('getenv = True\n')
    MakeSFTsFID.write('accounting_group = {}\n'.format(args.accounting_group))
    MakeSFTsFID.write('accounting_group_user = {}\n'.format(
        args.accounting_group_user))
    MakeSFTsFID.write('log = {}\n'.format(MakeSFTsLogFile))
    MakeSFTsFID.write('error = {}/MakeSFTs_$(tagstring).err\n'.format(
        args.log_path))
    MakeSFTsFID.write('output = {}/MakeSFTs_$(tagstring).out\n'.format(
        args.log_path))
    MakeSFTsFID.write('notification = never\n')
    MakeSFTsFID.write(f'request_memory = {args.request_memory}MB\n')
    MakeSFTsFID.write(f'request_disk = {args.request_disk}MB\n')
    MakeSFTsFID.write('RequestCpus = 1\n')
    MakeSFTsFID.write('queue 1\n')

# create the DAG file with the jobs to run
with open(args.dag_file, 'w') as dagFID:
    startTimeAllNodes = None
    firstSFTstartTime = 0 # need this for the synchronized start option
    nodeListIndex = 0

    # Loop over the segment list to generate the SFTs for each segment
    for seg in segList:
        # Each segment in the segList runs on one or more nodes;
        # initialize the number SFTs produced by the current node:
        numThisNode = 0
        numThisSeg = 0

        # Case 1: a segment file was given but the SFTs do not need their
        # start times to be synchronized
        if (adjustSegExtraTime and not args.synchronize_start):
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
                    segExtraTime = \
                      (segEndTime - segStartTime - args.time_baseline) \
                      % int((1.0 - args.overlap_fraction)*args.time_baseline)

            # We'll add half the extra time to the start of the SFTs to be
            # created in this segment and half at the end
            segExtraStart =  int(segExtraTime / 2)
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
        elif (args.synchronize_start):
            segStartTime = seg[0]
            segEndTime = seg[1]

            # If we haven't set the first SFT start time, then set it equal to
            # the start time of the first segment
            if firstSFTstartTime == 0:
                firstSFTstartTime = segStartTime

            # This is a tricky bit of math to set the start time based on when
            # the first SFT start time of all the segments
            args.analysis_start_time = \
              int(round(math.ceil((segStartTime - firstSFTstartTime)/((
                  1.0 - args.overlap_fraction)*args.time_baseline))*(
                      1.0 - args.overlap_fraction)*args.time_baseline)) + \
                      firstSFTstartTime

            # This shift may have pushed past the end time of the segment. In
            # that case, just fix the start time to the end time of the segment
            if args.analysis_start_time > segEndTime:
                args.analysis_start_time = segEndTime

            # This is a tricky bit of math to set the end time based on when
            # the first SFT start time of all the segments
            args.analysis_end_time = \
              int(round(math.floor((
                  segEndTime - args.analysis_start_time -
                  args.time_baseline)/((
                      1.0 - args.overlap_fraction)*args.time_baseline))*\
                            (1.0 - args.overlap_fraction)*\
                            args.time_baseline)) + args.time_baseline + \
                            args.analysis_start_time

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
        endTimeThisNode   = args.analysis_start_time
        endTimeAllNodes   = args.analysis_start_time
        while (endTimeAllNodes < args.analysis_end_time):
            # increment endTimeAllNodes by the args.time_baseline until we get
            # past the args.analysis_end_time
            if args.overlap_fraction != 0.0:
                # handle overlap
                if numThisSeg == 0:
                    endTimeAllNodes = endTimeAllNodes + args.time_baseline
                else:
                    endTimeAllNodes = endTimeAllNodes + \
                      int((1.0 - args.overlap_fraction)*args.time_baseline)
            else:
                # default case, no overlap
                endTimeAllNodes = endTimeAllNodes + args.time_baseline
            if (endTimeAllNodes <= args.analysis_end_time):
                # increment the number of SFTs output from this node, and
                # update the end time this node.
                numThisNode = numThisNode + 1
                numThisSeg = numThisSeg + 1
                endTimeThisNode = endTimeAllNodes
                if (numThisNode < args.max_num_per_node):
                    continue
                else:
                    # write jobs to dag for this node
                    nodeCount = nodeCount + 1

                    if (useNodeList):
                        args.output_sft_path = args.node_path + \
                          nodeList[nodeListIndex] + savedOutputSFTPath
                        if ((nodeCount % args.output_jobs_per_node) == 0):
                            nodeListIndex = nodeListIndex + 1
                        # END if ((nodeCount % args.output_jobs_per_node) == 0L)
                    # END if (useNodeList)

                    if (nodeCount == 1):
                        startTimeAllNodes = startTimeThisNode
                    writeToDag(dagFID, nodeCount, startTimeThisNode, endTimeThisNode, site, args)
                    # Update for next node
                    numThisNode       = 0
                    if args.overlap_fraction != 0.0:
                        # handle overlap
                        startTimeThisNode = endTimeThisNode - \
                          int((args.overlap_fraction)*args.time_baseline)
                    else:
                        # default case, no overlap
                        startTimeThisNode = endTimeThisNode
        else:
            # we are at or past the args.analysis_end_time; output job for last
            # node if needed.
            if (numThisNode > 0):
                # write jobs to dag for this node
                nodeCount = nodeCount + 1

                if (useNodeList):
                    args.output_sft_path = args.node_path + \
                      nodeList[nodeListIndex] + savedOutputSFTPath
                    if ((nodeCount % args.output_jobs_per_node) == 0):
                        nodeListIndex = nodeListIndex + 1
                    # END if ((nodeCount % args.output_jobs_per_node) == 0L)
                # END if (useNodeList)

                if (nodeCount == 1):
                    startTimeAllNodes = startTimeThisNode
                writeToDag(dagFID, nodeCount, startTimeThisNode, endTimeThisNode, site, args)
        # END while (endTimeAllNodes < args.analysis_end_time)
    # END for seg in segList
# Close the DAG file

# Update actual end time of the last job and print out the times all jobs will run on:
endTimeAllNodes = endTimeThisNode

if startTimeAllNodes is None:
  raise Exception('The startTimeAllNodes == none; \
  the DAG file contains no jobs!')

if (endTimeAllNodes <= startTimeAllNodes):
  raise Exception('The endTimeAllNodes <= startTimeAllNodes; \
  the DAG file contains no jobs!')

print(startTimeAllNodes, endTimeAllNodes)

