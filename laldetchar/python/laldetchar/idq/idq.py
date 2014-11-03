# Copyright (C) 2013 Lindy Blackburn, Reed Essick and Ruslan Vaulin.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## addtogroup pkg_py_laldetchar_idq
## Synopsis
# ~~~
# from laldetchar.idq import idq
# ~~~
# \author Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)

"""Module with basic functions and classes of idq pipeline.
"""

import glob
import os
import re
import sys
import time
import numpy
import math
from laldetchar.idq import event
from laldetchar.idq import ovl
from laldetchar.idq import auxmvc_utils
from laldetchar.idq import auxmvc
from laldetchar.idq import idq_tables
from pylal import frutils
import subprocess
import tempfile
import logging

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import ilwd

from laldetchar import git_version

__author__ = \
    'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

## addtogroup pkg_py_laldetchar_idq_idq
# @{

# common routines for idq codes
# idq-specific rountines should go here
# generic event processing should go in event.py

#=================================================
# timing utilities
#=================================================
gpsref = time.mktime(time.strptime('Tue Jan 06 00:00:00 1980')) \
    - time.timezone  # GPS ref in UNIX time


### return current GPS time, assume 16 leap-seconds
def nowgps():
    return time.time() - gpsref + 16  # add in leap-seconds

def cluster(
    glitches,
    columns,
    cluster_key='signif',
    cluster_window=1.0,
    ):
    """
        Clustering performed with the sliding window cluster_window keeping trigger with the highest rank;
        glitches is an array of glitches, and we sort it to be in ascending order (by GPS time)
        Clustering algorithm is borrowed from pylal.CoincInspiralUtils cluster method.
        """

    # sort glitches so they are in the propper order

    glitches.sort(key=lambda line: line[columns['GPS']])

    # initialize some indices (could work with just one)
    # but with two it is easier to read

    this_index = 0
    next_index = 1
    while next_index < len(glitches):

        # get the time for both indices

        thisTime = glitches[this_index][columns['GPS']]
        nextTime = glitches[next_index][columns['GPS']]

        # are the two coincs within the time-window?

        if nextTime - thisTime < cluster_window:

            # get the ranks

            this_rank = glitches[this_index][columns[cluster_key]]
            next_rank = glitches[next_index][columns[cluster_key]]

            # and remove the trigger which has the lower rank

            if next_rank > this_rank:
                del glitches[this_index]  # glitches = numpy.delete(glitches, this_index, 0)
            else:
                del glitches[next_index]  # glitches = numpy.delete(glitches, next_index, 0)
        else:

            # NOTE: we don't increment this_index and next_index if we've removed glitches
            # the two triggers are NOT in the time-window
            # so must increase index

            this_index += 1
            next_index += 1

    return glitches

#=================================================
# lock-file utilities
#=================================================
### try to set lockfile, die if cannot establish lock
def dieiflocked(lockfile='.idq.lock'):
    import fcntl
    global lockfp
    lockfp = open(lockfile, 'w')
    try:
        fcntl.lockf(lockfp, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except IOError:
        import sys
        sys.exit('ERROR: cannot establish lock on \'' + lockfile
                 + '\', possible duplicate process, exiting..')


def findCredential():
    """
....This function is a modifcation of frutils.findCredential() function to adopt to a new file naming convention for ligo
....certificates. Depricate it when the orginal function is modified accordingly. 
...."""

    rfc_proxy_msg = \
        """\
	Could not find a RFC 3820 compliant proxy credential.
Please run 'grid-proxy-init -rfc' and try again.
"""

    # use X509_USER_PROXY from environment if set

    if os.environ.has_key('X509_USER_PROXY'):
        filePath = os.environ['X509_USER_PROXY']
        try:
            if frutils.validateProxy(filePath):
                return (filePath, filePath)
        except:
            pass

        # ....raise RuntimeError(rfc_proxy_msg)

    # use X509_USER_CERT and X509_USER_KEY if set

    if os.environ.has_key('X509_USER_CERT'):
        if os.environ.has_key('X509_USER_KEY'):
            certFile = os.environ['X509_USER_CERT']
            keyFile = os.environ['X509_USER_KEY']
            return (certFile, keyFile)

    # search for proxy file on disk

    uid = os.getuid()

    # path = "/tmp/x509up_u%d" % uid

    certFiles = [path for path in glob.glob('/tmp/x509*')
                 if os.stat(path).st_uid == uid]
    certFiles.sort(key=lambda x: os.stat(x).st_mtime)

    # return the most recent valid certificate

    for path in reversed(certFiles):
        if os.access(path, os.R_OK):
            try:
                if frutils.validateProxy(path):
                    return (path, path)
            except:
                pass

    # if we get here could not find a credential

    raise RuntimeError(rfc_proxy_msg)


#=================================================
# logfile object for use with logging module
#=================================================
class LogFile(object):

    """File-like object to log text using the `logging` module."""

    def __init__(self, logger):
        self.logger = logger

    def write(self, msg, level=logging.INFO):
        msg = msg.strip('\n').strip()
        if msg and msg != '\n' and msg != '':
            self.logger.log(level, msg)

    def flush(self):
        for handler in self.logger.handlers:
            handler.flush()

#=================================================
# logfile status and file manipulation, general I/O
#=================================================
def start_of_line(file_obj, byte_stride=80):
    """
    moves back one line in the file and returns both that line and the number of bytes moved
    if you are the beginning of a line, this returns "" and leaves you where you are
    otherwise this moves you to the beginning of the line
    """
    current_byte = file_obj.tell()
    tell = file_obj.tell()
    while tell > 0:
        if tell < byte_stride: ### can't go before the beginning
            stride = tell
        else:
            stride = byte_stride

        ### move back one stride
        file_obj.seek(-stride, 1)

        ### read that stride
        line = file_obj.read(stride)

        ### is there at least one line break
        if "\n" in line: # yes!
            file_obj.seek( -len(line.split("\n")[-1]), 1) ### move back just enough to put you after the line break
            tell = file_obj.tell() ### update current position
            break ### exit loop
        else:
            ### move all back again and continue to next epoch
            file_obj.seek(-stride, 1)

        tell = file_obj.tell() ### update current position

    return current_byte - tell ### return the number of bytes we moved back

###
def most_recent_realtime_stride(file_obj):
    """
    works backward through file_obj until it finds realtime-stride information
    returns the stride information (start, stop) and returns to the position in the file at when
    this process was launched
    """
    starting_tell = file_obj.tell()
    tell = starting_tell
    bytes = 0
	
    while tell > 0: ### not at the beginning of the file
        ### move to the start of this line
        ### blocksize is the number of bytes moved
        blocksize = start_of_line(file_obj)
        bytes += blocksize ### count the number of bytes we moved back
        tell = file_obj.tell()

        ### read this line
        line = file_obj.readline()

        ### parse line for realtime stride information
        if "Begin: stride " in line: ### FIXME? this is fragile and could break if we change the logger statements
            ### find the stride start
            stride_start, stride_end = [float(l) for l in  line.strip().split("Begin: stride ")[-1].split("-")]
            break

        ### if we're still in the loop, we must have NOT found a summary stride
        if tell: ### we did not move to the beginning of the file
            ### move back to start of this line
            file_obj.seek(-(len(line)+1), 1)
            bytes += 1
            tell = file_obj.tell()

#        else:
#            file_obj.seek(0,0) ### we're at the beginnig of the file

    else: ### reached beginning of file without finding realtime stride. We now look forward
        stride_start = -np.infty ### we wait until we find a new stride or time out
        stride_end = -np.infty

    file_obj.seek(starting_tell, 0) ### go back to where we started

    return (stride_start, stride_end)

###
def block_until(t, file_obj, max_wait=600, timeout=600, wait=0.1):
    """
    this function blocks until the stride_start information from file_obj (which should be an open realtime logger file object) shows that the realtime process has passed "t"
    we wait a maximum of "max_wait" seconds for the realtime job to pass this t (only checked after we determine the first realtime_stride)
    we wait at most "timeout" seconds between lines printed to file_obj, otherwise we assume the job has died
    between requesting new lines from file_obj, we wait "wait" seconds
    """
    file_obj.seek(0,2) ### go to end of file

    stride_start, _ = most_recent_realtime_stride(file_obj) ### delegate to find most recent stride

    to = time.time()
    waited = 0.0
    while (stride_start < t) and (waited < max_wait) and (time.time()-to < timeout):
        line = file_obj.readline() ### get next line

        if not line: ### nothing new to report
            time.sleep( wait )
            waited += wait
            continue

        elif "Begin: stride " in line: ### FIXME? this is fragile and could break if we change the logger statements
            ### find the stride start
            stride_start = float( line.split("Begin: stride ")[-1].split("-")[-1] )

        waited = 0.0 ### we found a new line, so re-set this counter

    ### return statement conditioned on how the loop exited
    return stride_start > t, waited > max_wait, time.time()-to > timeout


###
def get_condor_dag_status(dag):
    """
....Check on status of the condor dag by parsing the last line of .dagman.out file.
...."""

    dag_out_file = dag + '.dagman.out'
    if not os.path.exists(dag_out_file):
        exit_status = 'incomplete'
        return exit_status
    lastline = open(dag_out_file, 'r').readlines()[-1].strip('\n')
    if 'EXITING WITH STATUS' in lastline:
        exit_status = int(lastline.split('EXITING WITH STATUS')[-1])
    else:
        exit_status = 'incomplete'

    return exit_status

# load *.dat files
def slim_load_datfile(file, skip_lines=1, columns=[]):
    """ 
....loads only the given columns from a *.dat file. assumes the variable names are given by the first line after skip_lines 
....left so that if an element in columns is not in the variables list, an error will be raised
...."""

    output = dict([(c, []) for c in columns])

    f_lines = open(file).readlines()[skip_lines:]

    if f_lines:

        # find variable names

        variables = dict([(line, i) for (i, line) in
                         enumerate(f_lines[0].strip('\n').split())])

        # fill in output

        for line in f_lines[1:]:
            if line[0] != '#' and line != '':
                line = line.strip('\n').split()
                for c in columns:
                    output[c].append(line[variables[c]])
        return output
    else:

        # file is empty, return empty dictionary

        return output

#=================================================
# subprocess wrappers
#=================================================
def submit_command(
    command,
    process_name='unspecified process',
    dir='.',
    verbose=False,
    ):
    """
....Utility function for submiting jobs via python's subprocess. 
....@param command is a list, first element of which is interpretted as an executable
....and all others as options and values for that executable.
....@param process_name is the name identifying the submitted job.
....@param dir is directory in which the job is to be executed. 
...."""

    process = subprocess.Popen(command, stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE, cwd=dir)
    (output, errors) = process.communicate()
    if process.returncode:
        print 'Error! ' + process_name + ' failed.'
        print 'Error! ' + process_name + ' ' + 'exit code: ' \
            + str(process.returncode)
        if errors:
            for line in errors.split('\n'):
                if line != '':
                    print 'Error! ' + process_name + ' error out: ' \
                        + line.rstrip('\n')
        if output:
            for line in output.split('\n'):
                if line != '':
                    print 'Error! ' + process_name + ' standard out: ' \
                        + line.rstrip('\n')
        exit_code = 1
        return exit_code
    else:
        exit_code = 0
    if verbose:
        for line in output.split('\n'):
            if line != '':
                print process_name + ': ' + line.rstrip('\n')
    return exit_code


def condor_command(dag, process_name):
    command = 'condor_submit_dag ' + dag.basename + '.dag'
    process = subprocess.Popen(command, shell=True,
                               stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE)  # , cwd=dir)
    (errors, output) = process.communicate()
    if process.returncode:
        print process_name + ' failed, Error:', [output]  # [process.returncode, errors, output]
        exit_code = 1
    else:
        exit_code = 0
    return exit_code


def execute_dag_create(cp, log_path, run_tag):

    return auxmvc.create_DAG(cp, log_path, run_tag)


def extract_dq_segments(xmlfile, dq_name):
    """
....Loads xmlfile containing dq segments and extracts segments of the type set by dq_name ( e.g. L1:DMT-SCIENCE:3).
....dq_name is what was given as an --include-segments option to ligolw_segment_query that generated the xmlfileobj"
...."""

    if type(xmlfile) == str:
        xmldoc = ligolw_utils.load_filename(xmlfile)  # load as filename
    else:
        xmldoc = ligolw_utils.load_fileobj(xmlfile)[0]  # laod as file object

    # get segment tables

    sdef = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
    ssum = table.get_table(xmldoc, lsctables.SegmentSumTable.tableName)
    seg = table.get_table(xmldoc, lsctables.SegmentTable.tableName)

    # segment definer ID corresponding to dq_name, which was use in segment data base query

    id = next(a.segment_def_id for a in sdef if a.name
              == dq_name.split(':')[1])

    # list of covered segments (e.g. where SCIENCE mode is defined)

    covered = [[a.start_time, a.end_time] for a in ssum
               if a.segment_def_id == id]

    # final RESULT segment from query

    result_id = next(a.segment_def_id for a in sdef if a.name
                     == 'RESULT')
    good = [[a.start_time, a.end_time] for a in seg if a.segment_def_id
            == result_id]
    return (good, covered)

#=================================================
# segment manipulations
#=================================================
def extract_lldq_segments(xmlfiles, lldq_name, hoft_name):
    """
....Loads the low-latency science and hoft segments from xml files.
....lldq_name is the segment definer name for the low-latency science segments,
....hoft_name is the segment definer name for the hoft segments.
...."""

    good = []
    covered = []
    for file in xmlfiles:
        xmldoc = ligolw_utils.load_filename(file)  # load file

        # get segment tables

        sdef = table.get_table(xmldoc,
                               lsctables.SegmentDefTable.tableName)
        seg = table.get_table(xmldoc, lsctables.SegmentTable.tableName)

        # segment definer ID corresponding to  hoft_name
        # FIXME: a.ifos should be replace by a.name when gstlal_segment_parser files are fixed

        hoft_seg_def_id = next(a.segment_def_id for a in sdef if a.name
                               == hoft_name)

        # get list of covered hoft segments (e.g. where SCIENCE mode is defined)

        covered.extend([[a.start_time, a.end_time] for a in seg
                       if a.segment_def_id == hoft_seg_def_id])

        # segmen definer ID correspodning to low-latency science segments
        # FIXME: a.ifos should be replace by a.name when gstlal_segment_parser files are fixed

        lldq_seg_def_id = next(a.segment_def_id for a in sdef if a.name
                               == lldq_name)

        # get list of science segments

        good.extend([[a.start_time, a.end_time] for a in seg
                    if a.segment_def_id == lldq_seg_def_id])

    return (good, covered)


def extract_dmt_segments(xmlfiles, dq_name):
    """
....Loads the segments from the dmt xml files. Determines segments covered by dmt process from the file names
....dq_name is the segment definer name.
...."""

    good = []
    covered = []
    for file in xmlfiles:
        xmldoc = ligolw_utils.load_filename(file)  # load file

        # get segment tables

        sdef = table.get_table(xmldoc,
                               lsctables.SegmentDefTable.tableName)
        seg = table.get_table(xmldoc, lsctables.SegmentTable.tableName)

        # segment definer ID correspodning to dq_name
        # FIX ME: in case of multiple segment versions this matching becomes ambiguous

        dq_seg_def_id = next(a.segment_def_id for a in sdef if a.name
                             == dq_name)

        # get list of  segments

        good.extend([[a.start_time, a.end_time] for a in seg
                    if a.segment_def_id == dq_seg_def_id])

    # coalesce segments....

    good = event.fixsegments(good)

    # convert file names into segments

    covered = map(filename_to_segment, xmlfiles)

    # coalesce segments

    covered = event.fixsegments(covered)

    return (good, covered)


def filename_to_segment(filename):
    """
....Get start and duration from the file name, and convert them into a segment
...."""

    start = int(os.path.split(filename)[1].split('.')[0].split('-')[-2])
    duration = int(os.path.split(filename)[1].split('.')[0].split('-'
                   )[-1])
    end = start + duration
    return [start, end]


def segment_query(
    config,
    start,
    end,
    url=None,
    ):
    """
....Performes ligolw_segment query for segments in the given gps time range and returns segment xml file.
...."""

# ....from glue.ligolw import table,lsctables,utils
    # ligolw_segment_query -t file:///gds-l1/dmt/triggers/DQ_Segments -q --include-segments=L1:DMT-SCIENCE -s 1042257690 -e 1042257715

    program = config.get('get_science_segments', 'program')
    if url is None:
        url = config.get('get_science_segments', 'xmlurl')
    include = config.get('get_science_segments', 'include')
    args = '%s -t %s -q -a %s -s %d -e %d' % (program, url, include,
            start, end)

    # query to database is faster
    # need to set export S6_SEGMENT_SERVER=https://segdb-er.ligo.caltech.edu
    # and then run with -d option instead of -t
    # args....= '%s -d -q -a %s -s %d -e %d' % (program, include, start, end)

    xmlfileobj = tempfile.TemporaryFile()
    print 'segment query: ', args.split()
    p = subprocess.Popen(args.split(), stdout=xmlfileobj,
                         stderr=sys.stderr)
    p.wait()
    xmlfileobj.seek(0)  # go back to beginning of file

    return xmlfileobj

def get_idq_segments(
    realtimedir,
    start,
    stop,
    suffix='.pat',
    ):
    import re
    matchfile = re.compile('.*-([0-9]*)-([0-9]*)\%s$' % suffix)

       # (dirname, starttime, endtime, pad=0, suffix='.xml')

    ls = get_all_files_in_range(realtimedir, start, stop, suffix=suffix)
    fileseg = []
    for file in ls:
        m = matchfile.match(file)
        (st, dur) = (int(m.group(1)), int(m.group(2)))
        fileseg.append([st, st + dur])
    fileseg = event.fixsegments(fileseg)
    return event.andsegments([[start, stop]], fileseg)

#=================================================
# trigger generation algorithm interfaces
#=================================================

# load in KW configuration file as dictionary
# note that nothing will be converted to numerical values, they will all be saved as string values
# also produce list of channel_flow_fhigh names
# original order of channels is preserved
def loadkwconfig(file):
    table = event.loadstringtable(file)
    kwconfig = dict()
    singleopt = \
        'stride basename segname threshold significance maxDist decimateFactor transientDuration'.split()
    multiopt = 'vetoseg rateflag channel'
    for line in table:
        (key, par) = (line[0], (line[1] if len(line)
                      == 2 else line[1:]))
        if key in singleopt:
            if key in kwconfig:
                raise '%s is defined twice in %s' % (key, file)
            kwconfig[key] = par
        if key in multiopt:
            if key in kwconfig:
                kwconfig[key].append(par)
            else:
                kwconfig[key] = [par]
    if 'channel' in kwconfig:
        kwconfig['names'] = ['_'.join(cpar).replace(':', '_')
                             for cpar in kwconfig['channel']]
    return kwconfig


#=================================================
# file discovery methods
#=================================================
# this is a modified (bugfix, efficiency fix, suffix option) version from the original in glue
def get_all_files_in_range(
    dirname,
    starttime,
    endtime,
    pad=0,
    suffix='.xml',
    ):
    """Returns all files in dirname and all its subdirectories whose
....names indicate that they contain segments in the range starttime
....to endtime"""

    ret = []
    matchfile = re.compile('.*-([0-9]*)-([0-9]*)\%s$' % suffix)

    # Maybe the user just wants one file...

    if os.path.isfile(dirname):
        m = matchfile.match(dirname)
        if m:
            (start, dur) = (int(m.group(1)), int(m.group(2)))
            if start + dur + pad > starttime and start - pad < endtime:
                return [dirname]
            else:
                return ret
        return ret

    first_four_start = int(starttime) / 100000
    first_four_end = int(endtime) / 100000

    for filename in os.listdir(dirname):
        m = matchfile.match(filename)
        if m:
            (start, dur) = (int(m.group(1)), int(m.group(2)))
            if start + dur + pad > starttime and start - pad < endtime:
                ret.append(os.path.join(dirname, filename))
        elif re.match('.*-[0-9]{5}$', filename):
            dirtime = int(filename[-5:])
            if dirtime >= first_four_start and dirtime \
                <= first_four_end:
                ret += get_all_files_in_range(os.path.join(dirname,
                        filename), starttime, endtime, pad=pad,
                        suffix=suffix)
        elif re.match('.*-[0-9]{4}$', filename):
            dirtime = int(filename[-4:])
            if dirtime >= first_four_start and dirtime \
                <= first_four_end:
                ret += get_all_files_in_range(os.path.join(dirname,
                        filename), starttime, endtime, pad=pad,
                        suffix=suffix)
        elif '.' not in filename and filename != 'latest':

            # FIX for dealing with dmt segments files. We need to avoid 'latest' which is a peropidically updated file.
            # It causes transient failures.
            # Keep recursing, we may be looking at directories of
            # ifos, each of which has directories with times

            ret += get_all_files_in_range(os.path.join(dirname,
                    filename), starttime, endtime, pad=pad,
                    suffix=suffix)

    return ret


def gather_trig_files(
    source_dir,
    gps_start_time,
    gps_end_time,
    glob_pattern='*',
    file_tag='',
    extension='',
    ):

    dirs = glob.glob(source_dir + '/' + glob_pattern)
    dirs.sort()
    dirs_in_segment = []
    dir_times = []
    for dir in dirs:
        dir_times.append(int(dir.split('-')[-1] + '00000'))
    dir_times = numpy.asarray(dir_times)
    first_time_in_segment = numpy.searchsorted(dir_times,
            gps_start_time, 'right')
    last_time_in_segment = numpy.searchsorted(dir_times, gps_end_time,
            'left')

    # expand by one directory on each end just to be safe

    dirs_in_segment = dirs[first_time_in_segment
        - 2:last_time_in_segment + 1]

    files = []
    for dir in dirs_in_segment:
        files.extend(glob.glob(dir + '/*' + file_tag + '*.'
                     + extension))

    files.sort()
    file_times = []
    for file in files:
        file_times.append(int(file.split('-')[-2]))
    file_times = numpy.asarray(file_times)

    first_time_in_segment = numpy.searchsorted(file_times,
            gps_start_time, 'right')
    last_time_in_segment = numpy.searchsorted(file_times, gps_end_time,
            'right')
    files_in_segment = files[first_time_in_segment
        - 1:last_time_in_segment]

    return files_in_segment


# find *.dat files within a certain time range. looks at just the files in the specified directory
def gather_datfiles(
    gpsstart,
    gpsstop,
    classifier=False,
    source_dir='./',
    ):
    """ 
....looks for *.dat files in source_dir that fall within GPSstart and GPSstop and are generated by a given classifier. Assumes the files have the following structure:
....  *-GPSstart-DUR.dat
....where DUR = GPSstop - GPSstart

...."""

# ....if classifiers:
# ........classifiers = list(classifiers)

    # find only those files that overlap [GPSstart, GPSstop]

    files = []
    for f in glob.glob(source_dir + '/*.dat'):  # glob.glob(source_dir + "/*.dat"): # iterate through all *.dat files

        if not classifier or classifier in f:
            name = f.strip('.dat').split('/')[-1].split('-')
            start = int(name[-2])
            stop = start + int(name[-1])
            if gpsstart <= start < gpsstop or gpsstart < stop \
                <= gpsstop:
                files.append(f)
            else:
                pass

# ....print gpsstart, gpsstop

    return files

#=================================================
# methods that convert datfiles into other data products
#=================================================
def datarray_to_roc(
    datarray,
    columns,
    filename,
    cluster_win=False,
    unsafe_win=False,
    gw_thr=False,
    ):
    """ 
....takes a 'datarray' and computes the appropriate standard roc file 
....includes steps like clustering and filtering cleans for times that are too close to glitches
....datarray is defined in datfiles_to_roc() (essentially a list of events, columns defines the indexies for values)
...."""

    if (cluster_win or unsafe_win) and ('GPS' not in columns or 'signif'
             not in columns):
        sys.exit('must supply a GPS column in order to cluster or determine unsafe cleans'
                 )
    if cluster_win and 'signif' not in columns:
        sys.exit('must supply signif column to cluster')
    if gw_thr and 'signif' not in columns:
        sys.exit('must supply signif column to threshold')

    roc = dict([(k, ind) for (ind, k) in enumerate(columns)])

    print 'No. events ' + str(len(datarray))

    # filter into cleans and glitches

    gw_gps = []
    cl_gps = []
    for l in datarray:
        if float(l[roc['i']]) == 1:
            gw_gps.append([float(line) for line in l])  #  gw events
        else:
            cl_gps.append([float(line) for line in l])

    # threshold

    if gw_thr:
        print '\tNo. gw_gps before thresholding ' + str(len(gw_gps))
        gw_gps = [gps for gps in gw_gps if gps[roc['signif']] >= gw_thr]
        print '\tNo. gw_gps after thresholding ' + str(len(gw_gps))

    # cluster

    if cluster_win:  # cluster gps times for glitches
        print '\tNo. gw_gps before clustering ' + str(len(gw_gps))
        gw_gps = cluster(gw_gps, roc, cluster_window=cluster_win)
        print '\tNo. gw_gps after clustering = ' + str(len(gw_gps))

    # remove cleans that are too close to glitches

    if unsafe_win:
        gw_seg = event.vetosegs(gw_gps, unsafe_win, tcent=roc['GPS'])  # define the "unclean" segments
        print '\tNo. cl_gps before windowing ' + str(len(cl_gps))
        cl_gps = event.exclude(cl_gps, gw_seg, tcent=roc['GPS'])  # only keep the cleans that fall outside of gw_seg
        print '\tNo. cl_gps after windowing ' + str(len(cl_gps))

    # compute cumulative statistics and write file
# ....f = open(filename, "w")

    datarray = gw_gps + cl_gps
    datarray.sort(key=lambda line: line[roc['rank']], reverse=True)  # sort roc list by ranks (High to Low)
    print 'No. events ' + str(len(datarray))

    if len(datarray) == 0:  # no data

# ........print >>f, "0\n0\n0 0 0"
# ........f.close()

        from laldetchar.idq.idq_summary_plots import rcg_to_ROC
        return rcg_to_ROC(
            filename,
            [0],
            [0],
            [0],
            0,
            0,
            )

# ........return filename

    i_vec = [int(l[roc['i']]) for l in datarray]
    r_vec = [float(l[roc['rank']]) for l in datarray]
    tot_n_g = sum(i_vec)  # total number of glitches
    tot_n_c = len(i_vec) - tot_n_g  # total number of clean samples

# ....print >>f, tot_n_g
# ....print >>f, tot_n_c

    rank = r_vec[0]
    c_n_g = 0  # cumulative number of glitches
    c_n_c = 0  # cumulative number of cleans
    c_cln = []
    c_gch = []
    c_ranks = []
    for (ind, i) in enumerate(i_vec):
        if rank != r_vec[ind]:
            c_ranks.append(rank)
            c_cln.append(c_n_c)
            c_gch.append(c_n_g)

# ............print >>f, rank, c_n_c, c_n_g

            rank = r_vec[ind]
        c_n_g += i
        c_n_c += 1 - i

    c_ranks.append(rank)
    c_cln.append(c_n_c)
    c_gch.append(c_n_g)

# ....print >>f, rank, c_n_c, c_n_g
# ....f.close()

    from laldetchar.idq.idq_summary_plots import rcg_to_ROC
    return rcg_to_ROC(
        filename,
        c_ranks,
        c_cln,
        c_gch,
        tot_n_c,
        tot_n_g,
        )


# ....return filename

# gather short-stride *.dat output files and write a combined *.roc file

def datfiles_to_roc(
    gpsstart,
    gpsstop,
    columns=False,
    classifiers=False,
    basename=False,
    source_dir='./',
    output_dir='./',
    cluster_win=False,
    unsafe_win=False,
    gw_thr=False,
    switch_snr_signif=[],
    ):
    """ 
....looks for all *.dat files in a directory structure and combines their contents into a single ROC file. 
....assumes the directory structue is as follows:
....  source_dir / *-time / *.dat
....where time = GPSstart % 1e5
....if classifiers:
....  only returns the *.dat files that have the structre *-classifier-gpsstart-dur.dat if (classifier in classifiers)
....  where dur = gpsstop - gpsstart
...."""

    if not columns:
        columns = ['rank', 'i']
    elif 'rank' not in columns or 'i' not in columns:
        sys.exit("columns must contain  fields 'rank' and 'i', exiting...."
                 )
    elif switch_snr_signif and ('SNR' not in columns or 'signif'
                                not in columns):
        sys.exit("columns must contain fields 'SNR' and 'signif' to switch snr and signif"
                 )

    roc_paths = []

    # if not classifier is supplied, accept all classifiers

    if not classifiers:
        classifiers = ['']

    # iterate for each classifier

    for classifier in classifiers:
        print classifier

        # find all appropriate subdirectories of source_dir and gather datfiles

        datfiles = []
        for dir in os.listdir(source_dir):
            try:
                d = dir.strip('/').split('-')
                start = int(d[-1])  # last element of d must be castable as an integer
                if start * 1e5 <= gpsstart < (start + 1) * 1e5 or start \
                    * 1e5 < gpsstop <= (start + 1) * 1e5:
                    if not basename or basename in dir:  # check basename
                        datfiles += gather_datfiles(gpsstart, gpsstop,
                                classifier=classifier,
                                source_dir=source_dir + '/' + dir)
            except:
                pass

        # gather roc information as tuples: (rank, i)
        # i=1 --> glitch
        # i=0 --> clean

        roc = dict([(k, ind) for (ind, k) in enumerate(columns)])
        dat = []
        livetime = 0
        for f in datfiles:
            try:
                output = slim_load_datfile(f, skip_lines=0,
                        columns=columns)
            except:
                print 'ERROR loading ' + f
            else:
                dat += [[output[k][ind] for k in columns] for ind in
                        range(len(output['rank']))]

                # get live time from the file name and add it to the total livetime.

                livetime += float(f.split('-')[-1].split('.')[0])

        # switch snr and signif

        if classifier in switch_snr_signif:
            print "  --> switching 'SNR' and 'signif' column labels for " \
                + classifier
            snr_ind = columns.index('SNR')
            signif_ind = columns.index('signif')
            columns[snr_ind] = 'signif'
            columns[signif_ind] = 'SNR'

        # compute cumulative statistics and write file

        roc_path = output_dir + '/'
        if basename:
            roc_path += basename + '-'
        if classifier != '':
            roc_path += classifier + '-'
        roc_path += repr(gpsstart) + '-' + repr(gpsstop - gpsstart) \
            + '.roc'

        roc_paths.append((datarray_to_roc(
            dat,
            columns,
            filename=roc_path,
            cluster_win=cluster_win,
            unsafe_win=unsafe_win,
            gw_thr=gw_thr,
            ), classifier))

        # read number of glitches and cleans for roc file

        roc_file = open(roc_path, 'r')
        lines = roc_file.readlines()
        n_glitches = int(lines[0].strip())
        n_cleans = int(lines[1].strip())
        roc_file.close()

        # append livetime, total number of glitches and cleans to the stat.txt file

        stat_path = output_dir + '/'
        if basename:
            stat_path += basename + '_'
        stat_path += 'stat-'
        stat_path += repr(gpsstart) + '-' + repr(gpsstop - gpsstart) \
            + '.stat'

        # check if summary stat file already exist
        # if it does append to it, if not creat one and add period analyzed on the first line

        if os.path.exists(stat_path):
            stat_file = open(stat_path, 'a')
        else:
            stat_file = open(stat_path, 'w')
            stat_file.write('Period ' + repr(gpsstart) + '-'
                            + repr(gpsstop - gpsstart) + '\n')
        stat_file.write(classifier + '\n')
        stat_file.write('Livetime = ' + str(livetime) + '\n')
        stat_file.write('Glitches = ' + str(n_glitches) + '\n')
        stat_file.write('Cleans = ' + str(n_cleans) + '\n')
        stat_file.close()
    return roc_paths


def datfilename_to_xmldocs(
    datfilename,
    ifo,
    FAPmap,
    Lmap=False,
    Effmap=False,
    classifier=None,
    ):
    """
....reads in a datfile and generates a idq_tables.GlitchTable containing the equivalent information, but adds FAP and likelihood as computed by FAPmap and Lmap (expected to be callable)
........if Lmap is not supplied, attempts to use L = Eff/FAP as estimate with Effmap
........==> must supply either Lmap or Effmap

....return gchTable, clnTable
...."""

    if not (Lmap or Effmap):
        raise StandardError('must supply either Lmap or Effmap in idq.datfilename_to_xmldocs()'
                            )

    gchxml_doc = ligolw.Document()  # generate the xml document object
    gch_ligolw_element = ligolw.LIGO_LW()  # generate the table-tree element
    gchxml_doc.appendChild(gch_ligolw_element)  # put table-tree into xml document

    clnxml_doc = ligolw.Document()  # repeat for clnTable
    cln_ligolw_element = ligolw.LIGO_LW()
    clnxml_doc.appendChild(cln_ligolw_element)

    gchTable = lsctables.New(idq_tables.GlitchTable)  # instantiate table objects
    clnTable = lsctables.New(idq_tables.GlitchTable)

    gch_ligolw_element.appendChild(gchTable)  # put gchTable into table-tree
    cln_ligolw_element.appendChild(clnTable)

    # iterate over dat and place info into appropriate tables

    gch_event_id = ilwd.ilwdchar('glicth:event_id:0')
    cln_event_id = ilwd.ilwdchar('glitch:event_id:0')

        # load data and fill out xml tables

    if classifier == 'ovl':  # ovl records more info than the other classifiers

        # setup ovl specific tables

        gchOVLTable = lsctables.New(idq_tables.OVLDataTable)  # instantiate table objects
        clnOVLTable = lsctables.New(idq_tables.OVLDataTable)
        gch_ligolw_element.appendChild(gchOVLTable)  # add tables to document tree
        cln_ligolw_element.appendChild(clnOVLTable)
        gch_ovl_id = ilwd.ilwdchar('ovl_data:event_id:0')
        cln_ovl_id = ilwd.ilwdchar('ovl_data:event_id:0')

        # read in data

        dat = slim_load_datfile(datfilename, skip_lines=0, columns=[  # read in relevant information from datfile
            'GPS',
            'i',
            'rank',
            'vchan',
            'vthr',
            'vwin',
            ])
    else:
        dat = slim_load_datfile(datfilename, skip_lines=0,
                                columns=['GPS', 'i', 'rank'])  # read in relevant information from datfile

    no_events = len(dat['GPS'])
    for ind in range(no_events):
        GPS = float(dat['GPS'][ind])
        i = float(dat['i'][ind])
        rank = float(dat['rank'][ind])

        # ## add info to GlitchTable() objects. Commmon for all classifiers

        row = idq_tables.Glitch()  # define row object

        row.ifo = ifo
        row.gps = int(GPS)  # fill in relevant data
        row.gps_ns = (GPS - row.gps) * 1e9
        row.rank = float(rank)
        row.fap = FAPmap(row.rank)
        if Lmap:
            row.likelihood = Lmap(row.rank)
        else:
            if row.fap == 0:
                if Effmap(row.rank) == 0:
                    row.likelihood = 0
                else:
                    row.likelihood = numpy.infty
            else:
                row.likelihood = Effmap(row.rank) / row.fap

                # ## add info to OVLDataTable() objects

        if classifier == 'ovl':
            ovl_row = idq_tables.OVLData()  # define row object
            ovl_row.ifo = ifo
            ovl_row.aux_channel = dat['vchan'][ind]
            ovl_row.veto_thr = float(dat['vthr'][ind])
            ovl_row.veto_win = float(dat['vwin'][ind])

        # ## update rows, etc. Do this all at once so the index counters don't get confused

        if i:  # glitch sample
            row.event_id = gch_event_id
            gchTable.append(row)
            if classifier == 'ovl':
                ovl_row.event_id = gch_ovl_id
                gchOVLTable.append(ovl_row)
                gch_ovl_id += 1
            gch_event_id += 1
        else:
            row.event_id = cln_event_id
            clnTable.append(row)
            if classifier == 'ovl':
                ovl_row.event_id = cln_ovl_id
                clnOVLTable.append(ovl_row)
                cln_ovl_id += 1
            cln_event_id += 1

    return (gchxml_doc, clnxml_doc)


##################################################

def datfiles_to_chanlist(
    gpsstart,
    gpsstop,
    columns=False,
    classifiers=['ovl'],
    basename=False,
    source_dir='./',
    output_dir='./',
    cluster_win=False,
    unsafe_win=False,
    gw_thr=False,
    switch_snr_signif=False,
    ):
    """
....computes the performance of each auxiliary channel based off datfile output.
....currently only supports ovl datfile output, but could be adapted to other classifiers assuming they have a column 'vchan' in their output.
...."""

    if not columns:
        columns = ['i', 'vchan']
    elif 'i' not in columns or 'vchan' not in columns:
        sys.exit("columns must contain fields 'i', 'vchan'. exiting...")
    elif cluster_win and 'signif' not in columns:
        sys.exit('must supply a signif column to cluster')
    elif gw_thr and 'signif' not in columns:
        sys.exit('must supply a signif column to cluster')
    elif switch_snr_signif and ('SNR' not in columns or 'signif'
                                not in columns):
        sys.exit("must supply 'SNR' and 'signif' columns to switch them"
                 )

    chanlist = []

    # pick up each classifier separately

    for classifier in classifiers:
        print classifier
        datfiles = []
        for dir in os.listdir(source_dir):
            try:
                d = dir.strip('/').split('-')
                start = int(d[-1])  # expect a particular directory structure
                if start * 1e5 <= gpsstart < (start + 1) * 1e5 or start \
                    * 1e5 < gpsstop <= (start + 1) * 1e5:
                    if not basename or basename in dir:
                        datfiles += gather_datfiles(gpsstart, gpsstop,
                                classifier=classifier,
                                source_dir=source_dir + '/' + dir)
            except:
                pass

        # gather chan performance measures

        dat = []
        for f in datfiles:
            try:
                output = slim_load_datfile(f, skip_lines=0,
                        columns=columns)
            except:
                print 'ERROR loading ' + f
            else:
                dat += [[output[k][ind] for k in columns] for ind in
                        range(len(output['vchan']))]

        # switch SNR and significance

        if switch_snr_signif:
            print '  --> switching SNR and signif column labels'
            snr_ind = columns.index('SNR')
            signif_ind = columns.index('signif')
            columns[snr_ind] = 'signif'
            columns[signif_ind] = 'SNR'

        columnD = dict([(c, ind) for (ind, c) in enumerate(columns)])

        # compute cumulative statistics and write file

        print 'No. events ' + str(len(dat))

        # filter into cleans and glitches

        gw_gps = []
        cl_gps = []
        for l in dat:
            L = []
            for c in columns:
                if c == 'vchan':
                    L.append(l[columnD['vchan']])
                elif c == 'i':
                    L.append(int(float(l[columnD['i']])))
                else:
                    L.append(float(l[columnD[c]]))
            if L[columnD['i']] == 1:
                gw_gps.append(L)  #  gw events
            else:
                cl_gps.append(L)

        # threshold

        if gw_thr:
            print '\tNo. gw_gps before thresholding ' + str(len(gw_gps))
            gw_gps = [gps for gps in gw_gps if gps[columnD['signif']]
                      >= gw_thr]
            print '\tNo. gw_gps after thresholding ' + str(len(gw_gps))

        # cluster

        if cluster_win:  # cluster gps times for glitches
            print '\tNo. gw_gps before clustering ' + str(len(gw_gps))
            gw_gps = cluster(gw_gps, columnD,
                             cluster_window=cluster_win)
            print '\tNo. gw_gps after clustering = ' + str(len(gw_gps))

        # remove cleans that are too close to glitches

        if unsafe_win:
            gw_seg = event.vetosegs(gw_gps, unsafe_win,
                                    tcent=columnD['GPS'])  # define the "unclean" segments
            print '\tNo. cl_gps before windowing ' + str(len(cl_gps))
            cl_gps = event.exclude(cl_gps, gw_seg, tcent=columnD['GPS'])  # only keep the cleans that fall outside of gw_seg
            print '\tNo. cl_gps after windowing ' + str(len(cl_gps))

        n_gch = float(len(gw_gps))
        n_cln = float(len(cl_gps))
        dat = gw_gps + cl_gps
        dat.sort(key=lambda line: line[columnD['GPS']], reverse=True)
        print 'No. events ' + str(len(dat))

        # associate glitches with a channel, configuration, etc

        chan_perform = {}
        for d in dat:

            # associate events with channel name

            chan = d[columnD['vchan']]
            if chan not in chan_perform.keys():
                chan_perform[chan] = {'c_cln': 0, 'c_gch': 0}
            if int(d[columnD['i']]) == 1:
                chan_perform[chan]['c_gch'] += 1
            else:
                chan_perform[chan]['c_cln'] += 1

        # write channel list file

        path = output_dir + '/'
        if basename:
            path += basename + '-'
        if classifier != '':
            path += classifier + '-'
        path += repr(gpsstart) + '-' + repr(gpsstop - gpsstart) \
            + '_channel_summary.txt'

        f = open(path, 'w')
        print >> f, '%-44s %7s %7s %9s %9s %12s' % (
            'channel',
            'n_cln',
            'n_gch',
            'c_fap',
            'c_eff',
            'eff/fap',
            )

        # compute statistics for each channel

        for chan in sorted(chan_perform.keys()):
            this_chanD = chan_perform[chan]
            c_cln = this_chanD['c_cln']
            c_gch = this_chanD['c_gch']
            if n_cln != 0:
                c_fap = c_cln / n_cln
            else:
                c_fap = 0
            if n_gch != 0:
                c_eff = c_gch / n_gch
            else:
                c_eff = 0
            chan_perform[chan]['c_fap'] = c_fap
            chan_perform[chan]['c_eff'] = c_eff
            if c_fap != 0:
                eff_fap = c_eff / c_fap
            elif c_eff == 0:
                eff_fap = 0
            else:
                eff_fap = numpy.inf

            # write statisitics to file

            print >> f, '%-44s %7.1f %7.1f %8.7f %8.7f %12.5f' % (
                chan,
                c_cln,
                c_gch,
                c_fap,
                c_eff,
                eff_fap,
                )

        f.close()
        chanlist.append((chan_perform, classifier, path))

    return chanlist


##################################################

def datfiles_to_ranklist(
    gpsstart,
    gpsstop,
    columns=False,
    classifiers=['ovl'],
    basename=False,
    source_dir='./',
    output_dir='./',
    cluster_win=False,
    unsafe_win=False,
    gw_thr=False,
    switch_snr_signif=False,
    ):
    """
....generates a list that associates gps times with each rank for all datfiles.
...."""

    if not columns:
        columns = ['i', 'GPS', 'rank']
    elif 'GPS' not in columns or 'i' not in columns or 'rank' \
        not in columns:
        sys.exit("columns must contain fields 'GPS', 'i', and 'rank'. exiting..."
                 )
    elif cluster_win and 'signif' not in columns:
        sys.exit('must supply a signif column to cluster')
    elif gw_thr and 'signif' not in columns:
        sys.exit('must supply a signif column to cluster')
    elif switch_snr_signif and ('SNR' not in columns or 'signif'
                                not in columns):
        sys.exit("must supply 'SNR' and 'signif' columns to switch them"
                 )

    ranklist = []

    # pick up each classifier separately

    for classifier in classifiers:
        print classifier
        datfiles = []
        for dir in os.listdir(source_dir):
            try:
                d = dir.strip('/').split('-')
                start = int(d[-1])  # expect a particular directory structure
                if start * 1e5 <= gpsstart < (start + 1) * 1e5 or start \
                    * 1e5 < gpsstop <= (start + 1) * 1e5:
                    if not basename or basename in dir:
                        datfiles += gather_datfiles(gpsstart, gpsstop,
                                classifier=classifier,
                                source_dir=source_dir + '/' + dir)
            except:
                pass

        # gather chan performance measures

        dat = []
        for f in datfiles:
            try:
                output = slim_load_datfile(f, skip_lines=0,
                        columns=columns)
            except:
                print 'ERROR loading ' + f
            else:
                dat += [[output[k][ind] for k in columns] for ind in
                        range(len(output['rank']))]

        # switch SNR and significance

        if switch_snr_signif:
            print '  --> switching SNR and signif column labels'
            snr_ind = columns.index('SNR')
            signif_ind = columns.index('signif')
            columns[snr_ind] = 'signif'
            columns[signif_ind] = 'SNR'

        columnD = dict([(c, ind) for (ind, c) in enumerate(columns)])

        # compute cumulative statistics and write file

        print 'No. events ' + str(len(dat))

        # filter into cleans and glitches

        gw_gps = []
        cl_gps = []
        for l in dat:
            L = []
            for c in columns:
                if c == 'vchan':
                    L.append(l[columnD['vchan']])
                elif c == 'i':
                    L.append(int(float(l[columnD['i']])))
                else:
                    L.append(float(l[columnD[c]]))
            if L[columnD['i']] == 1:
                gw_gps.append(L)  #  gw events
            else:
                cl_gps.append(L)

        # threshold

        if gw_thr:
            print '\tNo. gw_gps before thresholding ' + str(len(gw_gps))
            gw_gps = [gps for gps in gw_gps if gps[columnD['signif']]
                      >= gw_thr]
            print '\tNo. gw_gps after thresholding ' + str(len(gw_gps))

        # cluster

        if cluster_win:  # cluster gps times for glitches
            print '\tNo. gw_gps before clustering ' + str(len(gw_gps))
            gw_gps = cluster(gw_gps, columnD,
                             cluster_window=cluster_win)
            print '\tNo. gw_gps after clustering = ' + str(len(gw_gps))

        # remove cleans that are too close to glitches

        if unsafe_win:
            gw_seg = event.vetosegs(gw_gps, unsafe_win,
                                    tcent=columnD['GPS'])  # define the "unclean" segments
            print '\tNo. cl_gps before windowing ' + str(len(cl_gps))
            cl_gps = event.exclude(cl_gps, gw_seg, tcent=columnD['GPS'])  # only keep the cleans that fall outside of gw_seg
            print '\tNo. cl_gps after windowing ' + str(len(cl_gps))

        n_gch = float(len(gw_gps))
        n_cln = float(len(cl_gps))
        dat = gw_gps + cl_gps
        dat.sort(key=lambda line: line[columnD['GPS']], reverse=True)
        print 'No. events ' + str(len(dat))

        # associate glitches with a channel, configuration, etc

        rank_gps = {}
        for d in dat:

            # associate events with rank

            rank = d[columnD['rank']]
            if rank not in rank_gps.keys():
                rank_gps[rank] = {'cln': [], 'gch': []}
            if int(d[columnD['i']]) == 1:
                rank_gps[rank]['gch'].append(d[columnD['GPS']])
            else:
                rank_gps[rank]['cln'].append(d[columnD['GPS']])

        # write channel list file

        path = output_dir + '/'
        if basename:
            path += basename + '-'
        if classifier != '':
            path += classifier + '-'
        path += repr(gpsstart) + '-' + repr(gpsstop - gpsstart) \
            + '_rank_summary.txt'

        f = open(path, 'w')
        print >> f, '''# rank
# gch_gps
# cln_gps'''
        for rank in sorted(rank_gps.keys(), reverse=True):
            gch_gps = ''
            for g in rank_gps[rank]['gch']:
                gch_gps += '%f ' % g
            cln_gps = ''
            for c in rank_gps[rank]['cln']:
                cln_gps += '%f ' % c
            print >> f, '''%10.9f
%s
%s
#''' % (rank, gch_gps, cln_gps)

        f.close()
        ranklist.append((rank_gps, classifier, path))

        return ranklist


##################################################

def datfiles_to_configlist(
    gpsstart,
    gpsstop,
    columns=False,
    classifiers=['ovl'],
    basename=False,
    source_dir='./',
    output_dir='./',
    cluster_win=False,
    unsafe_win=False,
    gw_thr=False,
    switch_snr_signif=False,
    ):
    """
....generates an output structure similar to ovl vetolists using (channel name, threshold, window) and associates with it statistics and a rank.
...."""

    if not columns:
        columns = ['i', 'vchan', 'vthr', 'vwin', 'rank']
    elif 'i' not in columns or 'vchan' not in columns or 'vthr' \
        not in columns or 'vwin' not in columns or 'rank' \
        not in columns:
        sys.exit("columns must contain fields 'i', 'vchan', 'vthr', 'vwin', 'rank'. exiting..."
                 )
    elif cluster_win and 'signif' not in columns:
        sys.exit('must supply a signif column to cluster')
    elif gw_thr and 'signif' not in columns:
        sys.exit('must supply a signif column to cluster')
    elif switch_snr_signif and ('SNR' not in columns or 'signif'
                                not in columns):
        sys.exit("must supply 'SNR' and 'signif' columns to switch them"
                 )

    configlist = []

    # pick up each classifier separately

    for classifier in classifiers:
        print classifier
        datfiles = []
        for dir in os.listdir(source_dir):
            try:
                d = dir.strip('/').split('-')
                start = int(d[-1])  # expect a particular directory structure
                if start * 1e5 <= gpsstart < (start + 1) * 1e5 or start \
                    * 1e5 < gpsstop <= (start + 1) * 1e5:
                    if not basename or basename in dir:
                        datfiles += gather_datfiles(gpsstart, gpsstop,
                                classifier=classifier,
                                source_dir=source_dir + '/' + dir)
            except:
                pass

        # gather chan performance measures

        dat = []
        for f in datfiles:
            try:
                output = slim_load_datfile(f, skip_lines=0,
                        columns=columns)
            except:
                print 'ERROR loading ' + f
            else:
                dat += [[output[k][ind] for k in columns] for ind in
                        range(len(output['vchan']))]

        # switch SNR and significance

        if switch_snr_signif:
            print '  --> switching SNR and signif column labels'
            snr_ind = columns.index('SNR')
            signif_ind = columns.index('signif')
            columns[snr_ind] = 'signif'
            columns[signif_ind] = 'SNR'

        columnD = dict([(c, ind) for (ind, c) in enumerate(columns)])

        # compute cumulative statistics and write file

        print 'No. events ' + str(len(dat))

        # filter into cleans and glitches

        gw_gps = []
        cl_gps = []
        for l in dat:
            L = []
            for c in columns:
                if c == 'vchan':
                    L.append(l[columnD['vchan']])
                elif c == 'i':
                    L.append(int(float(l[columnD['i']])))
                else:
                    L.append(float(l[columnD[c]]))
            if L[columnD['i']] == 1:
                gw_gps.append(L)  #  gw events
            else:
                cl_gps.append(L)

        # threshold

        if gw_thr:
            print '\tNo. gw_gps before thresholding ' + str(len(gw_gps))
            gw_gps = [gps for gps in gw_gps if gps[columnD['signif']]
                      >= gw_thr]
            print '\tNo. gw_gps after thresholding ' + str(len(gw_gps))

        # cluster

        if cluster_win:  # cluster gps times for glitches
            print '\tNo. gw_gps before clustering ' + str(len(gw_gps))
            gw_gps = cluster(gw_gps, columnD,
                             cluster_window=cluster_win)
            print '\tNo. gw_gps after clustering = ' + str(len(gw_gps))

        # remove cleans that are too close to glitches

        if unsafe_win:
            gw_seg = event.vetosegs(gw_gps, unsafe_win,
                                    tcent=columnD['GPS'])  # define the "unclean" segments
            print '\tNo. cl_gps before windowing ' + str(len(cl_gps))
            cl_gps = event.exclude(cl_gps, gw_seg, tcent=columnD['GPS'])  # only keep the cleans that fall outside of gw_seg
            print '\tNo. cl_gps after windowing ' + str(len(cl_gps))

        n_gch = float(len(gw_gps))
        n_cln = float(len(cl_gps))
        dat = gw_gps + cl_gps
        dat.sort(key=lambda line: line[columnD['GPS']], reverse=True)
        print 'No. events ' + str(len(dat))

        # associate glitches with a channel, configuration, etc

        config_perform = {}
        for d in dat:

            # associate events with channel name

            chan = d[columnD['vchan']]
            vthr = d[columnD['vthr']]
            vwin = d[columnD['vwin']]
            rank = d[columnD['rank']]
            config = (chan, vthr, vwin, rank)
            if config not in config_perform.keys():
                config_perform[config] = {'c_cln': 0, 'c_gch': 0}
            if int(d[columnD['i']]) == 1:
                config_perform[config]['c_gch'] += 1
            else:
                config_perform[config]['c_cln'] += 1

        # write channel list file

        path = output_dir + '/'
        if basename:
            path += basename + '-'
        if classifier != '':
            path += classifier + '-'
        path += repr(gpsstart) + '-' + repr(gpsstop - gpsstart) \
            + '_config_summary.txt'

        f = open(path, 'w')
        print >> f, '%-44s %5s %5s %11s %7s %7s %9s %9s %12s' % (
            'channel',
            'vthr',
            'vwin',
            'rank',
            'n_cln',
            'n_gch',
            'c_fap',
            'c_eff',
            'eff/fap',
            )

        # compute statistics for each channel

        sorted_keys = config_perform.keys()
        sorted_keys.sort(key=lambda line: line[3], reverse=True)
        for config in sorted_keys:
            (chan, vthr, vwin, rank) = config
            this_configD = config_perform[config]
            c_cln = this_configD['c_cln']
            c_gch = this_configD['c_gch']
            if n_cln != 0:
                c_fap = c_cln / n_cln
            else:
                c_fap = 0
            if n_gch != 0:
                c_eff = c_gch / n_gch
            else:
                c_eff = 0
            config_perform[config]['c_fap'] = c_fap
            config_perform[config]['c_eff'] = c_eff
            if c_fap != 0:
                eff_fap = c_eff / c_fap
            elif c_eff == 0:
                eff_fap = 0
            else:
                eff_fap = numpy.inf

            # write statisitics to file

            print >> f, \
                '%-44s %5.1f %5.3f %10.9f %7.1f %7.1f %8.7f %8.7f %12.5f' \
                % (
                chan,
                vthr,
                vwin,
                rank,
                c_cln,
                c_gch,
                c_fap,
                c_eff,
                eff_fap,
                )

        f.close()
        configlist.append((config_perform, classifier, path))

    return configlist


### convert a list of ranks into a set of FAP estimates

def rank_to_FAPmap(roc_filename, kind='linear'):
    """
    generates mapping for FAP estimates using the data stored in 'roc_filename'
        expect 'roc_filename' to be a string corresponding to the location of an roc file
    """

    from laldetchar.idq.idq_summary_plots import ROC_to_rcg
    from scipy.interpolate import interp1d

    (rs, c_cln, _, tot_cln, _) = ROC_to_rcg(roc_filename)

    # ensure interpolation range covers allowable range for rank

    if rs[0] != 1.0:
        rs.insert(0, 1.0)
        c_cln.insert(0, 0.0)
    if rs[-1] != 0.0:
        rs.append(0.0)
        c_cln.append(tot_cln)

    if tot_cln == 0:
        return interp1d(rs[::-1], numpy.ones((len(rs), )), kind=kind)
    else:
        return interp1d(rs[::-1], 1.0 * numpy.array(c_cln)[::-1]
                        / tot_cln, kind=kind)


def rank_to_Effmap(roc_filename, kind='linear'):
    """
    generates mapping for Efficiency estimates using the data stored in 'roc_filename'
        expect 'roc_filename' to be a string corresponding to the location of an roc file
    """

    from laldetchar.idq.idq_summary_plots import ROC_to_rcg
    from scipy.interpolate import interp1d

    (rs, _, c_gch, _, tot_gch) = ROC_to_rcg(roc_filename)

    # ensure interpolation range covers allowable range for rank

    if rs[0] != 1.0:
        rs.insert(0, 1.0)
        c_gch.insert(0, 0.0)
    if rs[-1] != 0.0:
        rs.append(0.0)
        c_gch.append(tot_gch)

    if tot_gch == 0:
        return interp1d(rs[::-1], numpy.zeros((len(rs), )), kind=kind)
    else:
        return interp1d(rs[::-1], 1.0 * numpy.array(c_gch)[::-1]
                        / tot_gch, kind=kind)


def rank_to_EffFAPmap(roc_filename, kind='linear'):
    """
    generates both FAP and Eff maps from roc_filename
      return Effmap, FAPmap
    """

    from laldetchar.idq.idq_summary_plots import ROC_to_rcg
    from scipy.interpolate import interp1d

    (rs, c_cln, c_gch, tot_cln, tot_gch) = ROC_to_rcg(roc_filename)

    # ensure interpolation range covers allowable range for rank

    if rs[0] != 1.0:
        rs.insert(0, 1.0)
        c_gch.insert(0, 0.0)
        c_cln.insert(0, 0.0)
    if rs[-1] != 0.0:
        rs.append(0.0)
        c_gch.append(tot_gch)
        c_cln.append(tot_cln)

    if tot_gch == 0:
        Effmap = interp1d(rs[::-1], numpy.zeros((len(rs), )), kind=kind)
    else:
        Effmap = interp1d(rs[::-1], 1.0 * numpy.array(c_gch)[::-1]
                          / tot_gch, kind=kind)

    if tot_cln == 0:
        FAPmap = interp1d(rs[::-1], numpy.ones((len(rs), )), kind=kind)
    else:
        FAPmap = interp1d(rs[::-1], 1.0 * numpy.array(c_cln)[::-1]
                          / tot_cln, kind=kind)

    return (Effmap, FAPmap)


############################ Single channel KW trigger collection functions ######################################################

### collect and generate single-channel KW trig files

def collect_sngl_chan_kw(
    gpsstart,
    gpsstop,
    kw_config,
    width=False,
    source_dir='./',
    output_dir='./',
    verbose=False,
    ):
    """ 
....replicates the action of collect.py in that it moves KW triggers from multi-channel (low-latency) files to single channel files. These files will be written in to output_dir in a subdirecitory titled `gpsstart`_`gpsstop`
...."""

    # find the start time of the day which includes start

    if width:
        t0 = 847497600  # 00:00 GPS from a reference day?
        t = t0 + int((gpsstart - t0) / width) * width  # beginning time from reference time with strides of width
    else:
        t = gpsstart
        width = gpsstop - gpsstart

    # wait for previous jobs to finish

    if t + 128 > nowgps():  # wait 128sec after the start of a day before processing jobs
        if verbose:
            print 'sleeping ' + `wait` \
                + ' seconds to allow jobs to complete...'
        time.sleep(128)

    gps_start = []

    # add periods to process until we run into the end time
    # if a directory does not exist, we create it.

    if verbose:
        print t  # start of first day to be processed
    while t + width <= gpsstop:
        gps_start.append(t)
        range = repr(t) + '_' + repr(t + width)
        if verbose:
            print 'adding period to process: ' + range
        if not os.path.isdir(output_dir + '/' + range):
            os.makedirs(output_dir + '/' + range)
        t += width

    # process channels (listed in confg file)
    # generate a list of channels in the config file

    if verbose:
        print 'processing configuration file ' + kw_config
    stride = False
    basename = False
    channels = []
    f = open(kw_config, 'r')
    for line in f:
        fields = line.strip().split()
        if fields[0] == 'stride':
            stride = float(fields[1])
        elif fields[0] == 'basename':
            basename = fields[1]
        elif fields[0] == 'channel':
            channels.append([fields[1], (fields[1])[:2] + '_'
                            + (fields[1])[3:] + '_' + fields[2] + '_'
                            + fields[3]])

    if not stride:
        print 'stride not defined in ' + config.kwconfig
        sys.exit(1)
    if not basename:
        print 'basename not defined in ' + config.kwconfig
        sys.exit(1)
    f.close()

    if verbose:
        print ' found ' + repr(len(channels)) + ' channels'

    # generate single-channel summary files

    new_dirs = []  # record which new directories are created (for training jobs)

    # iterate through new data sets

    for day_start in gps_start:
        triggers = dict()  # dictionary holding triggers
        for (channel, tag) in channels:
            triggers[tag] = []
        day_end = day_start + width

        # in case stride does not fit neatly into GPS day

        t = day_start - day_start % stride  # beginning time for this day
        range = repr(day_start) + '_' + repr(day_end)
        stride_S = int(stride)
        lastt = 0

        new_dirs.append(output_dir + '/' + range)  # add new directory to list for training jobs (created in output_dir)

        # gather all the KW triggers and sort them by channel

        segments = []  # segment list for this day
        while t < day_end:
            t_S = int(t + 0.0000001)  # make sure we don'ti have a rounding error?
            t_dir = t_S / 100000  # the digits in the GPS time beyond 1e5, used for structure of KW output
            file = source_dir + '/' + basename + '-' + repr(t_dir) \
                + '/' + basename + '-' + repr(t_S) + '-' \
                + repr(stride_S) + '.trg'

            # build segment list for this day

            if os.path.exists(file):
                if len(segments) == 0:  # first segment
                    segments.append([t, t + stride])
                    lastt = t
                elif segments[-1][1] == t:

                                           # continuous data segment already generated

                    segments[-1][1] = t + stride
                    lastt = t
                else:

                      # discontinuous data, we skip one section of data to eschew filter transients

                    if lastt + stride == t:  # not continuous with anything in segments, but we've already skipped one section of data
                        if verbose:
                            print 'forming new segment at ' + `t`
                        segments.append([t, t + stride])
                    else:

                          # skip the first section of data in a new segment because of 'filter transients'

                        lastt = t
                        t += stride
                        continue

                # read in triggers from file and sort by channel

                if verbose:
                    print ' -> reading ' + file
                try:
                    f = open(file, 'r')

                    # for stride which overlaps day boundary, will get some extra triggers
                    # not worth extra processing to parse GPS times to sort them out

                    for line in f:
                        fields = line.strip().split()
                        if len(fields) == 9 and fields[-1] != '':
                            tag = fields[-1]
                            if not triggers.has_key(tag):
                                if verbose:
                                    print ' -> WARNING: triggers ' \
    + tag + ' not in configuration file (t=' + repr(t) \
    + '). Adding it to the list of channels.'
                                triggers[tag] = []
                            triggers[tag].append(line.strip()[:-len(tag)
                                    - 2])
                        else:
                            if verbose:
                                print '  --> ERROR: problem with line in trigger file ' \
                                    + file
                    f.close()
                except IOError:
                    print ' -> ERROR: problems opening file ' + file
            else:

                if verbose:
                    print ' -> WARNING: missing ' + file

            t += stride

        # write the single-channel summary files

        for tag in triggers.keys():
            file = output_dir + '/' + range + '/' + tag + '.trg'
            if verbose:
                print 'writing out triggers ' + file

            if len(triggers[tag]) == 0:
                if verbose:
                    print ' -> WARNING: no triggers found for ' + tag

            try:
                f = open(file, 'w')
                for line in triggers[tag]:
                    print >> f, line
                f.close()
            except:
                print ' -> ERROR: problem writing out trigger file ' \
                    + file

        # write segment file

        file = output_dir + '/' + range + '/' + basename + '.seg'
        f = open(file, 'w')
        for segment in segments:
            print >> f, '%.0f %.0f' % tuple(segment)
        f.close()

    return new_dirs


############  Generic functions for AuxMVC algorithms  ###########################

def build_auxmvc_vectors(
    trigger_dict,
    main_channel,
    time_window,
    signif_threshold,
    output_file_name,
    gps_start_time,
    gps_end_time,
    channels=None,
    unsafe_channels=None,
    science_segments=None,
    clean_times=None,
    clean_samples_rate=None,
    clean_window=None,
    filter_out_unclean=False,
    max_clean_samples=None,
    max_glitch_samples=None,
    ):
    """
....Given dictionary of triggers from multiple channels, the function constructs auxmvc
....vectors out of them. Result is saved into output_file_name. 
...."""

    if trigger_dict:
        print 'Number of triggers in the main channel before thresholding:', \
            len(trigger_dict[main_channel])
    else:

        # empty trig-dictionary, raise an error

        raise StandardError('Empty trig-dictionary. Can not build auxmvc vectors.'
                            )

    # use only channels from the channels file, if provided

    if channels:
        selected_channels = event.read_channels_from_file(channels)
        trigger_dict.keep_channels(selected_channels)

        # to ensure consistency in dimensionality of auxmvc vectors
        # add the channels from the selected channels that are absent in trigger_dict

        for channel in selected_channels:
            if not channel in trigger_dict.channels():
                trigger_dict[channel] = []

    # get rid of unsafe channels if provided

    if unsafe_channels:
        unsafe_channels = event.read_channels_from_file(unsafe_channels)
        trigger_dict.remove_channels(unsafe_channels)

    # keep only the triggers from the [gps_start_time, gps_end_time] segment

    # first keep all triggers from the segment expanded by the time concidence window, so that not to loose coincidences

    trigger_dict.include([[gps_start_time - time_window, gps_end_time
                         + time_window]])

    # then in the main channel keep triggers that fall within the segment.

    trigger_dict.include([[gps_start_time, gps_end_time]],
                         channels=[main_channel])

    # keep only triggers from the science segments if given........

    if science_segments:
        science_segments = event.andsegments([[gps_start_time
                - time_window, gps_end_time + time_window]],
                science_segments)
        trigger_dict.include(science_segments)

    # apply significance threshold to the triggers from the main channel

    trigger_dict.apply_signif_threshold(channels=[main_channel],
            threshold=signif_threshold)
    print 'Number of triggers in the main channel after thresholding:', \
        len(trigger_dict[main_channel])

    # construct glitch auxmvc vectors

    aux_glitch_vecs = event.build_auxmvc_vectors(trigger_dict,
            main_channel=main_channel,
            coincidence_time_window=time_window)

    # apply upper limit on the number of glitch samples if given

    if max_glitch_samples:
        if len(aux_glitch_vecs) > max_glitch_samples:
            aux_glitch_vecs = aux_glitch_vecs[-max_glitch_samples:]

    if not clean_times:
        if clean_samples_rate:

            # generate random times for clean samples

            if science_segments:
                clean_times = event.randomrate(clean_samples_rate,
                        event.andsegments([[gps_start_time
                        + time_window, gps_end_time - time_window]],
                        science_segments))
            else:
                clean_times = event.randomrate(clean_samples_rate,
                        [gps_start_time + time_window, gps_end_time
                        - time_window])
        else:
            clean_times = []

    # construct clean auxmvc vectors

    aux_clean_vecs = event.build_auxmvc_vectors(
        trigger_dict,
        main_channel=main_channel,
        coincidence_time_window=time_window,
        build_clean_samples=True,
        clean_times=clean_times,
        clean_window=clean_window,
        )

    # get rid of clean samples that are near real triggers in the main channel.

    if filter_out_unclean:
        aux_clean_vecs = auxmvc_utils.get_clean_samples(aux_clean_vecs)

    # apply upper limit on the number of clean samples if given

    if max_clean_samples:
        if len(aux_clean_vecs) > max_clean_samples:
            aux_clean_vecs = aux_clean_vecs[-max_clean_samples:]

    # convert glitch and clean auxmvc vectors into MVSC evaluation set

    mvsc_evaluation_set = \
        auxmvc_utils.ConvertKWAuxToMVSC(KWAuxGlitchTriggers=aux_glitch_vecs,
            KWAuxCleanTriggers=aux_clean_vecs)

    # save MVSC evaluation set in file........

    auxmvc_utils.WriteMVSCTriggers(mvsc_evaluation_set,
                                   output_filename=output_file_name,
                                   Classified=False)

    return mvsc_evaluation_set


def build_auxmvc_vector_at_gps(
    gps,
    trigger_dict,
    main_channel,
    time_window,
    output_file_name,
    channels=None,
    unsafe_channels=None,
    ):
    """
        build auxmvc vectors at the specified gps times regardless of what is actually in the main channel
...."""

    if isinstance(gps, (int, float)):
        gps = [gps]
    gps_start_time = min(gps)
    gps_end_time = max(gps)

    print 'Number of triggers in the main channel before thresholding:', \
        len(gps)

        # use only channels from the channels file, if provided

    if channels:
        selected_channels = event.read_channels_from_file(channels)
        trigger_dict.keep_channels(selected_channels)

                # to ensure consistency in dimensionality of auxmvc vectors
                # add the channels from the selected channels that are absent in trigger_dict

        for channel in selected_channels:
            if not channel in trigger_dict.channels():
                trigger_dict[channel] = []

        # get rid of unsafe channels if provided

    if unsafe_channels:
        unsafe_channels = event.read_channels_from_file(unsafe_channels)
        trigger_dict.remove_channels(unsafe_channels)

        # keep only the triggers from the [gps_start_time, gps_end_time] segment

        # first keep all triggers from the segment expanded by the time concidence window, so that not to loose coincidences

    trigger_dict.include([[gps_start_time - time_window, gps_end_time
                         + time_window]])

        # construct glitch auxmvc vectors

    aux_glitch_vecs = event.build_auxmvc_vectors_at_gps(gps,
            trigger_dict, main_channel=main_channel,
            coincidence_time_window=time_window)

      # only here for syntactic reasons fucntion call to ConvertKWAuxToMVSC

    aux_clean_vecs = event.build_auxmvc_vectors_at_gps([],
            trigger_dict, main_channel=main_channel,
            coincidence_time_window=time_window)

        # convert glitch and clean auxmvc vectors into MVSC evaluation set

    mvsc_evaluation_set = \
        auxmvc_utils.ConvertKWAuxToMVSC(aux_glitch_vecs, aux_clean_vecs)

        # save MVSC evaluation set in file

    auxmvc_utils.WriteMVSCTriggers(mvsc_evaluation_set,
                                   output_filename=output_file_name,
                                   Classified=False)

    return mvsc_evaluation_set


def execute_build_auxmvc_vectors(
    cp,
    execute_dir,
    trigdir,
    main_channel,
    output_file,
    gps_start_time,
    gps_end_time,
    channels=None,
    unsafe_channels=None,
    dq_segments=None,
    dq_segments_name='',
    ):
    """
....Submits the job that builds auxmvc feature vectors. Vectors are saved in the output file with .pat extension.
....Waits until the job is finished, returns its exit status and the output file name.
....execute_dir is the absolute path of the directory in which the job should be executed. 
...."""

    # initiate build_auxmvc_vectors job object

    build_auxmvc_vectors_job = auxmvc.build_auxmvc_vectors_job(cp,
            main_channel, channels=channels,
            unsafe_channels=unsafe_channels)

    # create node for this job

    build_auxmvc_vectors_node = auxmvc.build_auxmvc_vectors_node(
        build_auxmvc_vectors_job,
        trigdir,
        gps_start_time,
        gps_end_time,
        output_file,
        dq_segments=dq_segments,
        dq_segments_name=dq_segments_name,
        )

    # get full command line for this job

    build_auxmvc_vectors_command = \
        auxmvc.construct_command(build_auxmvc_vectors_node)

    print build_auxmvc_vectors_command

    # submit process

    exit_status = submit_command(build_auxmvc_vectors_command,
                                 'build_auxmvc_vectors', execute_dir,
                                 verbose=True)

    return (exit_status,
            build_auxmvc_vectors_node.get_output_files()[0])


def execute_prepare_training_auxmvc_samples(
    execute_dir,
    realtime_dir,
    cp,
    gps_start_time,
    gps_end_time,
    output_file,
    dq_segments='',
    dq_segments_name='',
    ):
    """
....Submits the job that prepares training auxmvc samples for MLAs using auxmvc feature vectors built during realtime
....evaluation step. Training samples are stored in the output file with .pat extension.
....Waits until the jobs is completed, returns its exit status and the output file name.
....realtime_dir is the absolute path to the directory with output of realtime evaluation step.   
...."""

    # initiate prepare_training_auxmvc_samples job object

    prepare_training_auxmvc_samples_job = \
        auxmvc.prepare_training_auxmvc_samples_job(cp)

    # create node for this job

    prepare_training_auxmvc_samples_node = \
        auxmvc.prepare_training_auxmvc_samples_node(
        prepare_training_auxmvc_samples_job,
        realtime_dir,
        gps_start_time,
        gps_end_time,
        output_file,
        dq_segments=dq_segments,
        dq_segments_name=dq_segments_name,
        )

    # get full command line for this job

    prepare_training_auxmvc_samples_command = \
        auxmvc.construct_command(prepare_training_auxmvc_samples_node)

    # submit process

    exit_status = \
        submit_command(prepare_training_auxmvc_samples_command,
                       'prepare_training_auxmvc_samples', execute_dir,
                       verbose=True)

    return (exit_status,
            prepare_training_auxmvc_samples_node.get_output_files()[0])


##########.... FOREST FUNCTIONS........####################################

def execute_forest_evaluate(
    patfile,
    trainedforest,
    ranked_file,
    cp,
    gps_start_time,
    gps_end_time,
    dir,
    ):
    """
....Submits job that evaluates samples of auxmvc feature vectors using random forest (MVSC).........
...."""

    # initiate use forest job

    use_forest_job = auxmvc.use_forest_job(cp)

    # create node for this job

    use_forest_node = auxmvc.use_forest_node(use_forest_job,
            trainedforest, patfile, ranked_file)

    # get full command line for this job

    use_forest_command = auxmvc.construct_command(use_forest_node)

    # submit process

    exit_status = submit_command(use_forest_command, 'forest_evaluate',
                                 dir)

    if exit_status == 0:

        # run postscript

        forest_add_excluded_variables_job = \
            auxmvc.forest_add_excluded_vars_job(cp)
        forest_add_excluded_variables_node = \
            auxmvc.forest_add_excluded_vars_node(forest_add_excluded_variables_job,
                patfile, ranked_file)

        # get full command line for this job

        forest_add_excluded_variables_command = \
            auxmvc.construct_command(forest_add_excluded_variables_node)

        # submit process

        exit_status = \
            submit_command(forest_add_excluded_variables_command,
                           'forest_add_excluded_variables', dir)

        return (exit_status, use_forest_node.get_output_files())
    else:
        return (exit_status, use_forest_node.get_output_files())


def execute_forest_train(
    training_samples_file,
    cache,
    cp,
    submit_dir,
    ):
    """
....Builds small dag to train  random forest (MVSC) and condor submits it.
...."""

    # set current directory to submit_dir

    os.chdir(submit_dir)

    # set path to condor log

    logpath = cp.get('idq_train', 'condorlogs')

    # set basename for condor dag

    basename = cp.get('general', 'ifo') + '_train_mvsc_' \
        + cp.get('general', 'usertag') + '-' \
        + training_samples_file.split('-')[-2] + '-' \
        + training_samples_file.split('-')[-1].split('.')[0]

    # creat directory for jobs .err and .out files

    if not os.path.exists('logs'):
        os.makedirs('logs')

    # initiate dag

    dag = auxmvc.auxmvc_DAG(basename, logpath)

    # get dag file

    dag_file = dag.get_dag_file()

    # initiate train forest job

    train_forest_job = auxmvc.train_forest_job(cp)

    # construct name for trained forest file

    trained_forest_filename = os.path.split(training_samples_file)[0] \
        + '/mvsc/' \
        + os.path.split(training_samples_file)[1].replace('.pat', '.spr'
            )

    # create node for this job

    train_forest_node = auxmvc.train_forest_node(train_forest_job,
            training_samples_file, trained_forest_filename)

    # append the node to dag

    dag.add_node(train_forest_node)

    # initiate add file to cache job

    add_file_to_cache_job = auxmvc.add_file_to_cache_job(cp)

    # create node for this job

    add_file_to_cache_node = \
        auxmvc.add_file_to_cache_node(add_file_to_cache_job,
            [train_forest_node.trainedforest], cache,
            p_node=[train_forest_node])

    # add the node to dag

    dag.add_node(add_file_to_cache_node)

    # write dag

    dag.write_sub_files()
    dag.write_dag()
    dag.write_script()

    # condor dag submit command....

    dag_submit_cmd = ['condor_submit_dag', dag_file]

    # submit dag

    exit_status = submit_command(dag_submit_cmd,
                                 'submit_forest_train_dag', submit_dir)

    return (exit_status, dag_file, train_forest_node.get_output_files())


############################## OVL functions ########################################....
# many of these functions simply delegate to ovl.py after building the appropriate params objects, etc. Essentially, these will be parsing/handling scripts that let ovl.py interface more smoothly with the idq code.

def ovl_evaluate(
    vetolist,
    GPS_headers=False,
    GPStimes=False,
    allvtrg=False,
    kw_trgfiles=False,
    gwchan=False,
    patfiles=False,
    skip_lines=1,
    filename='ovl_predict',
    output_dir='./',
    ):
    """
....takes the last line of vetolist_cache file and uses that as a pointer to the most recent vetolist information.
....generates vetosegment from kw_trgfiles
....if patfiles:
........generates GPStimes from patfiles
....else:
........generates GPStimes from kw_trgfiles using gwchan (REQUIRED)

....expects vetolsit to be a file path
........kw_trgfiles, patfiles to be lists of file paths
........gwchan to be a string
...."""

    if not GPStimes:

        # generate list of GPStimes to classify

        GPStimes = []
        if patfiles:
            for patfile in patfiles:
                GPStimes += ovl.patfile_to_GPStimes(patfile,
                        skip_lines=skip_lines)
        elif gwchan:
            for kw_trgfile in kw_trgfiles:
                GPStimes += ovl.kw_trigfile_to_GPStimes(kw_trgfile,
                        gwchan)
        else:
            print 'insufficient information provided to classify. Please supply either patfiles or gwchan'
        GPStimes.sort(key=lambda line: line[0])
        GPS_headers = ['GPS', 'i']

    # generate ovl prediction output

    gps_tcent = GPS_headers.index('GPS')
    return ovl.predict(
        vetolist,
        GPS_headers,
        GPStimes,
        gps_tcent,
        allvtrg=allvtrg,
        kw_trgfiles=kw_trgfiles,
        predict_filename=filename,
        output_dir=output_dir,
        )


def ovl_train(
    gpsstart,
    gpsstop,
    cp,
    scisegs=False,
    vetosegs=False,
    output_dir='./',
    ):
    """ 
....builds an ovl.params object and launches ovl training jobs on the specified data.
....pulls many parameters from "cp" config object

...."""

    # build params object

    analysis_range = [gpsstart, gpsstop]

    # load from cp object

    auxdir = cp.get('general', 'snglchndir')
    gwdir = cp.get('general', 'snglchndir')
    gwchans = cp.get('general', 'gwchannel').split()
    gwthr = float(cp.get('general', 'gw_kwsignif_thr'))
    ifos = cp.get('ovl_train', 'ifos').split()
    gwsets = cp.get('ovl_train', 'gwsets').split()
    safety = cp.get('ovl_train', 'safety')
    windows = [float(l) for l in cp.get('ovl_train', 'windows').split()]
    thresholds = [float(l) for l in cp.get('ovl_train', 'thresholds'
                  ).split()]
    Psigthr = float(cp.get('ovl_train', 'Psigthr'))
    effbydtthr = float(cp.get('ovl_train', 'effbydtthr'))

# ....channels = False
# ....notused = []

    channels = cp.get('general', 'selected-channels')
    notused = cp.get('general', 'unsafe-channels')
    if channels == '':
        channels = False
    if notused != '':
        notused = [l.strip('\n') for l in open(notused, 'r'
                   ).readlines()]
    else:
        notused = []

    params = ovl.params(
        analysis_range,
        auxdir,
        gwdir,
        gwchans,
        gwthr,
        ifos,
        gwsets,
        scisegs=scisegs,
        vetosegs=vetosegs,
        channels=channels,
        notused=notused,
        windows=windows,
        thresholds=thresholds,
        Psigthr=Psigthr,
        effbydtthr=effbydtthr,
        safety=safety,
        )

    # double check that windows are not bigger than padding set in idq_realtime

    len_windows = len(windows)
    windows = [w for w in windows if w <= float(cp.get('idq_realtime',
               'padding'))]
    if len_windows > len(windows):
        print 'WARNING: ovl windows are not consistent with idq_realtime padding! %d windows were removed.' \
            % (len_windows - len(windows))

    # load training parameters from cp object

    num_runs = int(cp.get('ovl_train', 'num_runs'))
    incremental = int(cp.get('ovl_train', 'incremental'))

    # launch training job

    vetolists = ovl.train(
        params,
        num_runs=num_runs,
        incremental=incremental,
        output_dir=output_dir,
        verbose=False,
        write_channels=True,
        )

    # append ovl_train.cach with most recent vetolist
    # ovl_train_cache = cp.get("general", "ovl_cache")
    # f = open(ovl_train_cache, "a")
    # for vetolist in vetolists:
     #  print >>f, vetolist
    # f.close()

    return vetolists


# create veto timeseries from datefile glitch events
# will create square window timeseries with +/-window around each event

def datfile2timeseries(datfile, window=.1, fs=256):

    # slim_load_datfile(file, skip_lines=1, columns=[])
    # GPS i w unclean signif SNR rank
    # L-KW_TRIGGERS_mvsc_1043300000_1043400000-1043896320-64.dat

    datdata = slim_load_datfile(datfile, skip_lines=0,
                                columns='GPS i rank'.split())
    for (key, val) in datdata.items():
        datdata[key] = map(float, val)
    (ifo, name, tstart, dur) = os.path.basename(datfile)[:-4].split('-')
    livetime = float(dur)
    segment = [float(tstart), float(tstart) + livetime]
    ts = numpy.zeros(livetime * fs)
    sidx = numpy.argsort(datdata['rank'])
    for i in sidx:
        if datdata['i'][i] == 0:  # skip clean samples
            continue
        gps = datdata['GPS'][i]
        seg = [gps - window, gps + window]
        rank = datdata['rank'][i]
        istart = max(0, int(math.floor((seg[0] - segment[0]) * fs)))  # start index
        istop = min(len(ts), int(math.ceil((seg[1] - segment[0]) * fs)))  # 1 + stop index
        ts[istart:istop] = rank

    return ts


######################  SVM functions  #######################

def execute_svm_evaluate(
    cp,
    test_file,
    range_file,
    model_file,
    predict_file,
    dir,
    ):

    use_svm_job = auxmvc.use_svm_job(cp)
    use_svm_node = auxmvc.use_svm_node(
        use_svm_job,
        cp,
        test_file,
        range_file,
        model_file,
        predict_file,
        )
    use_svm_command = auxmvc.construct_command(use_svm_node)
    exit_status = submit_command(use_svm_command, 'svm_evaluate', dir)
    return exit_status


def execute_svm_train(
    cp,
    train_file,
    range_file,
    model_file,
    cache_file,
    submit_dir,
    ):
    """
....Builds small condor dag to train SVM and submits it.
...."""

    # set current directory to submit_dir

    os.chdir(submit_dir)

    # set path to condor log

    logpath = cp.get('idq_train', 'condorlogs')

    # set basename for condor dag

    basename = cp.get('general', 'ifo') + '_train_svm_' \
        + cp.get('general', 'usertag') + '-' + train_file.split('-'
            )[-2] + '-' + train_file.split('-')[-1].split('.')[0]

    # creat directory for jobs .err and .out files

    if not os.path.exists('logs'):
        os.makedirs('logs')

    # initiate dag

    dag = auxmvc.auxmvc_DAG(basename, logpath)

    # get dag file

    dag_file = dag.get_dag_file()

    train_data = os.path.split(train_file)[1]
    range_data = os.path.split(range_file)[1]
    model_data = os.path.split(model_file)[1]

    # cache_data = os.path.split(cache_file)[1]

    if not os.path.exists(range_file):
        os.mknod(range_file, 0644)
    if not os.path.exists(model_file):
        os.mknod(model_file, 0644)
    if not os.path.exists(train_data):
        os.symlink(os.path.abspath(train_file),
                   os.path.join(os.getcwd(), train_data))
    if not os.path.exists(range_data):
        os.symlink(os.path.abspath(range_file),
                   os.path.join(os.getcwd(), range_data))
    if not os.path.exists(model_data):
        os.symlink(os.path.abspath(model_file),
                   os.path.join(os.getcwd(), model_data))

    # initiate svm train job
    # if not os.path.exists(cache_data): os.symlink(os.path.abspath(cache_file),os.path.join(os.getcwd(),cache_data))

    train_svm_job = auxmvc.train_svm_job(cp)

    # create node for this job

    train_svm_node = auxmvc.train_svm_node(train_svm_job, cp,
            train_data, range_data, model_data)

    # append the node to dag

    dag.add_node(train_svm_node)

    # initiate add file to cache job

    add_file_to_cache_job = auxmvc.add_file_to_cache_job(cp)

    # create node for this job

    add_file_to_cache_node = \
        auxmvc.add_file_to_cache_node(add_file_to_cache_job,
            [submit_dir + train_svm_node.model_file, submit_dir
            + train_svm_node.range_file], cache_file,
            p_node=[train_svm_node])

    # add the node to dag

    dag.add_node(add_file_to_cache_node)

    # post_script = cp.get("svm_evaluate", "svm_train_post_script")
    # train_svm_node = auxmvc.train_svm_node(train_svm_job, dag, cp, train_data, range_data, model_data, post_script, cache_data, model_file, range_file)

    dag.write_sub_files()
    dag.write_dag()
    dag.write_script()

    # condor dag submit command....

    dag_submit_cmd = ['condor_submit_dag', dag_file]

    # submit dag

    exit_status = submit_command(dag_submit_cmd, 'submit_svm_train_dag'
                                 , submit_dir)

    return (exit_status, dag_file, train_svm_node.get_output_files())


##@}

