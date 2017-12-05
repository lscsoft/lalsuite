# Copyright (C) 2013 Lindy Blackburn, Reed Essick and Ruslan Vaulin
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

## \defgroup laldetchar_py_idq_idq iDQ Functions
## \ingroup laldetchar_py_idq
## Synopsis
# ~~~
# from laldetchar.idq import idq
# ~~~
# \author Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>), Lindy Blackburn (<lindy.blackburn@ligo.org>)


#=================================================

import os
import sys
import glob
import re

import logging
import time

import numpy
import math

import multiprocessing
import subprocess
import tempfile

from collections import defaultdict

import ConfigParser

from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
from glue.ligolw import table

from glue.ligolw import ilwd
from glue.ligolw import types as ligolwtypes
from glue.ligolw.utils import process
from glue.lal import LIGOTimeGPS

#from pylal import frutils
from pylal import Fr

from laldetchar.idq import event
from laldetchar.idq import ovl
from laldetchar.idq import auxmvc_utils
from laldetchar.idq import auxmvc
from laldetchar.idq import idq_tables
from laldetchar.idq import pdf_estimation as pdf_e

from laldetchar import git_version

#===================================================================================================

description = """ a placeholder for a cleaned up idq.py """

__prog__ = 'idq.py'

__author__ = \
    'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)'

__version__ = git_version.id
__date__ = git_version.date

## \addtogroup laldetchar_py_idq_idq
# @{

#===================================================================================================
### hard coded magic numbers dependent on classifiers, ETG, etc
traincache_nlines = { "ovl":1 , "forest": 1, "svm": 2, "ann": 3}

mla_flavors = ["forest", "svm", "ann"]

train_with_dag = ["forest", "svm", "ann"]

#===================================================================================================
# standardize uploads to GraceDb and tagging
#===================================================================================================

tagnames = ["idq"]

#===================================================================================================
# forking
#===================================================================================================

def double_fork( cmd, stdout=None, stderr=None, cwd='.' ):
    '''
    double forks your process so it is orphaned

    NOTE: does not allow you to communicate through STDIN. If you need to do that, you should handle the calls yourself.
    '''
    proc = multiprocessing.Process(target=fork, args=(cmd, stdout, stderr, cwd))
    proc.start() ### start the job
    proc.join() ### wait for it to finish, which should be quick because it only forks

def fork( cmd, stdout=None, stderr=None, cwd='.' ):
    '''
    forks the process via subprocess
    returns the pid

    NOTE: does not allow you to communicate through STDIN. If you need to do that, you should handle the calls to subprocess yourself
    '''
    return subprocess.Popen( cmd, stdout=stdout, stderr=stderr, cwd=cwd )

#===================================================================================================
# parsing the config file
#===================================================================================================

def config_to_combinersD( config ):
    if not config.has_option('general','combiners'):
        return {}, []

    combiners = sorted(set(config.get('general', 'combiners').split()))
    classifiers = sorted(set(config.get('general', 'classifiers').split()))
    for classifier in classifiers:
        if not config.has_section(classifier):
            raise ValueError("classifier=%s does not have a corresonding section in config file"%classifier)

    combinersDict = {}
    referenced_classifiers = []
    for combiner in combiners:
        if config.has_section(combiner):
            combinersDict[combiner] = dict( config.items(combiner) )
            combinersDict[combiner]['flavor'] = None

            if combinersDict[combiner].has_key('classifiers'): ### check classifiers
                these_classifiers = combinersDict[combiner]['classifiers'].split()
                for classifier in these_classifiers:
                    if classifier not in classifiers:
                        raise ValueError("classifier=%s referenced by combiner=%s, but is not in the list of classifiers to be run"%(classifier, combiner))
                combinersDict[combiner]['classifiers'] = these_classifiers

                columns = []
                for classifier in these_classifiers: ### add section for each classifier
                    combinersDict[combiner][classifier] = {}
                    columns += ["%s_rank"%classifier, "%s_p(g)"%classifier, "%s_p(c)"%classifier]
                combinersDict[combiner]['columns'] = columns               

                referenced_classifiers += these_classifiers
            else:
                raise ValueError("combiner=%s does not have an option=\"classifiers\"")

        else:
            raise ValueError("combiner=%s does not have a corresponding section in config file"%(combiner))

    return combinersDict, sorted(set(referenced_classifiers))

def config_to_classifiersD( config ):
    classifiers = sorted(set(config.get('general', 'classifiers').split()))

    classifiersDict = {}
    mla = False
    ovl = False
    for classifier in classifiers:
        if config.has_section(classifier):
            classifiersDict[classifier] = dict( config.items(classifier) )
            flavor = classifiersDict[classifier]['flavor']

            classifiersDict[classifier]['config'] = classifier_specific_config( config, classifier, flavor )
            classifiersDict[classifier]['mla'] = (flavor in mla_flavors)

            mla = mla or classifiersDict[classifier]['mla']
            ovl = ovl or (flavor == 'ovl')

        else:
            raise ValueError("classifier=%s does not have a corresponding section in config file"%(classifier))

    return classifiersDict, mla, ovl

def classifier_specific_config( config, classifier, flavor ):

    ### NOTE: we require all classifiers to have a "condor" option, even though not all of them may use condor...
    ### this is done for simplicity within this function and throughout, but may seem unnecessary when writing the config object

    if config.has_option('classifier', 'config'):
        raise ValueError("classifier=%s has option 'config' that will be overwritten. Please remove that option"%classifier)
   
    cp = ConfigParser.SafeConfigParser()

    if flavor == "ovl":
        eval_section = 'ovl_evaluate'
        train_section = 'train_ovl'

        ### fill in ovl sections here? Not needed because OVL doesn't train with a dag...
        ### this is just here as a "place-holder"
        cp.add_section(eval_section)
        cp.add_section(train_section)

    elif flavor == "forest":
        eval_section = 'forest_evaluate'
        train_section = 'train_forest'

        cp.add_section(eval_section)
        for option in "A a z v".split():
            if config.has_option(classifier, option):
                cp.set( eval_section, option, config.get(classifier, option) )

        cp.add_section(train_section)
        for option in "a n l s c g i d z".split():
            if config.has_option(classifier, option):
                cp.set( train_section, option, config.get(classifier, option) )

    elif flavor == "svm":
        eval_section = 'svm_evaluate'
        train_section = 'train_svm'

        cp.add_section(eval_section)
        for option in "scale predict rank".split():
            if config.has_option(classifier, option):
                cp.set( eval_section, option, config.get(classifier, option) )

        cp.add_section(train_section)
        for option in "scale train rank gamma cost".split():
            if config.has_option(classifier, option):
                cp.set( train_section, option, config.get(classifier, option) )

    elif flavor == "ann":
        convert_section = 'ann_convert'
        eval_section = 'ann_evaluate'
        train_section = 'train_ann'

        cp.add_section(convert_section)
        for option in "normalization-attributes transform-dt-function".split():
            if config.has_option(classifier, option):
                cp.set( convert_section, option, config.get(classifier, option) )

        cp.add_section(eval_section)
        for option in "training-machine".split():
            if config.has_option(classifier, option):
                cp.set( eval_section, option, config.get(classifier, option) )

        cp.add_section(train_section)
        for option in "training-machine hidden-neurons connection-rate steep-out max-epochs weights-min weights-max increase-factor desired-error".split():
            if config.has_option(classifier, option):
                cp.set( train_section, option, config.get(classifier, option) )


    else:
        raise ValueError('flavor=%s not understood for classifier=%s'%(flavor, classifier))

    if not config.has_option(classifier, 'condor'):
        raise ValueError('classifier=%s section has no option=condor'%classifier)
    condor_section = config.get(classifier, 'condor')
    if not config.has_section(condor_section):
        raise ValueError('config has no section=%s'%condor_section)

    cp.add_section('condor')
    for option, value in config.items(condor_section):
        cp.set('condor', option, value)

    ### add accounting tags if present
    if config.has_option('general', 'accounting_group'):
        cp.set('condor', 'accounting_group', config.get('general', 'accounting_group') )
    if config.has_option('general', 'accounting_group_user'):
        cp.set('condor', 'accounting_group_user', config.get('general', 'accounting_group_user') )

    cp.add_section('idq_train')
    cp.set('idq_train', 'condorlogs', config.get('train', 'condorlogs'))

    cp.add_section('general')
    for option, value in config.items('general'):
        cp.set('general', option, value)

    return cp

#=================================================
# standardized names
#=================================================

def channame(ifo, classifier, tag):
    return "%s:%s%s"%(ifo, classifier, tag)

def cache(directory, classifier, tag):
    return "%s/%s%s.cache"%(directory, classifier, tag)

def pat(directory, ifo, tag, t, stride):
    return "%s/%s_mla%s-%d-%d.pat"%(directory, ifo, tag, t, stride)

def dat(directory, classifier, ifo, trained, tag, t, stride):
    return "%s/%s_%s_%s%s-%d-%d.dat"%(directory, ifo, classifier, trained, tag, t, stride)

def xml(directory, classifier, ifo, trained, calib, tag, t, stride):
    return "%s/%s_idq_%s_%s-%s%s-%d-%d.xml.gz"%(directory, ifo, classifier, trained, calib, tag, t, stride)

def gdb_xml(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_idq_%s%s-%d-%d.xml.gz"%(directory, ifo, classifier, tag, t, stride)

def timeseries(directory, classifier, ifo, trained, calib, tag, t, stride):
    return "%s/%s_idq_%s_%s-%s%s-%d-%d.npy.gz"%(directory, ifo, classifier, trained, calib, tag, t, stride)

def gdb_timeseries(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_idq_%s%s-%d-%d.npy.gz"%(directory, ifo, classifier, tag, t, stride)

def timeseriesgwf(directory, classifier, ifo, trained, calib, tag, t, stride):
    return "%s/%s_idq_%s_%s-%s%s-%d-%d.gwf"%(directory, ifo, classifier, trained, calib, tag, t, stride)

def gdb_timeseriesgwf(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_idq_%s%s-%d-%d.gwf"%(directory, ifo, classifier, tag, t, stride)

def segxml(directory, tag, t, stride):
    return "%s/science_segments%s-%d-%d.xml.gz"%(directory, tag, t, stride)

def segascii(directory, tag, t, stride):
    return "%s/science_segments%s-%d-%d.seg"%(directory, tag, t, stride)

def idqsegascii(directory, tag, t, stride):
    return "%s/idq_segments%s-%d-%d.seg"%(directory, tag, t, stride)

def roc(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_%s%s-%d-%d.roc"%(directory, ifo, classifier, tag, t, stride)

def uroc(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_%s%s-%d-%d.uroc"%(directory, ifo, classifier, tag, t, stride)

def calib_check(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_%s_calib%s-%d-%d.txt"%(directory, ifo, classifier, tag, t, stride)

def chan_perf(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_%s_chan-perf%s-%d-%d.txt"%(directory, ifo, classifier, tag, t, stride)

def sumhtml(directory, tag, t, stride):
    return "%s/summary-index%s-%d-%d.html"%(directory, tag, t, stride)

def kdename(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_%s%s_KDE-%d-%d.npy.gz"%(directory, ifo, classifier, tag, t, stride)

def gdb_summary(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_idq_%s_summary%s-%d-%d.txt"%(directory, ifo, classifier, tag, t, stride)

def gdb_minFap_json(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_%s%s-%d-%d.json"%(directory, ifo, classifier, tag, t, stride)

def gdb_ovlstripchart_json(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_%s_chanlist%s-%d-%d.json"%(directory, ifo, classifier, tag, t, stride)

def gdb_roc_json(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_%s%s_ROC-%d-%d.json"%(directory, ifo, classifier, tag, t, stride)

def gdb_calib_json( directory, ifo, classifier, tag, t, stride):
    return "%s/%s_%s%s_calib-%d-%d.json"%(directory, ifo, classifier, tag, t, stride)

def useSummary_json( directory, ifo, classifier, tag, t, stride):
    return "%s/%s_%s%s-%d-%d.json"%(directory, ifo, classifier, tag, t, stride)

#def frame2segment( directory, classifier, ifo, FAPthr, tag, right_padding, left_padding, t_lag, widen, t, stride ):
#    return "%s/%s_%s_FAP-%.3e_rp-%.3f_lp-%.3f_tl-%.3f_wd-%.3f%s-%s-%d.seg"%(directory, ifo, classifier, FAPthr, right_padding, left_padding, t_lag, widen, tag, t, stride)
def frame2segment( directory, classifier, ifo, FAPthr, tag, right_padding, left_padding, t_lag, t, stride ):
    return "%s/%s_%s_FAP-%.3e_rp-%.3f_lp-%.3f_tl-%.3f%s-%d-%d.seg"%(directory, ifo, classifier, FAPthr, right_padding, left_padding, t_lag, tag, t, stride)

def frame2segmentxml( directory, tag, t, stride ):
    return "%s/laldetcharIDQFrame2Segment%s-%d-%d.xml.gz"%(directory, tag, t, stride)

#=================================================
# extract start/dur
#=================================================

def extract_start_stop(filename, suffix='.dat'):
    """
    picks of start, dur information from filename and returns [start, stop]
    """
    matchfile = re.compile('.*-([0-9]*)-([0-9]*)\%s$' % suffix)
    m = matchfile.match(filename)
    (st, dur) = (int(m.group(1)), int(m.group(2)))
    return [st, st+dur]

#=================================================
# extract classifier name
#=================================================

def extract_dat_name(datfilename):
    """
    WARNING: will break if classifier name contains "_"
    """
    return datfilename.split("/")[-1].split("_")[1]

def extract_fap_name(fapfilename):
    """
    WARNING: will break if classifier name contains "_"
    """
    return fapfilename.split("/")[-1].split("_")[2]

def extract_roc_name(rocfilename):
    """
    WARNING: will break if classifier name contains "_"
    """
    return rocfilename.split("/")[-1].split("_")[1]

def extract_xml_name(xmlfilename):
    """
    WARNING: will break if classifier name contains "_"
    """
    return xmlfilename.split("/")[-1].split("_")[2]

#=================================================
# trained and calibration ranges
#=================================================

def extract_timeseries_ranges( filename ):
    """
    returns trained_range, calib_range
    """
    ans = filename.split("/")[-1].split("_")
    tend, cstart = ans[4].split('-')
    trained = [int(ans[3]), int(tend)]
    calib = [int(cstart), int(ans[5])]
    return trained, calib

def extract_trained_ranges(lines, flavor):
    return [extract_trained_range( line, flavor ) for line in lines]

def extract_trained_range(lines, flavor):
    if flavor == "forest":
        trainedforest = lines[0]
        ### getting start and end of training period from the name of the trained forest file
        mvsc_start_training_period = int(trainedforest.split('-')[-2])
        mvsc_end_training_period = mvsc_start_training_period + int(trainedforest.split('-')[-1].split('.')[0])
        trained_range = "%d_%d"%(mvsc_start_training_period, mvsc_end_training_period)

    elif flavor == "ovl":
        trained_range = [v for v in lines[0].split('/') if v != ''][-4]  # we assume a particular directory structure

    elif flavor == "svm":
        svm_model, svm_range_file = lines

        ## start and end of training period from the name of the svm model file
        svm_start_training_period = int(svm_model.split('-')[-2])
        svm_end_training_period = svm_start_training_period + int(svm_model.split('-')[-1].split('.')[0])
        trained_range = "%d_%d"%(svm_start_training_period,svm_end_training_period)

    elif flavor == "ann":
        trained_ann = lines[0]
        ### getting start and end of training period from the name of the trained forest file
        ann_start_training_period = int(trained_ann.split('-')[-2])
        ann_end_training_period = ann_start_training_period + int(trained_ann.split('-')[-1].split('.')[0])
        trained_range = "%d_%d"%(ann_start_training_period, ann_end_training_period)


    else:
        raise ValueError("do not know how to extract trained range for flavor=%s"%flavor)

    return trained_range

def extract_calib_ranges(urocs):
    return [extract_calibration_range(uroc) for uroc in urocs]

def extract_calib_range(uroc):
    return "%d_%d"%tuple(extract_start_stop( uroc, suffix=".uroc" ))
#    calib = uroc.split("/")[-1].split("-")[-3].split("_")[-1]
#    return calib

def extract_kde_ranges( kde_names ):
    return [extract_kde_range( kde_name ) for kde_name in kde_names ]

def extract_kde_range( kde_name ):
    return "%d_%d"%tuple(extract_start_stop( kde_name, suffix=".npy.gz" ))

def best_range(t, ranges):
    """
    finds the best range for time t
    """
    dists = []
    for r in ranges:
        start, stop = [int(l) for l in r.split("_")]
        if start <= gps and gps <= stop:
            dists.append( (r, start-gps, start, stop-start) )

    dists.sort(key=lambda l: l[2])  # sort by start times, earlier times at the top
    dists.sort(key=lambda l: l[3], reverse=True)  # sort by durations, longer durations at the top
    dists.sort(key=lambda l: l[1])  # sort by distance, smaller distances at the top

    return dists[0][0]

#=================================================
# current-time and scheduling
#=================================================
gpsref = time.mktime(time.strptime('Tue Jan 06 00:00:00 1980')) - time.timezone  # GPS ref in UNIX time

def nowgps():
    """
    return current GPS time, assume 16 leap-seconds
    """
    return time.time() - gpsref + 16  # add in leap-seconds

#=================================================
# managing duplicate processes
#=================================================

def dieiflocked(lockfile='.idq.lock'):
    """
    try to set pid file, die iif file exists and pid within is still running
    Note, there's still the possibility of a small race condition here, but it shouldn't matter in practice
    """
    if os.path.exists(lockfile):
        lockfp = open(lockfile, 'r')
        pid = int(lockfp.readline().strip()) ### read the pid from the file
        lockfp.close()

        try: ### confirm the process is running
            os.kill(pid, 0) ### sends 0 to the process. does nothing if it's running. raises an error if it isn't
            ### we'll only execute the next line if the call to os.kill did not raise an error -> process exists
            sys.exit('ERROR: possible duplicate process. %s contains running pid=%d. exiting..'%(lockfile, pid))

        except OSError: ### process is not running
            pass
        
    ### we only get here if there is no currently running process
    lockfp = open(lockfile, 'w')
    lockfp.write( str(os.getpid()) )
    lockfp.close()

def release(lockfile='.idq.lock'):
    """
    try to release a pid file (just delete it)
    """
    os.remove(lockfile)
    
#=================================================
### reporting/logging utilities
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


def setup_logger(name, logfile, stdout, format='%(asctime)s %(message)s'):
    """
    sets up the standard logger for all iDQ processes
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter(format)

    hdlr1 = logging.StreamHandler(stdout)
    hdlr1.setFormatter(formatter)
    hdlr1.setLevel(logging.INFO)
    logger.addHandler(hdlr1)

    hdlr2 = logging.FileHandler(logfile)
    hdlr2.setFormatter(formatter)
    hdlr2.setLevel(logging.INFO)
    logger.addHandler(hdlr2)

    return logger

#========================
# log parsing functionality
#========================

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
        stride_start = -numpy.infty ### we wait until we find a new stride or time out
        stride_end = -numpy.infty

    file_obj.seek(starting_tell, 0) ### go back to where we started

    return (stride_start, stride_end)

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
            stride_start = float( line.split("Begin: stride ")[-1].split("-")[-2] )

        waited = 0.0 ### we found a new line, so re-set this counter

    ### return statement conditioned on how the loop exited
    return stride_start > t, waited > max_wait, time.time()-to > timeout

def get_condor_dag_status(dag):
    """
    Check on status of the condor dag by parsing the last line of .dagman.out file.
    """

    dag_out_file = dag + '.dagman.out'
    if not os.path.exists(dag_out_file):
        return 'incomplete'

    lastline = open(dag_out_file, 'r').readlines()[-1].strip('\n')

    if 'EXITING WITH STATUS' in lastline:
        exit_status = int(lastline.split('EXITING WITH STATUS')[-1])

    else:
        exit_status = 'incomplete'

    return exit_status

#=================================================
# cachefiles
#=================================================

class Cachefile():
    """
    interface between a cachefile and the pipeline
    """

    def __init__(self, name):
        self.name = name
        if not self.exists():
            self.touch()
        self.set_time()

    def exists(self):
        """
        checks to see whether the file exists
        """
        return os.path.exists(self.name)

    def is_empty(self):
        """
        checks whether this file is empty
        """
        return len(open(self.name, 'r').readlines()) == 0

    def touch(self):
        """
        touch the file
        """
        base = os.path.dirname(self.name)
        if not os.path.exists(base):
            os.makedirs(base)

        return open(self.name, 'a').close()
		
    def append(self, string):
        """
        append "string" to this cachefile
        """
        file_obj = open(self.name, 'a')
        print >> file_obj, string
        file_obj.close()

        self.set_time()

    def timestamp(self):
        """
        get the current timestamp of this file
        """
        return int(os.path.getmtime(self.name))

    def was_modified(self):
        """
        determines whether this file has been modified
        """
        return (self.timestamp() != self.time)

    def set_time(self):
        self.time = self.timestamp()

    def tail(self, nlines=1):
        """
        reads the last nlines of the cache file
        """
        return [line.strip() for line in open(self.name, 'r').readlines()[-nlines:]]

    def readlines(self):
        return [line.strip() for line in self.readlines()]

#=================================================
# general I/O
#=================================================
def slim_load_datfile(file, skip_lines=1, columns=[]):
    """ 
    loads only the given columns from a *.dat file. assumes the variable names are given by the first line after skip_lines left so that if an element in columns is not in the variables list, an error will be raised
    """

    output = dict([(c, []) for c in columns])
    f_lines = open(file).readlines()[skip_lines:]

    if f_lines:
        # find variable names
        variables = dict([(line, i) for (i, line) in enumerate(f_lines[0].strip('\n').split())])

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

def slim_load_datfiles(files, skip_lines=1, columns=[]):
    """
    delegates to slim_load_datfile, and combines the total output into a single structure
    """
    output = dict( (c, []) for c in columns )
    for file in files:
        o = slim_load_datfile(file, skip_lines=skip_lines, columns=columns)
        for c in columns:
            output[c] += o[c]

    return output

def output_to_datfile( output, dat ):
    """
    write the info stored in output into a file
    """
    columns = output.keys()
    file_obj = open(dat, "w")
    
    if columns:
        s = ""
        for c in columns:
            s += " %s"%str(c)
        print >> file_obj, s

        for i in xrange(len(output[columns[0]])):
            s = ""
            for c in columns:
                s += " %s"%str(output[c][i])
            print >> file_obj, s

    file_obj.close()

def sort_output( output, key ):
    if not output.has_key(key):
        raise ValueError("output does not have key=%s"%key)

    columns = output.keys()

    out = [ [output[key][i]] + [output[column][i] for column in columns] for i in xrange(len(output[key])) ]
    if not len(out):
        return output

    out.sort( key=lambda l: l[0] )
    out = numpy.array(out)
    for ind, column in enumerate(columns):
        output[column] = list(out[:,1+ind])

    return output

def filter_datfile_output( output, segs ):
    """
    filters the data in output by segs
    requres "GPS" to be a key in output
    """
    if not output.has_key("GPS"):
        raise ValueError("output must have key: GPS")

    if not output['GPS']: ### empty list
        return output

    columns = output.keys()

    output = sort_output( output , 'GPS' ) ### make sure these are sorted by GPS

    out = numpy.array(event.include( [ [ float(output['GPS'][i]) ] + [output[column][i] for column in columns] for i in xrange(len(output['GPS'])) ], segs, tcent=0 ))

    if not len(out): ### nothing survived
        return dict( (column, []) for column in columns )

    for ind, column in enumerate(columns):
        if column == 'GPS':
            output[column] = out[:,1+ind].astype(float)
        elif column == 'i':
            output['i'] = out[:,1+ind].astype(float)
        elif column == 'rank':
            output['rank'] = out[:,1+ind].astype(float)
        else:
            output[column] = out[:,1+ind]

    return output

def cluster_datfile_output( output, cluster_key='signif', cluster_win=1.0 ):
    if not output.has_key("GPS"):
        raise ValueError("output must have \"GPS\" as a key")
    if not output.has_key(cluster_key):
        raise ValueError("output must have cluster_key=\"%s\" as a key"%cluster_key)

    ### separate output into glitches and cleans
    columns, glitches, cleans = separate_output( output )

    ### cluster glitches and cleans separately
    glitches = cluster( glitches, columns, cluster_key=cluster_key, cluster_window=cluster_win )

    cleans   = cluster( cleans  , columns, cluster_key=cluster_key, cluster_window=cluster_win ) 

    ### combine into output format
    return combine_separated_output( columns, [glitches, cleans] )

def separate_output( output ):
    """
    returns lists of glitches and cleans in output
    """
    if not output.has_key('i'):
        raise ValueError("output must have \"i\" as a key")

    columns = output.keys()

    ### separate into glitches and cleans
    glitches = []
    cleans = []
    N = len(output['i'])
    for ind in xrange(N):
        i = output['i'][ind]
        if i:
            glitches.append( [output[c][ind] for c in columns] )
        else:
            cleans.append( [output[c][ind] for c in columns] ) 

    return dict( (c,j) for j, c in enumerate(columns) ), glitches, cleans

def combine_separated_output( columns, event_lists ):
    output = dict( (c,[]) for c in columns.keys() )
    for events in event_lists:
        for trg in events:
            for c, i in columns.items():
                output[c].append( trg[i] )
    return output

def cluster( samples, columns, cluster_key='signif', cluster_window=1.0 ):
    """
    Clustering performed with the sliding window cluster_window keeping trigger with the highest rank;
    glitches is an array of glitches, and we sort it to be in ascending order (by GPS time)
    Clustering algorithm is borrowed from pylal.CoincInspiralUtils cluster method.
    """

    # sort glitches so they are in the propper order
    samples.sort(key=lambda line: line[columns['GPS']])

    # initialize some indices (could work with just one)
    # but with two it is easier to read
    this_index = 0
    next_index = 1
    while next_index < len(samples):

        # get the time for both indices
        thisTime = samples[this_index][columns['GPS']]
        nextTime = samples[next_index][columns['GPS']]

        # are the two coincs within the time-window?
        if nextTime - thisTime < cluster_window:

            # get the ranks
            this_rank = samples[this_index][columns[cluster_key]]
            next_rank = samples[next_index][columns[cluster_key]]

            # and remove the trigger which has the lower rank
            if next_rank > this_rank:
                del samples[this_index]  # glitches = numpy.delete(glitches, this_index, 0)
            else:
                del samples[next_index]  # glitches = numpy.delete(glitches, next_index, 0)

        else:
            # NOTE: we don't increment this_index and next_index if we've removed glitches
            # the two triggers are NOT in the time-window
            # so must increase index
            this_index += 1
            next_index += 1

    return samples

def extract_dq_segments(xmlfile, dq_name):
    """
    Loads xmlfile containing dq segments and extracts segments of the type set by dq_name ( e.g. L1:DMT-SCIENCE:3).
    dq_name is what was given as an --include-segments option to ligolw_segment_query that generated the xmlfileobj"
    """

    lsctables.use_in(ligolw.LIGOLWContentHandler)
    if type(xmlfile) == str:
        xmldoc = ligolw_utils.load_filename(xmlfile, contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))  # load as filename
    else:
        xmldoc = ligolw_utils.load_fileobj(xmlfile, contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))[0]  # laod as file object

    # get segment tables
    sdef = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
    ssum = table.get_table(xmldoc, lsctables.SegmentSumTable.tableName)
    seg = table.get_table(xmldoc, lsctables.SegmentTable.tableName)

    # segment definer ID corresponding to dq_name, which was use in segment data base query
    id = next(a.segment_def_id for a in sdef if (a.name == dq_name.split(':')[1]) ) ### fragile because based on naming convention

    # list of covered segments (e.g. where SCIENCE mode is defined)
    covered = [[a.start_time, a.end_time] for a in ssum if (a.segment_def_id == id)]

    # final RESULT segment from query
    result_id = next(a.segment_def_id for a in sdef if (a.name == 'RESULT'))
    good = [[a.start_time, a.end_time] for a in seg if (a.segment_def_id == result_id) ]

    return (good, covered)

#=================================================
# ROC manipulation via output
#=================================================

def dat_to_rcg( output ):
    """
    returns ranks, cum_num_cln, cum_num_gch based on the data in output
    assumes output is like the outarg from slim_load_datfile
    output must have keys: rank, i
    if output.has_key('weight'): we use those
    else: we set all weights = 1
    """
    if not output.has_key('weight'):
        output['weight'] = numpy.ones(len(output['rank']), dtype='float')

    tmp = zip( output['rank'], output['i'], output['weight'] )
    tmp.sort(key=lambda l: l[0], reverse=True) ### sort by decending rank
    r = []
    c = []
    g = []
    n_c = 0
    n_g = 0
    _rank = 1.0
    for rank, i, w in tmp:
        if rank != _rank: ### new rank is different
            r.append( _rank )
            c.append( n_c )
            g.append( n_g )
            _rank = rank
        n_c += (1-i)*w
        n_g += i*w
    r.append( _rank )
    c.append( n_c )
    g.append( n_g )

    return numpy.array( (r, c, g) )

def rcg_to_file(ROCfile, r, c, g):
    """
    prints rcg to file
    """
    f = open(ROCfile, "w")
    print >> f, g[-1]
    print >> f, c[-1]
    for R, C, G in zip(r, c, g):
        print >> f, R, C, G
    f.close()

def file_to_rcg(ROCfile):
    """ 
    reads in a standard ROC file and returns lists of ranks, cumulative_cleans, cumulative_glitches 
    """
    ranks = []
    c_cln = []
    c_gch = []
    tot_cln = False
    tot_gch = False
    i = 0  # counter to find first lines not commented
    f = open(ROCfile, 'r')
    for line in f:
        if line[0] != '#':
            if i == 0:
                tot_gch = float(line)
                i += 1
            elif i == 1:
                tot_cln = float(line)
                i += 1
            else:
                fields = line.strip().split()
                ranks.append(float(fields[0]))
                c_cln.append(float(fields[1]))
                c_gch.append(float(fields[2]))

    return (ranks, c_cln, c_gch, tot_cln, tot_gch)

def resample_rcg(uranks, r, c, g):
    """
    returns the cum_num_cln and cum_num_gch sampled at uranks rather than r
    """
    N = len(r)
    ucln = []
    ugch = []
    ind = 0
    cc = 0  # current value of c_cln corresponding to 'ind'
    cg = 0
    for ur in uranks:
        while ind < N and r[ind] >= ur:
            cc = c[ind]
            cg = g[ind]
            ind += 1
        ucln.append(cc)
        ugch.append(cg)

    return uranks, ucln, ugch

def rcg_to_diff( c, g):
    """
    returns list corresponding to the number of glitches, cleans at each element in the list, instead of cumulative counts
    """
    dc = []
    dg = []
    _c = 0
    _g = 0
    for C, G in zip(c, g):
        dc.append( C-_c )
        dg.append( G-_g )
        _c = C
        _g = G
    return dc, dg

def rcg_to_EffFAPmap(rs, c_cln, c_gch, kind='linear'):
    """
    generates both FAP and Eff maps from r, c, g
      return Effmap, FAPmap
    """
    from scipy.interpolate import interp1d

    tot_cln = c_cln[-1]
    tot_gch = c_gch[-1]

    #ensure interpolation range covers allowable range for rank
    if rs[0] != 1.0:
        rs.insert(0, 1.0)
        c_gch.insert(0, 0.0)
        c_cln.insert(0, 0.0)
    if rs[-1] != 0.0:
        rs.append(0.0)
        c_gch.append(tot_gch)
        c_cln.append(tot_cln)

    ### construct Effmap object
    effmap = BinomialMap( rs, c_gch, kind=kind )

    ### construct FAPmap object
    fapmap = BinomialMap( rs, c_cln, kind=kind )

    return (effmap, fapmap)

def kde_pwg(eval, r, ds, scale=0.1, s=0.01):
    """
    delegates to pdf_e.point_wise_gaussian_kde after constructing the correct input values
    """
    observ = []
    for R, dS in zip(r, ds):
        observ += [R]*dS

    if numpy.sum(ds):
        return pdf_e.point_wise_gaussian_kde(eval, observ, scale=scale, s=s)
    else:
        return numpy.ones_like(eval)

def kde_to_ckde( s ):
    c = numpy.zeros_like(s)
    _s = 0.0
    cum = 0.0
    for i, S in enumerate(s[::-1]): ### integrate from the right
        cum += 0.5*(S+_s)
        c[i] = cum
        _s = S
    c /= c[-1]
    return c[::-1] ### reverse to get the order correct

def bin_by_rankthr(rankthr, output, columns=None):

    if columns==None:
        columns = sorted(output.keys())

    gch = []
    cln = []
    for i in xrange(len(output['i'])):
        if output['rank'][i] >= rankthr:
            d = dict( (column, output[column][i]) for column in columns )
            if output['i'][i]:
                gch.append( d )
            else:
                cln.append( d )

    return cln, gch

#=================================================
# channel performance data
#=================================================

def dat_to_perf( output ):
    vchans = sorted(set(output['vchan']))
    performance = dict( (vchan,{"cln":0, "gch":0}) for vchan in vchans )
    num_gch = 0.0
    num_cln = 0.0
    for ind in xrange(len(output['vchan'])):
        vchan = output['vchan'][ind]
        i = output['i'][ind]

        performance[vchan]['gch'] += i
        performance[vchan]['cln'] += 1-i

        num_gch += i
        num_cln += 1-i

    return performance, num_gch, num_cln

def file_to_perf( filename ):
    perf = {}
    ngch = 0
    ncln = 0
    file_obj = open(filename, "r")
    for line in file_obj:
        if line[0] != "#":
            chan, gch, ngch, cln, ncln, eff, fap = line.strip().split()
            perf[chan] = {'gch':float(gch), 'cln':float(cln)}

    return perf, float(ngch), float(ncln)

def perf_to_file( performance, num_gch, num_cln, filename ):
    vchans = sorted(performance.keys() )
    file_obj = open(filename, "w")
    print >> file_obj, "# %-48s %6s %6s %6s %5s %10s %10s"%("channel", "No.gch", "totgch", "No.cln", "totcln", "eff", "fap")
    for vchan in vchans:
        print >> file_obj, "%-50s %6d %6d %6d %6d %0.8f %0.8f"%(vchan, performance[vchan]['gch'], num_gch, performance[vchan]['cln'], num_cln, performance[vchan]['gch']/num_gch, performance[vchan]['cln']/num_cln)
    file_obj.close()


#=================================================
# mapping objects for r --> fraction of events
#=================================================
def binomialUL( k, n , conf=0.99, jefferys=True):
    """
    returns the binomial UL consistent with "conf" and observing k successes after n trials
    delegates to scipy.stats.beta, because we assume a jeffery's prior for the binomial process

    if we apply a prior: Beta(0.5, 0.5)

    the posterior for is : p ~ Beta( n*p + 0.5, n*(1-p) + 0.5 )

    where p = k/n (the MLE estimate for the binomial success rate)
    """
    from scipy.stats import beta
    if jefferys:
        ans = beta.ppf( conf, k + 0.5, n-k + 0.5 ) ### probably the fastest implementation available...
    else:
        ans = beta.ppf( conf, k, n-k )
    ans[k==n] = 1 ### get rid of Nan returned by beta.ppf
    return ans

def binomialCR( k, n, conf=0.99, jefferys=True):
    """
    returns the symetric CR consistent with "conf" and observing k successes after n trials
    delegates to scipy.stats.beta, because we assume a jeffery's prior for the binomial process

    if we apply a prior: Beta(0.5, 0.5)

    the posterior for is : p ~ Beta( n*p + 0.5, n*(1-p) + 0.5 )

    where p = k/n (the MLE estimate for the binomial success rate)
    """
    from scipy.stats import beta
    alpha = (1.0-conf)/2

    if jefferys:
        return beta.ppf( alpha, k + 0.5, n-k + 0.5) , beta.ppf( 1-alpha, k+0.5, n-k+0.5)
    else:
        return beta.ppf( alpha, k, n-k), beta.ppf( 1-alpha, k, n-k)

class BinomialMap(object):
    """
    a class that can stores the mapping between rank --> fraction of events.
    does this in a callable way
    can also return error bars on the point estimate (quantiles, how do we want to do this?)

    sorts ranks and stores them in that order, along with c_samples
    """

    def __init__(self, ranks, c_samples, kind='linear'):
        if len(ranks) != len(c_samples):
            raise ValueError("incompatible lenghts between ranks and c_samples")
        if not isinstance(ranks, numpy.ndarray):
            ranks = numpy.array( ranks )
        if not isinstance(c_samples, numpy.ndarray):
            c_samples = numpy.array( c_samples )

        argorder = ranks.argsort() ### make sure ranks are sorted in ascending order

        self.ranks = ranks[argorder]
        self.c_samples = c_samples[argorder] 

        self.n_samples = numpy.max(c_samples) ### the biggest number of cumulative samples is the total number of samples
        self.kind = kind

    def __call__(self, r):
        """
        returns the point estimate for the fraction of events
        """
        if not self.n_samples:
            return numpy.ones_like(r) ### no samples, so we have a big FAP

        if self.kind == "linear":
            return numpy.interp(r, self.ranks, self.c_samples) / self.n_samples
        else:
            raise ValueError("do not know how to interpolate with kind=%s"%kind)

    def ul(self, r, conf=0.99):
        """
        returns the upper-limit on the fraction of events corresponding to confidence="conf"
        delegates to binomialUL
        """
        p = self(r) ### get point estimate
        return binomialUL( self.n_samples*p , self.n_samples*(1-p), conf=conf)

    def cr(self, r, conf=0.99):
        """
        returns the symmetric conficence region on the fraction of events corresponding to confidence="conf"
        delegates to binomialCR
        """ 
        p = self(r)
        return binomialCR( self.n_samples*p , self.n_samples*(1-p), conf=conf)

#=================================================
# segment queries
#=================================================
def segment_query( config, start, end, url=None ):
    """
    Performes ligolw_segment query for segments in the given gps time range and returns segment xml file.
    """

#   from glue.ligolw import table,lsctables,utils
#   ligolw_segment_query -t file:///gds-l1/dmt/triggers/DQ_Segments -q --include-segments=L1:DMT-SCIENCE -s 1042257690 -e 1042257715

    program = config.get('get_science_segments', 'program')
    if url is None:
        url = config.get('get_science_segments', 'xmlurl')
    include = config.get('get_science_segments', 'include')
    args = '%s -t %s -q -a %s -s %d -e %d' % (program, url, include,
            start, end)

    print 'segment query: ', args

    # query to database is faster
    # need to set export S6_SEGMENT_SERVER=https://segdb-er.ligo.caltech.edu
    # and then run with -d option instead of -t
    # args....= '%s -d -q -a %s -s %d -e %d' % (program, include, start, end)

    xmlfileobj = tempfile.TemporaryFile()
    p = subprocess.Popen(args.split(), stdout=xmlfileobj) #, stderr=sys.stderr)
    p.wait()
    xmlfileobj.seek(0)  # go back to beginning of file

    return xmlfileobj

def get_idq_segments( realtimedir, start, stop, suffix='.pat'):

    fileseg = event.fixsegments( [extract_start_stop(filename, suffix=suffix) for filename in get_all_files_in_range(realtimedir, start, stop, suffix=suffix)] )
    return event.andsegments([[start,stop]], fileseg)

#    matchfile = re.compile('.*-([0-9]*)-([0-9]*)\%s$' % suffix)
#
#    # (dirname, starttime, endtime, pad=0, suffix='.xml')
#    ls = get_all_files_in_range(realtimedir, start, stop, suffix=suffix)
#    fileseg = []
#    for file in ls:
#        m = matchfile.match(file)
#        (st, dur) = (int(m.group(1)), int(m.group(2)))
#        fileseg.append([st, st + dur])
#    fileseg = event.fixsegments(fileseg)
#
#    return event.andsegments([[start, stop]], fileseg)

#=================================================
# data discovery
#=================================================

def get_all_files_in_range(dirname, starttime, endtime, pad=0, suffix='.xml' ):
    """
    Returns all files in dirname and all its subdirectories whose names indicate that they contain segments in the range starttime to endtime
    """

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

def extract_dmt_segments(xmlfiles, dq_name):
    """
....Loads the segments from the dmt xml files. Determines segments covered by dmt process from the file names
....dq_name is the segment definer name.
...."""

    good = []
    covered = []
    lsctables.use_in(ligolw.LIGOLWContentHandler)
    for file in xmlfiles:
        xmldoc = ligolw_utils.load_filename(file, contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))  # load file
        # get segment tables

        sdef = table.get_table(xmldoc,
                               lsctables.SegmentDefTable.tableName)
        seg = table.get_table(xmldoc, lsctables.SegmentTable.tableName)

        # segment definer ID correspodning to dq_name
        # FIX ME: in case of multiple segment versions this matching becomes ambiguous

        dq_seg_def_id = next(a.segment_def_id for a in sdef if a.name
                             == dq_name.split(":")[1])
                             #== dq_name)

        # get list of  segments

        good.extend([[a.start_time, a.end_time] for a in seg
                    if a.segment_def_id == dq_seg_def_id])

    # coalesce segments....
    good = event.fixsegments(good)

    # convert file names into segments

    #covered = map(filename_to_segment, xmlfiles)
    for xmlfile in xmlfiles:
        covered.append(extract_start_stop(xmlfile,suffix='.xml'))

    # coalesce segments

    covered = event.fixsegments(covered)

    return (good, covered)

def get_scisegs(segments_location, dq_name, start, end, pad=0, logger=None):
    """
    finds and returns dq segments in this range
    a helper function for retrieve_scisegs, which builds in waiting logic
    """
    try:
        dmtfiles = get_all_files_in_range(segments_location, start, end, pad=pad, suffix='.xml')
        (good, covered) = extract_dmt_segments(dmtfiles, dq_name)
    except Exception as e:
        if logger:
            logger.warning('error from DMT segment query: %s, %s, %s' % (e[0], e[1], e[2]))
        else:
            print 'error from DMT segment query: %s, %s, %s' % (e[0], e[1], e[2])
        (good, covered) = ([], [])

    return (good, covered)

def retrieve_scisegs(segments_location, dq_name, start, stride, pad=0, sleep=0, nretry=0, logger=None):
    """
    finds and returns dq segments in this range via delegation to get_scisegs.
    included logic about waiting and re-trying through optional arguments.
    will sleep for "sleep" seconds if there is not enough coverage, and will re-try "nretry" times.
    """
    if sleep < 0:
        sleep = 0
    end = start+stride

    ### get DMT xml segments (changing the query)
    good, covered = get_scisegs(segments_location, dq_name, start, end, pad=pad)

    ntrial = 0
    while (event.livetime(covered) < stride) and ntrial < nretry:
        if logger:
            logger.info(' unknown science coverage, waiting additional %d seconds' % sleep)
        else:
            print ' unknown science coverage, waiting additional %d seconds' % sleep

        time.sleep(sleep)
        if sleep < stride: ### automatically increase sleep to 1 stride if needed
            sleep = stride

        good, covered = get_scisegs(segments_location, dq_name, start, end, pad=pad)

    return good, covered

def retrieve_kwtrig(gdsdir, kwbasename, t, stride, sleep=0, ntrials=1, logger=None, delay=0, verbose=True):
    """
    looks for kwtriggers and includes logic about waiting.
    will wait "sleep" seconds between each check that the file has appeared. Will check "ntrials" times
    """
    if ntrials < 1:
        raise ValueError("ntrials must be >= 1, otherwise we don't do anything....")
    if sleep < 0:
        sleep = 0

    ### check if KW file is there
    kwfilename = '%s/%s/%s-%d/%s-%d-%d.trg' % ( gdsdir, kwbasename, kwbasename, t / 1e5, kwbasename, t, stride )

    for i in xrange(ntrials):
        if not os.path.exists(kwfilename):  ### kw files may not have appeared yet
            if verbose:
                if logger:
                    logger.info('  missing KW triggers, waiting additional %d seconds' % sleep)
                else:
                    print '  missing KW triggers, waiting additional %d seconds' % sleep
            if i < ntrials-1:
                time.sleep(sleep)
                if sleep < stride: ### automatically increase sleep to 1 stride if needed
                    sleep = stride
        else:
            break
    else:
        if verbose:
            if logger:
                logger.warning('  still missing KW triggers, skipping')
            else:
                print '  still missing KW triggers, skipping'
        return None

    if verbose:
        if logger:
            logger.info('  loading KW triggers : %s'%kwfilename)
        else:
            print '  loading KW triggers : %s'%kwfilename

    if delay > 0:
        if verbose:
            if logger:
                logger.info('    waiting %.3f seconds to ensure file is completely written'%delay)
            else:
                print '    waiting %.3f seconds to ensure file is completely written'%delay
        time.sleep(delay)

    return event.loadkwm(kwfilename)

def retrieve_kwtrigs(gdsdir, kwbasename, t, stride, kwstride, sleep=0, ntrials=1, logger=None, segments=None, verbose=True):
    t_start = (int(t) / kwstride) * kwstride

    trgdicts = []
    while t_start < t+stride:
        if (segments==None) or (event.livetime(event.andsegments([[[t_start, t_start+kwstride]], segments]))): ### only keep if we've got some overlap
            trigger_dict = retrieve_kwtrig(gdsdir, kwbasename, t_start, kwstride, sleep=sleep, ntrials=ntrials, logger=logger, verbose=verbose)
            if trigger_dict:
                trgdicts.append( trigger_dict )
        t_start += kwstride

    if trgdicts:
        trigger_dict = trgdicts[0]
        for td in trgdicts[1:]:
            trigger_dict.add( td )
        return trigger_dict
    else:
        return event.trigdict()
#        return None

def retrieve_DMTOmegaTrig(gdsdir, t, stride, channels, sleep=0, ntrials=1, logger=None, delay=0, verbose=True):
    """
    looks for DMT-Omega triggers and includes logic about waiting.
    will wait "sleep" seconds between each check that the file has appeared. Will check "ntrials" time
    """
    if ntrials<1:
        raise ValueError("ntrials must be >= 1, otherwise we don't do anything....")
    if sleep < 0:
        sleep = 0

    trgdict = event.trigdict()
    for channel in channels:

        filename = "%s/%d/%s-%d-%d.xml"%(gdsdir, t / 1e5, channel, t, stride)

        for i in xrange(ntrials):
            if not os.path.exists(filename):  ### kw files may not have appeared yet
                if verbose:
                    if logger:
                        logger.info('  missing DMT-Omega triggers for %s, waiting additional %d seconds' % (channel, sleep) )
                    else:
                        print '  missing DMT-Omega triggers for %s, waiting additional %d seconds' % (channel, sleep)
                if i < ntrials-1:
                    time.sleep(sleep)
                    if sleep < stride: ### automatically increase sleep to 1 stride if needed
                        sleep = stride
            else:
                break
        else:
            if verbose:
                if logger:
                    logger.warning('  still missing DMT-Omega triggers for %s, skipping' % channel )
                else:
                    print '  still missing DMT-Omega triggers for %s, skipping' % channel
            continue

        if verbose:
            if logger:
                logger.info('  loading DMT-Omega triggers for %s : %s'%(channel, filename) )
            else:
                print '  loading DMT-Omega triggers for %s : %s'%(channel, filename)

        if delay > 0:
            if verbose:
                if logger:
                    logger.info('    waiting %.3f seconds to ensure file is completely written'%delay)
                else:
                    print '    waiting %.3f seconds to ensure file is completely written'%delay
            time.sleep(delay)

        trgdict.add( event.loadSingleBurst( filename ) )

    return trgdict



def retrieve_DMTOmegaTrigs(gdsdir, basename, t, stride, dmtostride, channels, sleep=0, ntrials=1, logger=None, segments=None, verbose=True):
    """
    discover and read in files from DMT-Omega processes
    we require the argument "channels" because DMT-Omega triggers are stored in separate files for each file
      -> we can save significant I/O if we only ever load a subset of these
    also, the separate single-channel files may appear with different latencies and without a definitive list it is difficult to construct waiting logic
    """
    t_start = (int(t) / ostride) * ostride

    trgdicts = []
    while t_start < t+stride:
        if (segments==None) or (event.livetime(event.andsegments([[[t_start, t_start+dmtostride]], segments]))): ### only keep if we've got some overlap
            trigger_dict = retrieve_DMTOmegaTrig(gdsdir, basename, t_start, dmtostride, channels, sleep=sleep, ntrials=ntrials, logger=logger, verbose=verbose)
            if trigger_dict:
                trgdicts.append( trigger_dict )
        t_start += ostride

    if trgdicts:
        trigger_dict = trgdicts[0]
        for td in trgdicts[1:]:
            trigger_dict.add( td )
        return trigger_dict
    else:
        return event.trigdict()


def retrieve_OmicronTrig(gdsdir, ifo, t, stride, channels, sleep=0, ntrials=1, logger=None, delay=0, verbose=True):
    """
    looks for Omicron triggers and includes logic about wiating.
    will wait "sleep" seconds between each check that the file has appeared. Will check "ntrials" times
    """

    if ntrials < 1:
        raise ValueError("ntrials must be >= 1, otherwise we don't do anything....")
    if sleep < 0:
        sleep = 0

    trgdict = event.trigdict()
    for channel in channels:
        ### check if Omicron file is there
        # /home/reed.essick/Omicron/test/triggers/H-11242/H1:GDS-CALIB_STRAIN/H1-GDS_CALIB_STRAIN_Omicron-1124203561-30.xml

        ### fancy channel name used to predict filename directory structure...
        fancy_channel = channel.split("-")[-1].split("_")[:-1] ### get rid of ifo and trailing _Omicron
        fancy_channel = "%s-%s"%(fancy_channel[0], "_".join(fancy_channel[1:])) ### only "-" is after the subsys (for all channels?

        filename = '%s/%s-%d/%s:%s/%s-%d-%d.xml' % ( gdsdir, ifo[0], t / 1e5, ifo, fancy_channel, channel, t+1, stride ) ### +1 is due to the way Omicron pads segments... this is fragile but should be relatively safe for the online Omicron processes managed by R. Essick (reed.essick@ligo.org)
#        filename = '%s/%s-%d/%s:%s/%s-%d-%d.xml' % ( gdsdir, ifo[0], t / 1e5, ifo, fancy_channel, channel, t, stride )

        for i in xrange(ntrials):
            if not os.path.exists(filename):  ### omicron files may not have appeared yet
                if verbose:
                    if logger:
                        logger.info('  missing Omicron triggers for %s, waiting additional %d seconds' % (channel, sleep) )
                    else:
                        print '  missing Omicron triggers for %s, waiting additional %d seconds' % (channel, sleep)
                if i < ntrials-1:
                    time.sleep(sleep)
                    if sleep < stride: ### automatically increase sleep to 1 stride if needed
                        sleep = stride
            else:
                break
        else:
            if verbose:
                if logger:
                    logger.warning('  still missing Omicront triggers for %s, skipping' % channel )
                else:
                    print '  still missing Omicron triggers for %s, skipping' % channel
            continue

        if verbose:
            if logger:
                logger.info('  loading Omicron triggers for %s : %s'%(channel, filename) )
            else:
                print '  loading Omicron triggers for %s : %s'%(channel, filename)

        if delay > 0:
            if verbose:
                if logger:
                    logger.info('    waiting %.3f seconds to ensure file is completely written'%delay)
                else:
                    print '    waiting %.3f seconds to ensure file is completely written'%delay
            time.sleep(delay)

        trgdict.add( event.loadSingleBurst( filename ) )

    return trgdict

def retrieve_OmicronTrigs(gdsdir, ifo, t, stride, ostride, channels, sleep=0, ntrials=1, logger=None, segments=None, verbose=True):
    """
    discover and read in files from "online" omicron processes
    we require the argument "channels" because Omicron triggers are stored in separate files for each file
      -> we can save significant I/O if we only ever load a subset of these
    also, the separate single-channel files may appear with different latencies and without a definitive list it is difficult to construct waiting logic
    """
    t_start = (int(t) / ostride) * ostride
   
    trgdicts = []
    while t_start < t+stride:
        if (segments==None) or (event.livetime(event.andsegments([[[t_start, t_start+ostride]], segments]))): ### only keep if we've got some overlap
            trigger_dict = retrieve_OmicronTrig(gdsdir, ifo, t_start, ostride, channels, sleep=sleep, ntrials=ntrials, logger=logger, verbose=verbose)
            if trigger_dict:
                trgdicts.append( trigger_dict )
        t_start += ostride

    if trgdicts:
        trigger_dict = trgdicts[0]
        for td in trgdicts[1:]:
            trigger_dict.add( td )
        return trigger_dict
    else:
        return event.trigdict()

def retrieve_OfflineOmicronTrigs(gdsdir, ifo, t, stride, channels=None, logger=None, segments=None, verbose=True):
    """
    discover and read in files from "standard" /home/detchar/triggers directory structure
    we allow the keyword argument "channels" because OfflienOmicron triggers are stored in separate files for each file
      -> we can save significant I/O if we only ever load a subset of these
    if we request a specific set of channels, this method will only load the necessary files and will ensure that trgdict has a key,value pair for every requested channel, even if no triggers are found.
    """
    xmlfiles = defaultdict( list )
    for xmlfile in get_all_files_in_range("%s/%s/"%(gdsdir, ifo), t, t+stride, suffix=".xml.gz"):
        if (segments==None) or (event.livetime(event.andsegments([ extract_start_stop(xmlfile, suffix=".xml.gz"), segments ] ))):
            xmlfiles[ extract_OfflineOmicron_name( xmlfile ) ].append( xmlfile )

    if channels!=None: ### only keep the relevant channels
        for key in xmlfiles.keys():
            if key not in channels:
                xmlfiles.pop( key )

    trgdict = event.trigdict()    
    for chan in sorted(xmlfiles.keys()):
        if verbose:
            if logger:
                logger.info('  loading OfflineOmicron triggers for : %s'%chan)
            else:
                print '  loading OfflineOmicron triggers for : %s'%chan
        for xmlfile in sorted(xmlfiles[chan]):
            if verbose:
                if logger:
                    logger.info("    loading Omicron triggers : %s"%xmlfile)
                else:
                    print "    loading Omicron triggers : %s"%xmlfile
            trgdict.add( event.loadSingleBurst( xmlfile ) )

    if channels!=None: ### ensure that we have a key,value pair for every requested channel
        for chan in channels:
            if not trgdict.has_key(chan):
                trgdict[chan] = []

    return trgdict


def extract_Omicron_name( filename ):
    return "-".join(filename.split("/")[-1].split("-")[:1])

#=================================================
# trigger generation algorithm interfaces
#=================================================

def loadkwconfig(file):
    """
    load in KW configuration file as dictionary
    note that nothing will be converted to numerical values, they will all be saved as string values
    also produce list of channel_flow_fhigh names
    original order of channels is preserved
    """
    table = event.loadstringtable(file)
    kwconfig = dict()
    singleopt = 'stride basename segname threshold significance maxDist decimateFactor transientDuration'.split()
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

def collect_sngl_chan_kw(
    gpsstart,
    gpsstop,
    kw_config,
    width=False,
    source_dir='./',
    output_dir='./',
    verbose=False,
    chans=None,
    scisegs=None,
    include_all_found_channels=False
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
            channels.append([fields[1], "%s_%s_%d_%d"%(fields[1][:2], fields[1][3:], float(fields[2]), float(fields[3])) ])
    f.close()
    
    if chans: ### if supplied, only keep those channels that are specified
        channels = [chan for chan in channels if chan[1] in chans]

    if not stride:
        print 'stride not defined in ' + config.kwconfig
        sys.exit(1)
    if not basename:
        print 'basename not defined in ' + config.kwconfig
        sys.exit(1)

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
            t_S = int(t + 0.0000001)  # make sure we don't have a rounding error?
            t_dir = t_S / 100000  # the digits in the GPS time beyond 1e5, used for structure of KW output
            file = source_dir + '/' + basename + '-' + repr(t_dir) \
                + '/' + basename + '-' + repr(t_S) + '-' \
                + repr(stride_S) + '.trg'

            ### check sciseg overlap
            if scisegs!=None and (not event.livetime( event.andsegments([scisegs, [[t, t+stride]]]) ) ): ### check for overlap with scisegs
                if verbose:
                    print "%s has no overlap with scisegs. Skipping"%file
                t += stride
                continue

            # build segment list for this day
            if os.path.exists(file):
                if len(segments) == 0:  # first segment
                    segments.append([t, t + stride])
                    lastt = t
                elif segments[-1][1] == t: # continuous data segment already generated
                    segments[-1][1] = t + stride
                    lastt = t
                else: # discontinuous data, we skip one section of data to eschew filter transients
                    if lastt + stride == t:  # not continuous with anything in segments, but we've already skipped one section of data
                        if verbose:
                            print 'forming new segment at ' + `t`
                        segments.append([t, t + stride])
                    else: # skip the first section of data in a new segment because of 'filter transients'
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
                        if line[0]=="#": ### respect comments
                            continue

                        fields = line.strip().split()
                        if len(fields) == 9 and fields[-1] != '':
                            tag = fields[-1]
                            if triggers.has_key(tag):
                                triggers[tag].append(line.strip()[:-len(tag)
                                    - 2])
                            elif include_all_found_channels:
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


#=================================================
# auxmvc helper functions
#=================================================
def build_auxmvc_vectors( trigger_dict, main_channel, time_window, signif_threshold, output_file_name, gps_start_time, gps_end_time, channels=None, unsafe_channels=None, science_segments=None, clean_times=None, clean_samples_rate=None, clean_window=None, filter_out_unclean=False, max_clean_samples=None, max_glitch_samples=None, verbose=False):
    """
    Given dictionary of triggers from multiple channels, the function constructs auxmvc
    vectors out of them. Result is saved into output_file_name. 
    """

    if not trigger_dict: # empty trig-dictionary
        raise StandardError('Empty trig-dictionary. Can not build auxmvc vectors.')

    if verbose:
        print 'Number of triggers in the main channel before thresholding:', len(trigger_dict[main_channel])

    # use only channels from the channels file, if provided
    if channels:
        selected_channels = event.read_channels_from_file(channels)
        if main_channel not in selected_channels:
            selected_channels.append( main_channel )
        trigger_dict.keep_channels(selected_channels)

        # to ensure consistency in dimensionality of auxmvc vectors
        # add the channels from the selected channels that are absent in trigger_dict
        for channel in selected_channels:
            if not channel in trigger_dict.channels():
                trigger_dict[channel] = []

    # get rid of unsafe channels if provided
    if unsafe_channels:
        unsafe_channels = event.read_channels_from_file(unsafe_channels)
        trigger_dict.remove_channels([channel for channel in unsafe_channels if channel != main_channel]) ### never remove the main channel

    # keep only the triggers from the [gps_start_time, gps_end_time] segment
    # first keep all triggers from the segment expanded by the time concidence window, so that not to loose coincidences
    trigger_dict.include([[gps_start_time - time_window, gps_end_time + time_window]])

    # then in the main channel keep triggers that fall within the segment.
    trigger_dict.include([[gps_start_time, gps_end_time]], channels=[main_channel])

    # keep only triggers from the science segments if given........
    if science_segments!=None:
        science_segments = event.andsegments([[gps_start_time - time_window, gps_end_time + time_window]], science_segments)
        trigger_dict.include(science_segments)

    # apply significance threshold to the triggers from the main channel
    trigger_dict.apply_signif_threshold(channels=[main_channel], threshold=signif_threshold)
    print 'Number of triggers in the main channel after thresholding:', len(trigger_dict[main_channel])

    # construct glitch auxmvc vectors
    aux_glitch_vecs = event.build_auxmvc_vectors(trigger_dict, main_channel=main_channel, coincidence_time_window=time_window)

    # apply upper limit on the number of glitch samples if given
    if max_glitch_samples:
        if len(aux_glitch_vecs) > max_glitch_samples:
            aux_glitch_vecs = aux_glitch_vecs[-max_glitch_samples:]

    if not clean_times:
        if clean_samples_rate:

            # generate random times for clean samples
            if science_segments!=None:
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
    aux_clean_vecs = event.build_auxmvc_vectors( trigger_dict, main_channel=main_channel, coincidence_time_window=time_window, build_clean_samples=True, clean_times=clean_times, clean_window=clean_window )

    # get rid of clean samples that are near real triggers in the main channel.
    if filter_out_unclean:
        aux_clean_vecs = auxmvc_utils.get_clean_samples(aux_clean_vecs)

    # apply upper limit on the number of clean samples if given
    if max_clean_samples:
        if len(aux_clean_vecs) > max_clean_samples:
            aux_clean_vecs = aux_clean_vecs[-max_clean_samples:]

    # convert glitch and clean auxmvc vectors into MVSC evaluation set
    mvsc_evaluation_set = \
        auxmvc_utils.ConvertKWAuxToMVSC(KWAuxGlitchTriggers=aux_glitch_vecs, KWAuxCleanTriggers=aux_clean_vecs)

    # save MVSC evaluation set in file........
    auxmvc_utils.WriteMVSCTriggers(mvsc_evaluation_set, output_filename=output_file_name, Classified=False)

    return mvsc_evaluation_set

#def execute_build_auxmvc_vectors(
#    cp,
#    execute_dir,
#    trigdir,
#    main_channel,
#    output_file,
#    gps_start_time,
#    gps_end_time,
#    channels=None,
#    unsafe_channels=None,
#    dq_segments=None,
#    dq_segments_name='',
#    ):
#    """
#....Submits the job that builds auxmvc feature vectors. Vectors are saved in the output file with .pat extension.
#....Waits until the job is finished, returns its exit status and the output file name.
#....execute_dir is the absolute path of the directory in which the job should be executed. 
#...."""
#
#    # initiate build_auxmvc_vectors job object
#    build_auxmvc_vectors_job = auxmvc.build_auxmvc_vectors_job(cp,
#            main_channel, channels=channels,
#            unsafe_channels=unsafe_channels)
#
#    # create node for this job
#    build_auxmvc_vectors_node = auxmvc.build_auxmvc_vectors_node(
#        build_auxmvc_vectors_job,
#        trigdir,
#        gps_start_time,
#        gps_end_time,
#        output_file,
#        dq_segments=dq_segments,
#        dq_segments_name=dq_segments_name,
#        )
#
#    # get full command line for this job
#    build_auxmvc_vectors_command = auxmvc.construct_command(build_auxmvc_vectors_node)
#
#    print " ".join(build_auxmvc_vectors_command)
#
#    # submit process
##    exit_status = submit_command(build_auxmvc_vectors_command, 'build_auxmvc_vectors', execute_dir, verbose=True)
#    exit_status = subprocess.Popen(build_auxmvc_vectors_command, cwd=execute_dir).wait() ### block!
#
#    return (exit_status, build_auxmvc_vectors_node.get_output_files()[0])

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
    prepare_training_auxmvc_samples_job = auxmvc.prepare_training_auxmvc_samples_job(cp)

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
    prepare_training_auxmvc_samples_command =  auxmvc.construct_command(prepare_training_auxmvc_samples_node)

    # submit process
#    exit_status = submit_command(prepare_training_auxmvc_samples_command, 'prepare_training_auxmvc_samples', execute_dir, verbose=True)
    exit_status = subprocess.Popen(prepare_training_auxmvc_samples_command, cwd=execute_dir).wait() ### block!

    return (exit_status, prepare_training_auxmvc_samples_node.get_output_files()[0])



#===================================================================================================
#
# EXECUTION
#
#===================================================================================================
def evaluate(flavor, lines, dat, config, gps_start_time=-numpy.infty, gps_end_time=numpy.infty, dir=".", samples=[], trgdict=None, auxmvc_vectors=None, samples_header=['GPS','i','rank']):
    """
    a delegation function that calls the correct functions based on flavor
    """
    Ntotal = len(samples)

    if flavor == "ovl":
        vetolist = lines[0]
        if len(samples):
            (gw_predict, dat) = ovl_evaluate(vetolist, GPStimes=samples, GPS_headers=samples_header, allvtrg=trgdict, filename=dat.split("/")[-1], output_dir=dir )
            return 0
        else:
            ovl_dat_vars = ['rank', 'eff/dt', 'vchan', 'vthr', 'vwin']
            varline = ' '.join(samples_header) + ' ' + ' '.join(ovl_dat_vars) + '\n'
            file = open(dat, 'w')
            print >>file, "%s %s"%(' '.join(samples_header), ' '.join(ovl_dat_vars))
            file.close()
            return 0

    elif flavor == "forest":
        pat = trgdict ### silly way I'm passing this
        if len(auxmvc_vectors):
            trainedforest = lines[0]
            return forest_evaluate( pat, trainedforest, dat, config, gps_start_time=gps_start_time, gps_end_time=gps_end_time, dir=dir )[0]
        else:
            nonaux_vars = ['index', 'i', 'w', 'GPS_s', 'GPS_ms', 'signif', 'SNR', 'unclean', 'glitch-rank' ]
            file = open(dat, 'w')
            print >> file, 'GPS i w unclean signif SNR rank %s' % (' '.join([var for var in auxmvc_vectors.dtype.names if not var in nonaux_vars]))
            file.close()
            return 0

    elif flavor == "svm":
        pat = trgdict ### silly way I'm passing this
        if len(auxmvc_vectors):
            svm_model = lines[0].strip('\n')
            svm_range_file = lines[1].strip('\n')

            return svm_evaluate( config, pat, svm_range_file, svm_model, dat, dir=dir )
        else:
            nonaux_vars = ['index', 'i', 'w', 'GPS_s', 'GPS_ms', 'signif', 'SNR', 'unclean', 'glitch-rank' ]
            file = open(dat, 'w')
            print >> file, 'GPS i w unclean signif SNR rank %s' % (' '.join([var for var in auxmvc_vectors.dtype.names if not var in nonaux_vars]))
            file.close()
            return 0

    elif flavor == "ann":
        pat = trgdict ### silly way I'm passing this
        if len(auxmvc_vectors):
            trained_ann = lines[0]
            return ann_evaluate( pat, trained_ann, dat, config, gps_start_time=gps_start_time, gps_end_time=gps_end_time, dir=dir )[0]
        else:
            nonaux_vars = ['index', 'i', 'w', 'GPS_s', 'GPS_ms', 'signif', 'SNR', 'unclean', 'glitch-rank' ]
            file = open(dat, 'w')
            print >> file, 'GPS i w unclean signif SNR rank %s' % (' '.join([var for var in auxmvc_vectors.dtype.names if not var in nonaux_vars]))
            file.close()
            return 0

    else:
        raise ValueError("do not know how to execute for flavor=%s"%flavor)

#===================================================================================================
#
# TRAINING
#
#===================================================================================================
def dag_train(flavor, pat, cache, config, train_dir='.', cwd='.'):
    """
    a delegation funtion that call the correct functions based on flavor
    """
    if flavor == "forest":
        ### submit mvsc training job
        (submit_dag_exit_status, dag_file, trained_file) = execute_forest_train(pat, cache.name, config, train_dir)
        os.chdir(cwd) ### switch back to the current directory

        return submit_dag_exit_status, "%s/%s"%(train_dir, dag_file)

    elif flavor == "svm":
        ### auxiliary files required for SVM training
        svm_range_file = "%s/%s.range"%(train_dir, os.path.split(pat)[1])
        svm_model_file = "%s/%s.model"%(train_dir, os.path.split(pat)[1])

        ### submit SVM training job
        (submit_dag_exit_status, dag_file, output_files) = execute_svm_train( config, pat, svm_range_file, svm_model_file, cache.name, train_dir )
        os.chdir(cwd) ### switch back to the current directory

        return submit_dag_exit_status, "%s/%s"%(train_dir, dag_file)

    elif flavor == "ann":
        ### submit ann training job
        (submit_dag_exit_status, dag_file, trained_file) = execute_ann_train(pat, cache.name, config, train_dir)
        os.chdir(cwd) ### switch back to the current directory

        return submit_dag_exit_status, "%s/%s"%(train_dir, dag_file)


    else:
        raise ValueError("do not know how to dag_train flavor=%s"%flavor)

def blk_train(flavor, config, classifierD, start, stop, ovlsegs=False, vetosegs=False, train_dir='.', cwd='.', force=False, cache=None, min_num_gch=0, min_num_cln=0, padding=0, conn=False):
    """
    a delegation function that calls the correct function based on flavor
    """
    if flavor == "ovl":
        if ovlsegs: ### silly formatting thing
            ovlsegs = [ovlsegs]

        vetolists = ovl_train( start, stop, dict(config.items('general')), classifierD, scisegs=ovlsegs, vetosegs=vetosegs, output_dir=train_dir, padding=padding)
        vetolist = vetolists[0].replace("//","/")### make sure the path doesn't contain any "//" elements

        if cache: ### we're told where we might append
            ### append new vetolist to training cache
            if not force:
                go_ahead_and_append = False
                f = open(vetolist, 'r')
                for line in f:
                    if line[0] != '#': ### the first un-commented line will contain the total number of gwtriggers
                                   # require certain number of glitches for this to be a good training set
                        go_ahead_and_append = float(line.strip().split()[ovl.vD['#gwtrg']]) >= min_num_gch
                        break
                f.close()

            if force or go_ahead_and_append:
                cache.append( vetolist )

                if conn:
                    conn.send((0, vetolist))
                else:
                    return 0, vetolist ### there should be only one list...
            else:
                if conn:
                    conn.send((1, vetolist))
                else:
                    return 1, vetolist

    else:
        raise ValueError("do not know how to blk_train flavor=%s"%flavor)

#=================================================
# timeseries
#=================================================

def timeseries2frame( framename, chan_vect, start, dt):
    commonpar = {'start':start, 'dx':dt, 'type':1}
    outdict = [dict(name=channel_name, data=x, **commonpar) for channel_name, x in chan_vect.items()]

    Fr.frputvect(framename, outdict)

def datfile2timeseries(flavor, dat, lines, trgdict=None, start=0, stride=0, window=0.1, fs=256):
    """
    create veto timeseries from datefile glitch events
    will create square window timeseries with +/-window around each event
    """
    if flavor == "ovl":
        return ovl.vetolist2timeseries(lines[0], trgdict, [start, start + stride], fs=fs)

    elif flavor in mla_flavors:
        datdata = slim_load_datfile(dat, skip_lines=0, columns='GPS i rank'.split())
        for (key, val) in datdata.items():
            datdata[key] = map(float, val)

#        (ifo, name, tstart, dur) = os.path.basename(dat)[:-4].split('-')
#        livetime = float(dur)
#        segment = [float(tstart), float(tstart) + livetime]
        tstart, tstop = extract_start_stop(dat, suffix=".dat")
        segment = [tstart, tstop]
        livetime = tstop-tstart

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

    else:
        raise ValueError("don't know how to compute timeseries for flavor=%s"%flavor)

def timeseries_to_segments(t, ts, thr):
    """
    computes segments from t  = time stamps
                            ts = time series (values)
                            thr=threshold on time series
    so that t \f$\in\f$ segments iff ts(t) >= thr
    """    
    truth = ts >= thr ### determine which time samples are above the threshold
    if numpy.any(truth):
        edges = list(numpy.nonzero(truth[1:]-truth[:-1])[0]+1) ### find the edges corresponding to state changes
        if truth[0] and (edges[0]!=0): ### we start off in a segment and that edge is not included
            edges.insert(0, 0)
        if truth[-1] and (edges[-1]!=len(ts)-1): ### we end in a segment and that edge is not included
            edges.append( len(ts)-1 )

        if len(edges)%2:
            raise ValueError("odd number of edges...something is wrong")

        edges = numpy.array(edges)
        segs = numpy.transpose( numpy.array([t[edges[:-1:2]], t[edges[1::2]]]) )
        return [list(seg) for seg in segs], numpy.min(ts[truth])
    else:
        return [], None

    '''
    segs = []
    in_seg = False
    min_TS = numpy.infty
    for (T, TS) in zip(t, ts):
        if TS < thr:
            if in_seg:
                segs.append([start, T])
                in_seg = False
            else:
                pass
        else: ### TS >= thr

            if min_TS > TS:
                min_TS = TS
            if not in_seg:
                start = T
                in_seg = True
    if in_seg:
        segs.append([start, T])

    if len(segs):
        return (segs, min_TS)
    else:
        return ([], None)
    ''' 

def combine_gwf(filenames, channels):
    """
    combine mutliple frame files into a single time-series. Assumes filenames have the standard LIGO naming convention : *-start-dur.suffix
    Also assumes that filenames are sorted into chronological order

    returns a list of arrays, with each array consisting of only contiguous data
        return timeseries, times

    channels is a list of channel names that are to be read from the frames.
    channels will be stored in timeseries[i] in the order in which they appear in "channels"
    """

    n = len(channels)
    if n < 1:
        return numpy.array([]), numpy.array([])
    elif n > 1:
        ts = numpy.array([[]]*n)
    else:
        ts = numpy.array([[]])

    t = numpy.array([])

    times = []
    timeseries = []

    end = False
    for filename in filenames:
        _start, _end = extract_start_stop(filename, suffix='.gwf')
        _dur = _end-_start

        if not end or end==_start: # beginning of data
            end = _end
            
            _ts = []
            for ind, channel in enumerate(channels):
                _ts.append( Fr.frgetvect1d( filename, channel )[0] )
            _ts = numpy.array( _ts )

            ts = numpy.concatenate((ts, _ts), axis=1)
            len_ts = len(_ts[0])
            t = numpy.concatenate( (t, numpy.arange(_start, _end, 1.0*_dur/len_ts) ) )

        else:
            # gap in the data!
            times.append(t)  # put old continuous data into lists
            timeseries.append(ts)

            ts = []
            for ind, channel in enumerate(channels):
                ts.append( Fr.frgetvect1d( filename, channel )[0] )
            ts = numpy.array( ts )
            len_ts = len(ts[0])
            t = numpy.arange(_start, _end, 1.0*_dur/len_ts)
            end = _end

    times.append(t)
    timeseries.append(ts)

    if n > 1:
        return (times, timeseries)
    else:
        return (times, [ts[0] for ts in timeseries])

def combine_ts(filenames, n=1):
    """ 
    combine multiple files into a single time-series. Assumes filenames have the standard LIGO naming covention: *-start-dur.suffix
    Also assumes that filenames are sorted into chronological order

    returns lists of arrays, with each array consisting of only contiguous data
        return timeseries, times

    n=number of columns expected in the ts arrays -> len(ts) = n
    """
    t = numpy.array([])  # the array storing continuous data
    if n > 1:
        ts = numpy.array([[]]*n)
    else:
        ts = numpy.array([])
    times = []  # a list that will contain stretches of continuous data
    timeseries = []

#    matchfile = re.compile('.*-([0-9]*)-([0-9]*).*$')
    end = False
    for filename in filenames:
        _start, _end = extract_start_stop(filename, suffix='.npy.gz')
        _dur = _end-_start
#        m = matchfile.match(filename)
#        (_start, _dur) = (int(m.group(1)), int(m.group(2)))

        ### check to see if we have continuous data
        if not end or end == _start:  # beginning of data
            end = _start + _dur

            _file = event.gzopen(filename)
            _ts = numpy.load(_file)
            _file.close()

            ts = numpy.concatenate((ts, _ts), axis=1) ### axis argument takes care of multi-dimensional arrays
            if n > 1:
                len_ts = len(_ts[0])
            else:
                len_ts = len(_ts)
            t = numpy.concatenate((t, numpy.arange(_start, _start + _dur, 1.0* _dur / len_ts)))

        else:
            # gap in the data!
            times.append(t)  # put old continuous data into lists
            timeseries.append(ts)

            _file = event.gzopen(filename)  # start new continuous data
            ts = numpy.load(_file)
            _file.close()
            if n > 1:
                len_ts = len(ts[0])
            else:
                len_ts = len(ts)
            t = numpy.arange(_start, _start + _dur, 1.0 * _dur / len_ts)
            end = _start + _dur

    times.append(t)
    timeseries.append(ts)

    return (times, timeseries)

def stats_ts(ts):
    """ 
    compute basic statistics about ts 

    return min(ts), max(ts), mean(ts), stdv(ts)
    """

    return (numpy.min(ts), numpy.max(ts), numpy.mean(ts), numpy.std(ts))

def timeseries_in_segments(t, ts, segs):
    """
    returns only those elements of the time-series which are contained in the segments
    """
    if not len(t): ### no sample points...
        return t, ts

    truth = numpy.zeros_like(t).astype(bool)

    for start, end in segs:
        truth += (start <= t)*(t <= end) ### keep only times within this segment

    return t[truth], numpy.transpose( numpy.transpose(ts)[truth] ) ### fancy ts stuff will handle multiple rows

#=================================================
# idq's xmltables
#=================================================
def datfile2xmldocs(datfilename, ifo, FAPmap, Lmap=False, Effmap=False, flavor=None, gwchan = None, gwtrigs=None, prog = None, options = None, version = None ):
    """
    reads in a datfile and generates xml tables containing the equivalent information, but adds FAP and likelihood as computed by FAPmap and Lmap (expected to be callable)
        if Lmap is not supplied, attempts to use L = Eff/FAP as estimate with Effmap
        ==> must supply either Lmap or Effmap

    return gchxmldoc, clnxmldoc
    """

    if not (Lmap or Effmap):
        raise StandardError('must supply either Lmap or Effmap')

    gchxml_doc = ligolw.Document()  # generate the xml document object
    gch_ligolw_element = ligolw.LIGO_LW()  # generate the table-tree element
    gchxml_doc.appendChild(gch_ligolw_element)  # put table-tree into xml document

    clnxml_doc = ligolw.Document()  # repeat for clnTable
    cln_ligolw_element = ligolw.LIGO_LW()
    clnxml_doc.appendChild(cln_ligolw_element)

    # add process and process_params tables to xml docs
    gchproc = process.register_to_xmldoc(gchxml_doc, prog, options, version = version)
    clnproc = process.register_to_xmldoc(clnxml_doc, prog, options, version = version)

    gchTable = lsctables.New(idq_tables.IDQGlitchTable)  # instantiate glitch table objects
    clnTable = lsctables.New(idq_tables.IDQGlitchTable)

    gchCoincDefTable = lsctables.New(lsctables.CoincDefTable) # instantiate coinc_def table objects
    clnCoincDefTable = lsctables.New(lsctables.CoincDefTable)

    gchCoincTable = lsctables.New(lsctables.CoincTable) # instantiate coinc_event table objects
    clnCoincTable = lsctables.New(lsctables.CoincTable)

    gchCoincMapTable = lsctables.New(lsctables.CoincMapTable) # instantiate coinc_event_map table objects
    clnCoincMapTable = lsctables.New(lsctables.CoincMapTable)

    if gwtrigs:
        gchSnglBurstTable = lsctables.New(lsctables.SnglBurstTable) # instantiate sngl_burst table for gw triggers

    gch_ligolw_element.appendChild(gchTable)  # put gchTable into table-tree
    cln_ligolw_element.appendChild(clnTable)

    gch_ligolw_element.appendChild(gchCoincDefTable)  # put coinc_def table into table-tree
    cln_ligolw_element.appendChild(clnCoincDefTable)

    gch_ligolw_element.appendChild(gchCoincTable)  # put coinc_event table into table-tree
    cln_ligolw_element.appendChild(clnCoincTable)

    gch_ligolw_element.appendChild(gchCoincMapTable)  # put coinc_map table into table-tree
    cln_ligolw_element.appendChild(clnCoincMapTable)

    if gwtrigs:
        gch_ligolw_element.appendChild(gchSnglBurstTable)  # put sngl_burst table into table-tree

    # iterate over dat and place info into appropriate tables
    # load data and fill out xml tables
    if flavor == 'ovl':  # ovl records more info than the other classifiers
        # setup ovl specific tables
        gchOVLTable = lsctables.New(idq_tables.OVLDataTable)  # instantiate table objects
        clnOVLTable = lsctables.New(idq_tables.OVLDataTable)
        gch_ligolw_element.appendChild(gchOVLTable)  # add tables to document tree
        cln_ligolw_element.appendChild(clnOVLTable)

        # read in data
        dat = slim_load_datfile(datfilename, skip_lines=0, columns=['GPS', 'i', 'rank', 'vchan', 'vthr', 'vwin' ])
    else:
        dat = slim_load_datfile(datfilename, skip_lines=0, columns=['GPS', 'i', 'rank'])  # read in relevant information from datfile

    no_events = len(dat['GPS'])
    no_glitch_events = len([a for a in dat['i'] if float(a) == 1])
    no_clean_events = len([a for a in dat['i'] if float(a) == 0])

    #generate coinc_def tables
    if no_glitch_events > 0:
        if gwtrigs:
            gchCoincDefTable.get_coinc_def_id( search=idq_tables.IDQCoincDef['idq_glitch<-->sngl_burst'][0], search_coinc_type=idq_tables.IDQCoincDef['idq_glitch<-->sngl_burst'][1], create_new = True, description = 'idq_glitch<-->sngl_burst')
        if flavor == 'ovl':
            gchCoincDefTable.get_coinc_def_id( search=idq_tables.IDQCoincDef['idq_glitch<-->ovl_data'][0], search_coinc_type=idq_tables.IDQCoincDef['idq_glitch<-->ovl_data'][1], create_new = True, description = 'idq_glitch<-->ovl_data')

    if no_clean_events > 0:
        if flavor == 'ovl':
            clnCoincDefTable.get_coinc_def_id( search=idq_tables.IDQCoincDef['idq_glitch<-->ovl_data'][0], search_coinc_type=idq_tables.IDQCoincDef['idq_glitch<-->ovl_data'][1], create_new = True, description = 'idq_glitch<-->ovl_data')


    #WARNING: Below we associate glitches with GW triggers relying on the fact that both lists are time ordered and one-to-one.
    # It is OK as long as dat files are generated from the same GW triggers and order of elements was not changed.
    # check that the number of glitches in dat file is the same as number of gwtrigs
    if gwtrigs:
        if not len(gwtrigs) == no_glitch_events:
            print "WARNING: Length of gwtrigs does not match number of glitch events in the dat file."
            print "Glitch <-> GW trigger association will probably be incorrect"
    glitch_index = 0
    for ind in range(no_events):
        GPS = float(dat['GPS'][ind])
        i = float(dat['i'][ind])
        rank = float(dat['rank'][ind])

        ### add info to GlitchTable() objects. Commmon for all classifiers
        idq_glitch_row = idq_tables.IDQGlitch()  # define row object

        idq_glitch_row.ifo = ifo
        idq_glitch_row.gps = int(GPS)  # fill in relevant data
        idq_glitch_row.gps_ns = int( round( GPS % 1.0, ndigits=9 ) * 10**9 )
        idq_glitch_row.rank = float(rank)
        idq_glitch_row.fap = FAPmap(idq_glitch_row.rank)
        if Lmap:
            idq_glitch_row.likelihood = Lmap(idq_glitch_row.rank)
        else:
            if idq_glitch_row.fap == 0:
                if Effmap(idq_glitch_row.rank) == 0:
                    idq_glitch_row.likelihood = 0
                else:
                    idq_glitch_row.likelihood = numpy.infty
            else:
                idq_glitch_row.likelihood = Effmap(idq_glitch_row.rank) / idq_glitch_row.fap

                # ## add info to OVLDataTable() objects

        if flavor == 'ovl':
            ovl_row = idq_tables.OVLData()  # define row object
            ovl_row.ifo = ifo
            ovl_row.aux_channel = dat['vchan'][ind]
            ovl_row.veto_thr = float(dat['vthr'][ind])
            ovl_row.veto_win = float(dat['vwin'][ind])

        ### define coincs, update rows, etc. Do this all at once so the index counters don't get confused
        if i == 1:  # glitch sample
            idq_glitch_row.event_id = gchTable.get_next_id()
            gchTable.append(idq_glitch_row)
            if gwtrigs: # add gw trigger to sngl_burst table and add idq_glitch<-->sngl_burst coinc
                sngl_burst_row = lsctables.SnglBurst()
                sngl_burst_columns = lsctables.SnglBurstTable.validcolumns

                gw_trig = gwtrigs[glitch_index] # columns for kw triggers are defined in event.py

                sngl_burst_row.process_id = gchproc.process_id
                sngl_burst_row.ifo = ifo
                sngl_burst_row.search = "kleineWelle"
                sngl_burst_row.channel = gwchan
                sngl_burst_row.set_start(LIGOTimeGPS(gw_trig[0]))
                #sngl_burst_row.start_time = int(gw_trig[0])
                #sngl_burst_row.start_time_ns = int( round( gw_trig[0] % 1.0, ndigits=9 ) * 10**9 )
                sngl_burst_row.set_peak(LIGOTimeGPS(gw_trig[2]))
                #sngl_burst_row.peak_time = int(gw_trig[2])
                #sngl_burst_row.peak_time_ns = int( round( gw_trig[2] % 1.0, ndigits=9 ) * 10**9 )
                sngl_burst_row.set_stop(LIGOTimeGPS(gw_trig[1]))
                sngl_burst_row.duration = (sngl_burst_row.stop_time - sngl_burst_row.start_time) + (sngl_burst_row.stop_time_ns - sngl_burst_row.start_time_ns)*10**(-9)
                sngl_burst_row.central_freq = gw_trig[3]
                sngl_burst_row.snr = math.sqrt(gw_trig[-3] - gw_trig[-2])
                sngl_burst_row.confidence = gw_trig[7]
                sngl_burst_row.event_id = gchSnglBurstTable.get_next_id()

                # define list of columns which we filled in
                used_valid_columns = ['process_id', 'ifo', 'search', 'channel', 'start_time', 'start_time_ns', 'stop_time', 'stop_time_ns',
                                      'peak_time', 'peak_time_ns', 'duration', 'central_freq',
                                      'snr', 'confidence', 'event_id']

                # fill all other columns with None
                for (column, type) in sngl_burst_columns.iteritems():
                    if not column in used_valid_columns: setattr(sngl_burst_row, column, None)

                gchSnglBurstTable.append(sngl_burst_row)

                # creat idq_glitch<-->sngl_burst coinc event
                coinc_row = lsctables.Coinc()
                coinc_row.process_id = gchproc.process_id
                coinc_row.coinc_def_id = gchCoincDefTable.get_coinc_def_id(
                    search=idq_tables.IDQCoincDef['idq_glitch<-->sngl_burst'][0],\
                    search_coinc_type=idq_tables.IDQCoincDef['idq_glitch<-->sngl_burst'][1],\
                    description = 'idq_glitch<-->sngl_burst'
                    )
                coinc_row.time_slide_id = lsctables.TimeSlideID(0)
                coinc_row.instruments = ifo
                coinc_row.nevents = 1
                coinc_row.likelihood = idq_glitch_row.likelihood
                coinc_row.coinc_event_id = gchCoincTable.get_next_id()

                gchCoincTable.append(coinc_row)

                # add entries to coinc_map table
                coinc_map_row = lsctables.CoincMap()
                coinc_map_row.coinc_event_id = coinc_row.coinc_event_id
                coinc_map_row.table_name = idq_tables.IDQGlitchTable.tableName
                coinc_map_row.event_id = idq_glitch_row.event_id

                gchCoincMapTable.append(coinc_map_row)

                coinc_map_row = lsctables.CoincMap()
                coinc_map_row.coinc_event_id = coinc_row.coinc_event_id
                coinc_map_row.table_name = lsctables.SnglBurstTable.tableName
                coinc_map_row.event_id = sngl_burst_row.event_id

                gchCoincMapTable.append(coinc_map_row)

            if flavor == 'ovl':
                ovl_row.event_id = gchOVLTable.get_next_id()
                gchOVLTable.append(ovl_row)

                # creat idq_glitch<-->ovl_data coinc event
                coinc_row = lsctables.Coinc()
                coinc_row.process_id = gchproc.process_id
                coinc_row.coinc_def_id = gchCoincDefTable.get_coinc_def_id(
                    search=idq_tables.IDQCoincDef['idq_glitch<-->ovl_data'][0],\
                    search_coinc_type=idq_tables.IDQCoincDef['idq_glitch<-->ovl_data'][1],\
                    description = 'idq_glitch<-->ovl_data'
                    )
                coinc_row.time_slide_id = lsctables.TimeSlideID(0)
                coinc_row.instruments = ifo
                coinc_row.nevents = 1
                coinc_row.likelihood = idq_glitch_row.likelihood
                coinc_row.coinc_event_id = gchCoincTable.get_next_id()

                gchCoincTable.append(coinc_row)

                # add entries to coinc_map table
                coinc_map_row = lsctables.CoincMap()
                coinc_map_row.coinc_event_id = coinc_row.coinc_event_id
                coinc_map_row.table_name = idq_tables.IDQGlitchTable.tableName
                coinc_map_row.event_id = idq_glitch_row.event_id

                gchCoincMapTable.append(coinc_map_row)

                coinc_map_row = lsctables.CoincMap()
                coinc_map_row.coinc_event_id = coinc_row.coinc_event_id
                coinc_map_row.table_name = idq_tables.OVLDataTable.tableName
                coinc_map_row.event_id = ovl_row.event_id

                gchCoincMapTable.append(coinc_map_row)

            glitch_index += 1

        else:
            idq_glitch_row.event_id = clnTable.get_next_id()
            clnTable.append(idq_glitch_row)
            if flavor == 'ovl':
                ovl_row.event_id = clnOVLTable.get_next_id()
                clnOVLTable.append(ovl_row)

                # creat idq_glitch<-->ovl_data coinc event
                coinc_row = lsctables.Coinc()
                coinc_row.process_id = clnproc.process_id
                coinc_row.coinc_def_id = clnCoincDefTable.get_coinc_def_id(
                    search=idq_tables.IDQCoincDef['idq_glitch<-->ovl_data'][0],\
                    search_coinc_type=idq_tables.IDQCoincDef['idq_glitch<-->ovl_data'][1],\
                    description = 'idq_glitch<-->ovl_data'
                    )
                coinc_row.time_slide_id = lsctables.TimeSlideID(0)
                coinc_row.instruments = ifo
                coinc_row.nevents = 1
                coinc_row.likelihood = idq_glitch_row.likelihood
                coinc_row.coinc_event_id = clnCoincTable.get_next_id()

                clnCoincTable.append(coinc_row)

                # add entries to coinc_map table
                coinc_map_row = lsctables.CoincMap()
                coinc_map_row.coinc_event_id = coinc_row.coinc_event_id
                coinc_map_row.table_name = idq_tables.IDQGlitchTable.tableName
                coinc_map_row.event_id = idq_glitch_row.event_id

                clnCoincMapTable.append(coinc_map_row)

                coinc_map_row = lsctables.CoincMap()
                coinc_map_row.coinc_event_id = coinc_row.coinc_event_id
                coinc_map_row.table_name = idq_tables.OVLDataTable.tableName
                coinc_map_row.event_id = ovl_row.event_id

                clnCoincMapTable.append(coinc_map_row)

    return (gchxml_doc, clnxml_doc)

#===================================================================================================
# OVL WRAPPERS
#===================================================================================================
def ovl_evaluate( vetolist, GPS_headers=False, GPStimes=False, allvtrg=False, kw_trgfiles=False, gwchan=False, patfiles=False, skip_lines=1, filename='ovl_predict', output_dir='./'):
    """
    takes the last line of vetolist_cache file and uses that as a pointer to the most recent vetolist information.
    generates vetosegment from kw_trgfiles
    if patfiles:
        generates GPStimes from patfiles
    else:
        generates GPStimes from kw_trgfiles using gwchan (REQUIRED)

    expects vetolsit to be a file path
        kw_trgfiles, patfiles to be lists of file paths
        gwchan to be a string
    """
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
    return ovl.predict( vetolist, GPS_headers, GPStimes, gps_tcent, allvtrg=allvtrg, kw_trgfiles=kw_trgfiles, predict_filename=filename, output_dir=output_dir )

def ovl_train(gpsstart, gpsstop, generalD, classifierD, scisegs=False, vetosegs=False, output_dir='./' , padding=numpy.infty ):
    """ 
    builds an ovl.params object and launches ovl training jobs on the specified data.
    pulls many parameters from "cp" config object
    """
    # build params object
    analysis_range = [gpsstart, gpsstop]

    # load from cp object
    auxdir = generalD['snglchndir']
    gwdir = generalD['snglchndir']
    gwchans = generalD['gwchannel'].split()
    gwthr = float(generalD['gw_kwsignif_thr'])
    ifos = generalD['ifo'].split()

    metric = classifierD['metric']

    gwsets = classifierD['gwsets'].split()
    safety = classifierD['safety']
    windows = [float(l) for l in classifierD['windows'].split()]
    thresholds = [float(l) for l in classifierD['thresholds'].split()]
    Psigthr = float(classifierD['psigthr'])
    effbydtthr = float(classifierD['effbydtthr'])

    #channels = False
    #notused = []

    channels = generalD['selected-channels']
    notused = generalD['unsafe-channels']
    if channels:
        channels = False
    if notused:
        notused = [l.strip('\n') for l in open(notused, 'r').readlines() if l.strip('\n')]
    else:
        notused = []

    params = ovl.params( analysis_range, auxdir, gwdir, gwchans, gwthr, ifos, gwsets, scisegs=scisegs, vetosegs=vetosegs, channels=channels, notused=notused, windows=windows, thresholds=thresholds, Psigthr=Psigthr, effbydtthr=effbydtthr, safety=safety, metric=metric )

    # double check that windows are not bigger than padding set in idq_realtime
    len_windows = len(windows)
    windows = [w for w in windows if w <= padding]
    if len_windows > len(windows):
        print 'WARNING: ovl windows are not consistent with idq_realtime padding! %d windows were removed.' % (len_windows - len(windows))

    # load training parameters from cp object
    num_runs = int(classifierD['num_runs'])
    incremental = int(classifierD['incremental'])

    # launch training job
    vetolists = ovl.train(params, num_runs=num_runs, incremental=incremental, output_dir=output_dir, verbose=False, write_channels=True )

    ### run safety to get a vetolist file for "independent" application
    ovl.vetolist_safety(params, output_dir=output_dir+"/ovl/", source_dir=output_dir+"/ovl/", verbose=False)

    return vetolists

#===================================================================================================
# FOREST WRAPPERS
#===================================================================================================
def forest_evaluate( patfile, trainedforest, ranked_file, cp, gps_start_time, gps_end_time, dir):
    """
    Submits job that evaluates samples of auxmvc feature vectors using random forest (MVSC)
    """

    # initiate use forest job
    use_forest_job = auxmvc.use_forest_job(cp)

    # create node for this job
    use_forest_node = auxmvc.use_forest_node(use_forest_job, trainedforest, patfile, ranked_file)

    # get full command line for this job
    use_forest_command = auxmvc.construct_command(use_forest_node)

    # submit process
#    exit_status = submit_command(use_forest_command, 'forest_evaluate', dir)
    exit_status = subprocess.Popen(use_forest_command, cwd=dir).wait() ### block!

    if exit_status == 0:
        # run postscript
        forest_add_excluded_variables_job = auxmvc.forest_add_excluded_vars_job(cp)
        forest_add_excluded_variables_node = auxmvc.forest_add_excluded_vars_node(forest_add_excluded_variables_job, patfile, ranked_file)

        # get full command line for this job
        forest_add_excluded_variables_command = auxmvc.construct_command(forest_add_excluded_variables_node)

        # submit process
#        exit_status = submit_command(forest_add_excluded_variables_command, 'forest_add_excluded_variables', dir)
        exit_status = subprocess.Popen(forest_add_excluded_variables_command, cwd=dir).wait() ### block!

        return (exit_status, use_forest_node.get_output_files())

    else:
        return (exit_status, use_forest_node.get_output_files())


def execute_forest_train(training_samples_file, cache, cp, submit_dir ):
    """
    Builds small dag to train  random forest (MVSC) and condor submits it.
    """

    # set current directory to submit_dir
    os.chdir(submit_dir)

    # set path to condor log
    logpath = cp.get('idq_train', 'condorlogs')

    # set basename for condor dag
    basename = cp.get('general', 'ifo') + '_train_mvsc_' + cp.get('general', 'usertag') + '-' + training_samples_file.split('-')[-2] + '-' + training_samples_file.split('-')[-1].split('.')[0]

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
    #trained_forest_filename = os.path.split(training_samples_file)[0] + '/mvsc/' + os.path.split(training_samples_file)[1].replace('.pat', '.spr')
    trained_forest_filename = submit_dir +'/' + os.path.split(training_samples_file)[1].replace('.pat', '.spr')

    # create node for this job
    train_forest_node = auxmvc.train_forest_node(train_forest_job, training_samples_file, trained_forest_filename)

    # append the node to dag
    dag.add_node(train_forest_node)

    # initiate add file to cache job
    add_file_to_cache_job = auxmvc.add_file_to_cache_job(cp)

    # create node for this job
    add_file_to_cache_node = auxmvc.add_file_to_cache_node(add_file_to_cache_job, [train_forest_node.trainedforest], cache, p_node=[train_forest_node])

    # add the node to dag
    dag.add_node(add_file_to_cache_node)

    # write dag
    dag.write_sub_files()
    dag.write_dag()
    dag.write_script()

    # condor dag submit command....
    dag_submit_cmd = ['condor_submit_dag', dag_file]

    # submit dag
#    exit_status = submit_command(dag_submit_cmd, 'submit_forest_train_dag', submit_dir)
    exit_status = subprocess.Popen(dag_submit_cmd, cwd=submit_dir).wait() ### block!

    return (exit_status, dag_file, train_forest_node.get_output_files())

#===================================================================================================
# SVM WRAPPERS
#===================================================================================================

def svm_evaluate( cp, test_file, range_file, model_file, predict_file, dir ):

    use_svm_job = auxmvc.use_svm_job(cp)
    use_svm_node = auxmvc.use_svm_node( use_svm_job, cp, test_file, range_file, model_file, predict_file )
    use_svm_command = auxmvc.construct_command(use_svm_node)

#    exit_status = submit_command(use_svm_command, 'svm_evaluate', dir)
    exit_status = subprocess.Popen(use_svm_command, cwd=dir).wait() ### block!

    return exit_status


def execute_svm_train( cp, train_file, range_file, model_file, cache_file, submit_dir ):
    """
    Builds small condor dag to train SVM and submits it.
    """
    # set current directory to submit_dir
    os.chdir(submit_dir)

    # set path to condor log
    logpath = cp.get('idq_train', 'condorlogs')

    # set basename for condor dag
    basename = cp.get('general', 'ifo') + '_train_svm_'  + cp.get('general', 'usertag') + '-' + train_file.split('-')[-2] + '-' + train_file.split('-')[-1].split('.')[0]

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
    train_svm_node = auxmvc.train_svm_node(train_svm_job, cp, train_data, range_data, model_data)

    # append the node to dag
    dag.add_node(train_svm_node)

    # initiate add file to cache job
    add_file_to_cache_job = auxmvc.add_file_to_cache_job(cp)

    # create node for this job
    add_file_to_cache_node = auxmvc.add_file_to_cache_node(add_file_to_cache_job, [submit_dir + train_svm_node.model_file, submit_dir + train_svm_node.range_file], cache_file, p_node=[train_svm_node])

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
#    exit_status = submit_command(dag_submit_cmd, 'submit_svm_train_dag', submit_dir)
    exit_status = subprocess.Popen(dag_submit_cmd, cwd=submit_dir).wait() ### block!

    return (exit_status, dag_file, train_svm_node.get_output_files())

#===================================================================================================
# ANN WRAPPERS
#===================================================================================================
def ann_evaluate( patfile, trained_ann, ranked_file, cp, gps_start_time, gps_end_time, dir):
    """
    Submits job that evaluates samples of auxmvc feature vectors using random ann (ANN)
    """

    # initiate converting pat_file to ann file(FANN library type file) job
    convert_annfile_job = auxmvc.convert_annfile_job(cp)

    # create node for convert_annfile_job
    convert_annfile_node = auxmvc.convert_annfile_node(convert_annfile_job, patfile)

    # initiate use ann job
    use_ann_job = auxmvc.use_ann_job(cp)

    # create node for this job
    use_ann_node = auxmvc.use_ann_node(use_ann_job, trained_ann, patfile, ranked_file)

    # get full command line for convert_annfilejob
    convert_annfile_command = auxmvc.construct_command(convert_annfile_node)

    # submit process
    #exit_status = submit_command(convert_annfile_command, 'ann_convert', dir)
    exit_status = subprocess.Popen(convert_annfile_command, cwd=dir).wait() ### block

    if exit_status == 0:
        # get full command line for this job
        use_ann_command = auxmvc.construct_command(use_ann_node)

        # submit process
        #exit_status = submit_command(use_ann_command, 'ann_evaluate', dir)
        exit_status = subprocess.Popen(use_ann_command, cwd=dir).wait() ### block!

        return (exit_status, use_ann_node.get_output_files())
    else:
        print "File conversion for ann is failed."
        return (exit_status, use_ann_node.get_output_files())

def execute_ann_train( training_samples_file, cache, cp, submit_dir ):
    """
    Builds small dag to train ANN and condor submits it.
    """

    # set current directory to submit_dir
    os.chdir(submit_dir)

    # set path to condor log
    logpath = cp.get('idq_train', 'condorlogs')

    # set basename for condor dag
    basename = cp.get('general', 'ifo') + '_train_ann_' + cp.get('general', 'usertag') + '-' + training_samples_file.split('-')[-2] + '-' + training_samples_file.split('-')[-1].split('.')[0]

    # creat directory for jobs .err and .out files
    if not os.path.exists('logs'):
        os.makedirs('logs')

    # initiate dag
    dag = auxmvc.auxmvc_DAG(basename, logpath)

    # get dag file
    dag_file = dag.get_dag_file()

    # initiate convert_annfile_job
    convert_annfile_job = auxmvc.convert_annfile_job(cp)

    # create node for convert_annfile_job
    convert_annfile_node = auxmvc.convert_annfile_node(convert_annfile_job, training_samples_file)

    # append the node to dag
    dag.add_node(convert_annfile_node)

    # construct name for trained ann file
    #training_samples_file_redirect = os.path.split(training_samples_file)[0] + '/ann/' + os.path.split(training_samples_file)[1].replace('.pat', '.ann')
    training_samples_file_redirect = submit_dir +'/' + os.path.split(training_samples_file)[1].replace('.pat', '.ann')

    trained_ann_filename = training_samples_file_redirect.replace('.ann','.net')

    # initiate train ann job
    train_ann_job = auxmvc.train_ann_job(cp)

    # create node for this job
    train_ann_node = auxmvc.train_ann_node(train_ann_job, training_samples_file_redirect, trained_ann_filename, p_node=[convert_annfile_node])

    # append the node to dag
    dag.add_node(train_ann_node)

    # initiate add file to cache job
    add_file_to_cache_job = auxmvc.add_file_to_cache_job(cp)

    # create node for this job
    add_file_to_cache_node = auxmvc.add_file_to_cache_node(add_file_to_cache_job, [train_ann_node.trained_ann], cache, p_node=[train_ann_node])

    # add the node to dag
    dag.add_node(add_file_to_cache_node)

    # write dag

    dag.write_sub_files()
    dag.write_dag()
    dag.write_script()

    # condor dag submit command....

    dag_submit_cmd = ['condor_submit_dag', dag_file]

    # submit dag
    #exit_status = submit_command(dag_submit_cmd, 'submit_ann_train_dag', submit_dir)
    exit_status = subprocess.Popen(dag_submit_cmd, cwd=submit_dir).wait() ### block!

    return (exit_status, dag_file, train_ann_node.get_output_files())

##@}
