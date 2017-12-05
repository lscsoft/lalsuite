# Copyright (C) 2013 Lindy Blackburn
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

# Last modified: Jul 11, 2013 LLB
#     - added loadkwconfig for reading kw configuration file as dictionary
#     - added random times routine
#     - added random times from rate routine
#     - added loadkwsingles routine
#
# original version from
# Python routines for basic trigger and segment processing
# Lindy Blackburn (lindy@ligo.mit.edu)
# Last modified: Jan 9, 2013 L. Blackburn
#     - trigger lists and segment lists are stored as two-dimensional lists,
#       even for a single trigger/segment
#     - it is important to use the correct column for central time
#     - andsegments/orsegments work on a list of segment lists, a 3-D list
#     - segment lists assumed to be sorted, continuous and non-overlapping
#       segments are treated as half-open intervals [start, stop)
#     - triggers are not automatically sorted when they are read from disk
#       functions will create a sorted copy of the triggers when called
#     - added routines for loading segments (including XML)
#     - added routine to load triggers with cuts
#     - added ncoinc and cluster from backup version
#     - added loadkwm for reading multichannel KW triggers into a dict

## \defgroup laldetchar_py_idq_event Event Module
## \ingroup laldetchar_py_idq
## Synopsis
# ~~~
# from laldetchar.idq import event
# ~~~
# \author Lindy Blackburn (<lindy.blackburn@ligo.org>)

""" Module with classes and functions for construction and manipulation of composite glitch events.
"""

import math
import sys
import numpy
import random

from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
from glue.ligolw import table

from laldetchar import git_version

__author__ = 'Lindy Blackburn <lindy.blackburn@ligo.org>, Reed Essick <reed.essick@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

## \addtogroup laldetchar_py_idq_event
# @{

# Classes and functions for construction and manipulation of
# composite glitch events.

col_kw = {
    'tstart': 0,
    'tstop': 1,
    'tcent': 2,
    'fcent': 3,
    'uenergy': 4,
    'nenergy': 5,
    'npix': 6,
    'signif': 7,
    }
col_veto = {
    'tcent': 0, 
    'signif': 1,
}

col_snglBurst = {
    "tstart":0, 
    "duration":1, 
    "tcent":2, 
    "fpeak":3, 
    "fcent":4, 
    "bandwidth":5, 
    "amplitude":6, 
    "snr":7, 
    "conf":8, 
    "chi2":9, 
    "chi2_dof":10,
}

ifo = dict()
ifo['kwh1'] = {'H1': 2}
ifo['kwh2'] = {'H2': 2}
ifo['kwl1'] = {'L1': 2}
ifo['kwv1'] = {'V1': 2}
ifo['cwb'] = {'l1': 23, 'h1': 24, 'h2': 25}
ifo['cwbh1h2'] = {'h1': 20, 'h2': 21}
ifo['cwbhlv'] = {
    'h1': 27,
    'h2': 28,
    'l1': 26,
    'v1': 29,
    }
ifo['qh1h2'] = {'h1': 0, 'h2': 9}
ifo['qh1'] = {'h1': 0}
ifo['ql1'] = {'l1': 0}
ifo['qhfh1h2'] = {'h1': 0, 'h2': 1}

col = col_kw.copy()


# take subset of triggers which fall within segments
#     - segments are interpreted as [start, stop)
#     - tlag is added to the trigger times
#     - tcent defines the column for central time and is optional

def include(
    triggers,
    segments,
    tlag=0,
    tcent=None,
    ):
    if tcent == None:
        tcent = col['tcent']
    sortedtriggers = triggers[:]
    sortedtriggers.sort(key=lambda x: x[tcent])
    subset = []
    i = 0
    for seg in segments:
        while i < len(sortedtriggers) and sortedtriggers[i][tcent] \
            + tlag < seg[0]:
            i += 1
        while i < len(sortedtriggers) and sortedtriggers[i][tcent] \
            + tlag < seg[1]:
            subset.append(sortedtriggers[i])
            i += 1
    return subset


# take subset of triggers which fall outside of segments
#     - segments are interpreted as [start, stop)
#     - tlag is added to the trigger times
#     - tcent defines the column for central time and is optional

def exclude(
    triggers,
    segments,
    tlag=0,
    tcent=None,
    ):
    if tcent == None:
        tcent = col['tcent']
    sortedtriggers = triggers[:]
    sortedtriggers.sort(key=lambda x: x[tcent])
    subset = []
    i = 0
    for seg in segments:
        while i < len(sortedtriggers) and sortedtriggers[i][tcent] \
            + tlag < seg[0]:
            subset.append(sortedtriggers[i])
            i += 1
        while i < len(sortedtriggers) and sortedtriggers[i][tcent] \
            + tlag < seg[1]:
            i += 1
    while i < len(sortedtriggers):
        subset.append(sortedtriggers[i])
        i += 1
    return subset


# return [included, excluded] triggers given segments

def includeexclude(
    triggers,
    segments,
    tlag=0,
    tcent=None,
    ):
    if tcent == None:
        tcent = col['tcent']
    sortedtriggers = triggers[:]
    sortedtriggers.sort(key=lambda x: x[tcent])
    included = []
    excluded = []
    i = 0
    for seg in segments:
        while i < len(sortedtriggers) and sortedtriggers[i][tcent] \
            + tlag < seg[0]:
            excluded.append(sortedtriggers[i])
            i += 1
        while i < len(sortedtriggers) and sortedtriggers[i][tcent] \
            + tlag < seg[1]:
            included.append(sortedtriggers[i])
            i += 1
    while i < len(sortedtriggers):
        excluded.append(sortedtriggers[i])
        i += 1
    return [included, excluded]


# load triggers that pass threshold from ascii file into list
#     - file: filename string. '.gz' files are supported, and will be
#       tried automatically if file does not exist
#     - only triggers with significance >= threshold are loaded
#     - signif defines the column for significance and is optional

def loadtrg(file, threshold=0, signif=None):
    if signif == None:
        signif = col['signif']
    triggers = []
    f = gzopen(file, 'r')
    lines = f.readlines()
    f.close()
    lines = list(set(lines))  # remove duplicate lines
    for line in lines:
        line = line.strip()
        if line[0] == '#' or line == '':
            continue
        fields = line.split()
        if float(fields[signif]) >= threshold:
            triggers.append([float(num) for num in fields])
    return triggers


# load triggers that pass cuts
#     - file: filename string. '.gz' files are supported, and will be
#       tried automatically if file does not exist
#     - cuts is a list of lists of the form: [col, r1, r2]
#       if r1 < r2: r1 <= trg[col] AND trg[col] < r2 for the trigger to pass (inclusive)
#       may use interval with r1 > r2 to make an exclusion interval:
#       if r1 >= r2: trg[col] < r2 OR trg[col] <= r1 to pass
#     - each trigger must pass ALL cuts
#     - duplicate lines are removed
#     - does NOT preserve order of triggers (due to duplicate line removal)

def loadwithcuts(file, cuts=[]):
    if cuts != [] and dim(cuts) == 1:
        cuts = [cuts]
    triggers = []
    f = gzopen(file, 'r')
    lines = f.readlines()
    f.close()
    lines = list(set(lines))  # remove duplicate lines
    for line in lines:
        line = line.strip()
        if line[0] == '#' or line == '':
            continue
        fields = line.split()
        reject = False
        for cut in cuts:
            (val, r1, r2) = (float(fields[cut[0]]), cut[1], cut[2])
            if r1 < r2:
                reject = val < r1 or val >= r2  # inclusion interval [r1, r2)
            else:
                reject = val < r1 and val >= r2  # exclusion interval [r2, r1)

            # print (val, r1, r2, val < r1, val >= r2, reject)

            if reject:
                break  # break on first cut that fails
        if not reject:
            triggers.append([float(num) for num in fields])
    return triggers


# load triggers light version for veto work, only load tcent and signif
#     - file: filename string. '.gz' files are supported, and will be
#       tried automatically if file does not exist
#     - only triggers with significance >= threshold are loaded
#     - tcent defines the column for central time and is optional
#     - signif defines the column for significance and is optional

def loadvetotrg(
    file,
    threshold=0,
    tcent=None,
    signif=None,
    ):
    if tcent == None:
        tcent = col['tcent']
    if signif == None:
        signif = col['signif']
    triggers = []
    f = gzopen(file, 'r')
    lines = f.readlines()
    f.close()
    lines = list(set(lines))  # remove duplicate lines
    for line in lines:
        line = line.strip()
        if line[0] == '#' or line == '':
            continue
        fields = line.split()
        if float(fields[signif]) >= threshold:
            trg = [float(fields[tcent]), float(fields[signif])]
            triggers.append(trg)
    return triggers


# load in ASCII table as float values, '.gz' supported

def load(file):
    f = gzopen(file, 'r')

    # table = [[float(field) for field in line.strip().split()] for line in f if line[0] != '#' and line.strip() != '' and line[-1] == '\n'] // necessary if there may be incomplete lines

    table = [[float(field) for field in line.strip().split()]
             for line in f if line[0] != '#' and line.strip() != '']
    f.close()
    return table


def loadcwb(file):
    f = gzopen(file, 'r')
    table = [[1] + [float(field) for field in line.strip().split()[1:]]
             for line in f if line[0] == '+']
    f.close()
    return table


def loadsegments(file):
    f = gzopen(file, 'r')
    if file[-4:] == '.xml' or file[-7:] == '.xml.gz':
        from xml.dom import minidom
        xmldoc = minidom.parse(f)
        table = [map(float, line.split(',')[3:5]) for line in
                 xmldoc.lastChild.childNodes[-2].childNodes[-2].firstChild.data.strip().split('\n'
                 )]
        return table
    else:
        table = [map(float, line.strip().split()) for line in f
                 if line[0] != '#' and line.strip() != '']
        if table == []:
            return table
        else:
            idx = table[0].index(max(table[0]))
            segments = [[line[idx - 1], line[idx]] for line in table]
            return segments


# load in ASCII table as string values, '.gz' supported

def loadstringtable(file):
    f = gzopen(file, 'r')
    table = [line.strip().split() for line in f if line[0] != '#'
             and line.strip() != '']
    f.close()
    return table


# load in single string list

def loadlist(file):
    f = open(file, 'r')
    list = [line.strip() for line in f if line[0] != '#'
            and line.strip() != '']
    f.close()
    return list


def read_channels_from_file(file):
    """
    Reads channel names from file. Return list of channel names.
    """

    channels = open(file).readlines()
    channels = [ch.strip() for ch in channels if ch.strip()]
    return channels


class trigdict(dict):

    """
    Class to hold multi-channel triggers. It is a dictionary, with keys defined by the names of the channels. Each value is a list of triggers.
    triggers are lists themselves, with entries correposnding to kw_col, defined at the top.
    """

    def __init__(self, channels=None):

        # self = {}

        if channels:
            for channel in channels:
                self[channel] = []

        # self.channels = self.keys()
        # return self

    def channels(self):
        return self.keys()

    def add(self, otherdict):
        for key, value in otherdict.items():
            if self.has_key(key):
                self[key] += value
            else:
                self[key] = value

    def get_triggers_from_channel(self, channel):
        """
        Return triggers from a given channel
        """

        if channel in self.keys():
            return self[channel]
        else:
            raise NameError('The specified channel is not recognized. It is not in the dictionary keys'
                            )

    def apply_signif_threshold(
        self,
        channels=[],
        threshold=0,
        signif=None,
        ):
        """
        Aplly threshold on triggers significance, by getting rid of those triggers which significance is less than the threshold.
        channels is the lits of channels in which threshold is applied.
        """

        if signif == None:
            signif = col['signif']

        # check if all channels in the channels are recognized

        for channel in channels:
            if not channel in self.channels():
                raise NameError('One or more channels in the list are not recognized.'
                                )

        for channel in channels:
            triggers_above_thr = []
            for trigger in self[channel]:
                if trigger[signif] >= threshold:
                    triggers_above_thr.append(trigger)
            self[channel] = triggers_above_thr

    def remove_channels(self, channels):
        """
        Removes unsafe channels from dictionary of triggers.
        channels is the list of channel names.
        """

        for key in self.keys():
            if key in channels:
                del self[key]

    def keep_channels(self, channels):
        """
        Keeps triggers only from the designated channels.
        channels is the list of channel names.
        """

        for key in self.keys():
            if not key in channels:
                del self[key]

    def include(
        self,
        segments,
        channels=None,
        tlag=0,
        tcent=None,
        ):
        """
        Keeps triggers that fall within the segments.
        If channels are given, action is applied only to triggers from these channels.
        """

        if not channels:
            channels = self.channels()
        for channel in channels:
            self[channel] = include(self.get_triggers_from_channel(channel), segments, tlag, tcent)

    def resort(self, tcent=col['tcent']):
        """
        sort triggers by central time. optionally define key
        """

        for triglist in self.values():
            triglist.sort(key=lambda x: x[tcent])

def loadSingleBurst( files, trigs_dict=None):
    """
    loads snglburst tables (produced by Omicron) into trgdict object
    files - is the list of file names
    """
    if type(files) is str:
        files = [files]
    if trigs_dict is None:
        trigs_dict = trigdict()
    for file in files:
        for row in table.get_table( ligolw_utils.load_filename(file, contenthandler=lsctables.use_in(ligolw.LIGOLWContentHandler)), lsctables.SnglBurstTable.tableName ):
            channel = "%s-%s_%s"%(row.ifo, row.channel.replace("-","_"), row.search)
            tcent = row.peak_time + 1e-9*row.peak_time_ns
            tstart = row.start_time + 1e-9*row.start_time_ns
            dur = row.duration
            fpeak = row.peak_frequency
            fcent = row.central_freq
            bndwth = row.bandwidth
            amp = row.amplitude
            snr = row.snr
            conf = row.confidence
            chi2 = row.chisq
            chi2_dof = row.chisq_dof

            trigger = [tstart, dur, tcent, fpeak, fcent, bndwth, amp, snr, conf, chi2, chi2_dof]

            if channel in trigs_dict.channels():
                trigs_dict[channel].append( trigger ) ### SingleBurst trigger structure
            else:
                trigs_dict[channel] = [ trigger ]

    return trigs_dict

# load in KW multichannel trigger file and merge into trigs_dict
def loadkwm(files, trigs_dict=None):
    """
    Loads multi-channel KW trigger files. Returns a dictionary with channel names used as keys.
    files - is the list of file names.
    """

    if type(files) is str:
        files = [files]
    if trigs_dict is None:
        trigs_dict = trigdict()
    for file in files:
        table = loadstringtable(file)
        for line in table:
            (trigger, channel) = ([float(val) for val in line[:-1]],
                                  line[-1])
            if channel in trigs_dict.channels():
                trigs_dict[channel].append(trigger)
            else:
                trigs_dict[channel] = [trigger]
    return trigs_dict


# load in KW singles triggers from many channels into a trigs_dict
# if channels is specified, then only those channels are read (and [] list created for no file)
# if trigs_dict is specified, then triggers are merged into existing dict, otherwise new one is formed

def loadkwsingles(directory, channels=None, trigs_dict=None):
    import os
    import glob
    if channels is None:
        trgfiles = glob.glob(directory + '/*.trg')
    else:
        trgfiles = ['%s/%s.trg' % (directory, channel) for channel in
                    channels]
    if trigs_dict is None:
        trigs_dict = trigdict()
    for file in trgfiles:
        channel = file[1 + len(directory):-4]
        if channel in trigs_dict.channels():
            trigs_dict[channel].append(load(file))
        else:
            trigs_dict[channel] = load(file)
    return trigs_dict


def random_trigs(times):
    """
    Generates fake triggers for corresponding times.
    """

    return [[
        0.0,
        0.0,
        time,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        ] for time in times]

    {
        'tstart': 0,
        'tstop': 1,
        'tcent': 2,
        'fcent': 3,
        'uenergy': 4,
        'nenergy': 5,
        'npix': 6,
        'signif': 7,
        }


def build_auxmvc_vectors(
    trigdict,
    main_channel,
    coincidence_time_window,
    dqsegs=None,
    build_clean_samples=False,
    clean_times=None,
    clean_window=None,
    ):
    """
    Construct auxmvc feature vectors using triggers from multi-channel trig dictionary. Triggers from the channel designated as main are used as glitch indicators.
    """

    darmtrg = trigdict.get_triggers_from_channel(main_channel)
    darmtrg.sort(key=lambda x: x[col['tcent']])  # sort by central time
    if dqsegs == None:

        # if dqsegs are not given, central times of the triggers in main channel are used to define them

        if not clean_window:
            clean_window = coincidence_time_window
        dqsegs = vetosegs(darmtrg, clean_window)

    if build_clean_samples:

        # replace real triggers from the main channel by the fake ones.

        darmtrg = random_trigs(clean_times)

    chanlist = [channel for channel in trigdict.channels() if channel
                != main_channel]

    # must sort channel list to ensure the header is the same from file to file.

    chanlist.sort()

    header = 'GPS signif SNR unclean'
    fmt = ['g8', 'g8', 'g8', 'i']
    for chan in chanlist:
        fmt += ['g8', 'g8', 'g8', 'g8', 'g8']
        header += ' %s_signif %s_dt %s_dur %s_freq %s_npts' % (chan,
                chan, chan, chan, chan)

    variables = header.split(' ')
    auxmvctrigs = numpy.empty((len(darmtrg), ),
                              dtype={'names': variables,
                              'formats': fmt})

    trainingset = []
    idxa = dict()
    idxb = dict()
    for chan in chanlist:
        trigdict[chan].sort(key=lambda x: x[col['tcent']])  # sort by central time
        idxa[chan] = 0
        idxb[chan] = 0
    for (i, line) in enumerate(darmtrg):
        indq = []
        if len(include([line], dqsegs)) > 0:
            indq.append(1)
        else:
            indq.append(0)
        auxvars = []
        for chan in chanlist:
            subseg = [[line[col['tcent']] - coincidence_time_window,
                      line[col['tcent']] + coincidence_time_window]]
            if trigdict[chan]:  # if aux triggers exist for this channel
                while idxa[chan] + 1 < len(trigdict[chan]) \
                    and trigdict[chan][idxa[chan]][col['tcent']] \
                    < subseg[0][0]:
                    idxa[chan] += 1
                while idxb[chan] + 1 < len(trigdict[chan]) \
                    and trigdict[chan][idxb[chan] + 1][col['tcent']] \
                    <= subseg[0][1]:
                    idxb[chan] += 1
                subtrg = include((trigdict[chan])[idxa[chan]:idxb[chan]
                                 + 1], subseg)  # may still need to filter if idx=0
            else:
                subtrg = []
            if len(subtrg) > 0:
                besttrig = None
                bestsnr = 0
                for trg in subtrg:
                    dt = abs(trg[col['tcent']] - line[col['tcent']])
                    if trg[col['signif']] - dt / 1e3 > bestsnr:  # use dt as a tie-breaker
                        besttrig = trg
                        bestsnr = trg[col['signif']] - dt / 1e3
            else:
                besttrig = None

            if besttrig == None:
                auxvars += [0, random.uniform(-coincidence_time_window,
                            coincidence_time_window), 0, 0, 0]
            else:
                auxvars.append(besttrig[col['signif']])  # signif
                auxvars.append(besttrig[col['tcent']] - line[col['tcent'
                               ]])  # dt
                auxvars.append(besttrig[col['tstop']]
                               - besttrig[col['tstart']])  # duration
                auxvars.append(besttrig[col['fcent']])  # freq
                auxvars.append(besttrig[col['npix']])  # npts
        auxmvcvector = [line[col['tcent']], line[col['signif']],
                        math.sqrt(line[col['nenergy']] - line[col['npix'
                        ]])] + indq + auxvars
        for (var, val) in zip(variables, auxmvcvector):
            auxmvctrigs[var][i] = val
    return auxmvctrigs


def build_auxmvc_vectors_at_gps(
    gps,
    trigdict,
    main_channel,
    coincidence_time_window,
    ):
    """
        Construct auxmvc feature vectors using triggers from multi-channel trig dictionary at the specified gps time(s)
        """

    if isinstance(gps, (int, float)):
        gps = [gps]

    gps.sort()

    chanlist = [channel for channel in trigdict.channels() if channel
                != main_channel]

        # must sort channel list to ensure the header is the same from file to file.

    chanlist.sort()

    header = 'GPS signif SNR unclean'
    fmt = ['g8', 'g8', 'g8', 'i']
    for chan in chanlist:
        fmt += ['g8', 'g8', 'g8', 'g8', 'g8']
        header += ' %s_signif %s_dt %s_dur %s_freq %s_npts' % (chan,
                chan, chan, chan, chan)

    variables = header.split(' ')
    auxmvctrigs = numpy.empty((len(gps), ), dtype={'names': variables,
                              'formats': fmt})

    trainingset = []
    idxa = dict()
    idxb = dict()
    for chan in chanlist:
        trigdict[chan].sort(key=lambda x: x[col['tcent']])  # sort by central time
        idxa[chan] = 0
        idxb[chan] = 0

    for (i, line) in enumerate(gps):
        indq = [1]
        auxvars = []
        for chan in chanlist:
            subseg = [[line - coincidence_time_window, line
                      + coincidence_time_window]]
            if trigdict[chan]:  # if aux triggers exist for this channel
                while idxa[chan] + 1 < len(trigdict[chan]) \
                    and trigdict[chan][idxa[chan]][col['tcent']] \
                    < subseg[0][0]:
                    idxa[chan] += 1
                while idxb[chan] + 1 < len(trigdict[chan]) \
                    and trigdict[chan][idxb[chan] + 1][col['tcent']] \
                    <= subseg[0][1]:
                    idxb[chan] += 1
                subtrg = include((trigdict[chan])[idxa[chan]:idxb[chan]
                                 + 1], subseg)  # may still need to filter if idx=0
            else:
                subtrg = []
            if len(subtrg) > 0:
                besttrig = None
                bestsnr = 0
                for trg in subtrg:
                    dt = abs(trg[col['tcent']] - line)
                    if trg[col['signif']] - dt / 1e3 > bestsnr:  # use dt as a tie-breaker
                        besttrig = trg
                        bestsnr = trg[col['signif']] - dt / 1e3
            else:
                besttrig = None

            if besttrig == None:
                auxvars += [0, random.uniform(-coincidence_time_window,
                            coincidence_time_window), 0, 0, 0]
            else:
                auxvars.append(besttrig[col['signif']])  # signif
                auxvars.append(besttrig[col['tcent']] - line)  # dt
                auxvars.append(besttrig[col['tstop']]
                               - besttrig[col['tstart']])  # duration
                auxvars.append(besttrig[col['fcent']])  # freq
                auxvars.append(besttrig[col['npix']])  # npts
        auxmvcvector = [line, 0.0, 0.0] + indq + auxvars

        for (var, val) in zip(variables, auxmvcvector):
            auxmvctrigs[var][i] = val

    return auxmvctrigs


# write ASCII table with optional formatstring (printf style)
#     - filename ending in '.gz' will create gzip file

def write(table, file, formatstring=None):
    f = gzopen(file, 'w')
    if formatstring == None:
        table = [' '.join([repr(field) for field in line]) for line in
                 table]
    else:
        table = [formatstring % tuple(line) for line in table]
    for line in table:
        print >> f, line
    f.close()


# return list of non-overlapping segments by applying window about tcent
#     - triggers: list of triggers with tcent central time column
#     - window: scalar value gives a symmetric +/- window
#       list value gives segments [tcent+window[0], tcent+window[1])
#     - triggers with significance >= threshold contribute to veto segment
#     - tlag is added to the trigger times before making segments

def vetosegs(
    triggers,
    window,
    threshold=None,
    tlag=0,
    tcent=None,
    signif=None,
    ):
    if tcent == None:
        tcent = col['tcent']
    if signif == None:
        signif = col['signif']
    if isinstance(window, (int, long, float)):
        window = [-window, window]
    if threshold == None:
        sortedtriggers = triggers[:]
    else:
        sortedtriggers = [trg for trg in triggers if trg[signif]
                          >= threshold]
    sortedtriggers.sort(key=lambda x: x[tcent])
    segments = []
    for trigger in sortedtriggers:

        # segments is already populated

        if len(segments):

            # continuation of last segment

            if trigger[tcent] + tlag + window[0] \
                <= segments[len(segments) - 1][1]:
                segments[len(segments) - 1][1] = trigger[tcent] + tlag \
                    + window[1]
            else:

            # beginning of new segment

                segments.append([trigger[tcent] + tlag + window[0],
                                trigger[tcent] + tlag + window[1]])
        else:

        # segments is not yet populated

            segments.append([trigger[tcent] + tlag + window[0],
                            trigger[tcent] + tlag + window[1]])
    return segments


# take the intersection of segments
#     - segmentlists: list of segment lists that need to be intersected
#       e.g. gpstime = segmentlists[ifo][seg_number][0(start) or 1(end)]

def andsegments(segmentlists, wrongsyntax=None):
    if wrongsyntax != None:  # did not wrap two segment lists into a list
        return andsegments([segmentlists, wrongsyntax])
    if len(segmentlists) > 0 and isinstance(segmentlists[0], (int,
            long, float)):  # only one segment list
        return segmentlists
    if len(segmentlists) == 1:  # only one segment list in list
        return segmentlists[0]
    elif len(segmentlists) >= 2:
        lists = segmentlists[:]  # do not modify original list
        lists.sort(key=lambda x: len(x))  # loop over smallest lists first
        return reduce(lambda x, y: andtwosegmentlists(x, y), lists)  # (((a&b)&c)&d)& ...


# take the intersection of two segment lists: list1 and list2
#     - this is an internal function which is called by andsegments
#     - only use this instead if you know that the particular ordering of
#       list1/list2 is faster as andsegments(segmentlists) will sort the
#       segmentlists by their length

def andtwosegmentlists(list1, list2):
    newsegments = []
    index = 0
    for seg in list1:
        while index > 0 and list2[index][0] > seg[0]:
            index -= 1
        while index < len(list2) and list2[index][1] <= seg[0]:
            index += 1
        while index < len(list2) and list2[index][0] < seg[1]:
            newsegments.append([max(seg[0], list2[index][0]),
                               min(seg[1], list2[index][1])])
            index += 1
        if index > 0:
            index -= 1
    return newsegments


# take the n-fold intersection of segment lists
#  - segmentlists: list of segment lists
#  - n: overlap of n or more segments is included in output segments
#  - n should be <= the number of segment lists
#  - for n == number of segment lists, andsegments is faster

def nandsegments(segmentlists, n):
    starttimes = []
    stoptimes = []
    newsegments = []
    for segmentlist in segmentlists:
        for seg in segmentlist:
            starttimes.append([seg[0], 1])
            stoptimes.append([seg[1], -1])
    times = starttimes + stoptimes
    times.sort()
    coinc = 0
    seg = [0, 0]
    for time in times:
        coinc += time[1]
        if coinc == n and time[1] == 1:
            seg[0] = time[0]
        if coinc == n - 1 and time[1] == -1:
            seg[1] = time[0]

            # zero length segments can occur
            # if we avoid this by putting stop times first, segments can get broken up

            if seg[1] - seg[0] > 0:
                newsegments.append(seg[:])
    return newsegments


# take the union of segments
#     - segmentlists: list of segment lists that need to be merged.
#       e.g. gpstime = segmentlists[ifo][seg_number][0(start) or 1(end)]

def orsegments(segmentlists, wrongsyntax=None):
    if wrongsyntax != None:  # did not wrap two segment lists into a list
        return orsegments([segmentlists, wrongsyntax])
    if len(segmentlists) > 0 and isinstance(segmentlists[0], (int,
            long, float)):  # only one segment list
        return segmentlists
    if len(segmentlists) == 1:  # only one segment list in list
        return segmentlists[0]
    elif len(segmentlists) >= 2:

        # fixsegments turns out to be faster for this operation than a custom O(n) routine
        # because of Python's C sort which does profiling and O(n) mergesort

        return fixsegments(reduce(lambda x, y: x + y, segmentlists))  # fixsegments(a + b + c + ...)


# remove the "removesemgnets" from segments: segments - removesegments
#     - do not attempt to neglect to put in removesemgnets, use []
#       otherwise the function will calculate segments[0] - segments[1]

def removesegments(segments, removesegments=None):
    if removesegments == None:
        removesegments = segments[1]
        segments = segments[0]
    elif removesegments == []:
        return segments
    elif segments == []:
        return []
    notsegments = zip([a[1] for a in removesegments[:-1]], [a[0]
                      for a in removesegments[1:]])
    if segments[0][0] < removesegments[0][0]:
        notsegments.insert(0, [segments[0][0], removesegments[0][0]])
    if removesegments[-1][1] < segments[-1][1]:
        notsegments.append([removesegments[-1][1], segments[-1][1]])
    return andsegments([segments, notsegments])


# check segments for order and overlap, returns True for good segments

def checksegments(segments):
    return reduce(lambda x, y: [x[0] and y[0] > x[1] and y[1] - y[0]
                  > 0, y[1]], segments, [True, 0])[0]


# sort and merge an out-of-order or overlapping segment list

def fixsegments(segments):
    newsegments = []
    if segments == []:
        return newsegments
    segments.sort()
    nextseg = (segments[0])[:]
    for seg in segments:  # we want segments[1:] but this avoids a copy for large lists
        if seg[0] <= nextseg[1]:
            nextseg[1] = max(seg[1], nextseg[1])
        else:
            newsegments.append(nextseg)
            nextseg = seg[:]
    newsegments.append(nextseg)
    return newsegments


# cumulative livetime of segments

def livetime(segments):
    if len(segments) == 0:
        return 0
    elif isinstance(segments[0], (int, long, float)):
        segments = [segments]
    return sum(a[1] - a[0] for a in segments)


# cumulative sum of list

def cumsum(list):
    if list == []:
        return []
    else:
        return reduce(lambda x, y: x + [x[-1] + y], [list[0:1]]
                      + list[1:-1])


# dimensionality of a multidimensional list

def dim(multilist):
    list1 = multilist
    dim1 = 0
    while isinstance(list1, (list, tuple)):
        dim1 += 1
        if list1 == []:
            break
        list1 = list1[0]
    list2 = multilist
    dim2 = 0
    while isinstance(list2, (list, tuple)):
        dim2 += 1
        if list2 == []:
            break
        list2 = list2[-1]
    if dim1 == dim2:
        return dim1
    else:
        raise TypeError, \
            'first and last dimensions of the multidimensional list are not the same'


# open file or gzip file transparently, replaces open(file, mode)
#     - 'file.gz' is automatically tried for read modes if 'file' does
#        not exist. 'file.gz' in write mode will write a gzip file

def gzopen(file, mode='r'):
    import os

    # explicitly defined gzip file

    if file.endswith('.gz'):
        import gzip
        return gzip.GzipFile(file, mode)
    elif mode.startswith('r') and not os.path.exists(file):

    # if the file does not exist for read mode, try the gzip file
    # this will give a confusing error message for missing files

        import gzip
        return gzip.GzipFile(file + '.gz', mode)

    # this line should occur for 'r' modes where the file exists, or all 'w', 'a' modes

    return open(file, mode)


# ncoinc - n-fold coincidence, triggers must be sorted by central time
# channels: a list of trigger lists, each element will represent a different channel
# window:    time window for coidence (tcent is used +/- time window in seconds)
# triggers are assumed to be sorted by central time!!
# tlag: list of time lags for each trigger list in 'channels'
# Jan 11, 2008 L. Blackburn

def ncoinc(channels, window, tlag=[]):
    nchan = len(channels)

    # set time lags to zero if undefined

    if len(tlag) == 0:
        tlag = [0] * nchan

    # these are indexes for the triggers, and represent of the current range of triggers
    # being considered for coincidence with the first channel. we only consider a range
    # of trigger surrounding the first trigger at any one time for efficiency.

    head = [0] * nchan
    tail = [1] * nchan

    # trig is the index of the current set of triggers being considered for coincidence
    # it is a temporary variable and ultimately stores the indices of a coincident event

    trig = [0] * nchan

    # coinc is the matrix of coincident events. each element is a list of triggers corresponding to a coincident event

    coinc = []

    # loop over all triggers in the first channel. the first channel is 'special' and we consider them one at a time

    while head[0] < len(channels[0]):

        # increment pointers to cover range of first trigger. "tail" will always be one trigger outside window boundary.
        # tail index may fall one outside range of list (tail[i] = len(channels[i]), while head index will always be within list

        for i in range(1, nchan):
            while head[i] < len(channels[i]) \
                and channels[i][head[i]][col['tcent']] + tlag[i] \
                < channels[0][head[0]][col['tcent']] + tlag[0] - window:
                head[i] += 1
            while tail[i] < len(channels[i]) \
                and channels[i][tail[i]][col['tcent']] + tlag[i] \
                < channels[0][head[0]][col['tcent']] + tlag[0] + window:
                tail[i] += 1

        # run coincidence over all combinations by walking through local triggers, set starting point

        trig[0] = head[0]
        i = 1
        trig[i] = head[i]

        # once we decide to "move on" from the only trigger in the first channel, we are done

        while trig[0] == head[0]:

            # go back to previous channel if we've exhausted the possibilities for a channel

            if trig[i] == tail[i]:
                i -= 1
                trig[i] += 1
            else:

                # check backward for coincidence with everything up to current trigger

                coincident = True
                for j in range(0, i):
                    coincident &= abs(channels[i][trig[i]][col['tcent'
                            ]] + tlag[i]
                            - channels[j][trig[j]][col['tcent']]
                            - tlag[j]) <= window
                if coincident:

                    # we are at last channel (full n-fold coincident event)

                    if i == nchan - 1:

                        # save triggers and move on to next trigger in last channel

                        coinc.append([channels[n][trig[n]] for n in
                                range(0, nchan)])
                        trig[i] += 1
                    else:

                        # stay on current trigger and move on to next channel

                        i += 1
                        trig[i] = head[i]
                else:

                    # not coincident, go to next trigger in current channel

                    trig[i] += 1

        # go to next trigger in first channel

        head[0] += 1
        tail[0] += 1

    # finished looping over triggers in primary channel

    return coinc


# cluster nearby triggers or segments, triggers must be sorted by start time, average values saved
# triggers: list of triggers
# Jan 11, 2008 L. Blackburn

def cluster(triggers, window):
    clusters = []

    # [:] is used to get a copy of the original trigger so that the original is not changed

    cluster = (triggers[0])[:]

    # working with triggers

    if len(cluster) == len(col):
        for trigger in triggers:
            if trigger[col['tstart']] <= cluster[col['tstop']] + window:

                # merge and update cluster properties

                eratio = trigger[col['nenergy']] / cluster[col['nenergy'
                        ]]
                cluster[col['tstart']] = min(cluster[col['tstart']],
                        trigger[col['tstart']])
                cluster[col['tstop']] = max(cluster[col['tstop']],
                        trigger[col['tstop']])
                cluster[col['tcent']] = (cluster[col['tcent']] + eratio
                        * trigger[col['tcent']]) / (1 + eratio)
                cluster[col['fcent']] = (cluster[col['fcent']] + eratio
                        * trigger[col['fcent']]) / (1 + eratio)
                cluster[col['uenergy']] = cluster[col['uenergy']] \
                    + trigger[col['uenergy']]
                cluster[col['nenergy']] = cluster[col['nenergy']] \
                    + trigger[col['nenergy']]
                cluster[col['npix']] = cluster[col['npix']] \
                    + trigger[col['npix']]
                cluster[col['signif']] = cluster[col['signif']] \
                    + trigger[col['signif']]
            else:

                # save cluster and move on

                clusters.append(cluster)
                cluster = trigger[:]

        # write out final cluster

        clusters.append(cluster)
    else:

    # working with segments or something else (only assume start and stop times)

        for segment in triggers:
            if segment[0] <= cluster[1] + window:

                # min() call should never do anything for a properly sorted list

                cluster[0] = min(cluster[0], segment[0])
                cluster[1] = max(cluster[1], segment[1])
            else:
                clusters.append(cluster)
                cluster = segment[:]
        clusters.append(cluster)
    return clusters


# random times from segment list

def randomtimes(ntimes, segments):
    import random
    randomtimes = []
    if isinstance(segments[0], (int, long, float)):
        segments = [segments]
    lt = livetime(segments)
    for i in range(0, ntimes):
        bglt = lt * random.random()
        bgseg = 0
        while bglt > segments[bgseg][1] - segments[bgseg][0]:
            bglt -= segments[bgseg][1] - segments[bgseg][0]
            bgseg += 1
        randomtimes.append(segments[bgseg][0] + bglt)
    randomtimes.sort()
    return randomtimes


# random (poisson) times corresponding to a particular rate

def randomrate(rate, segments):
    if not livetime(segments): ### no time in which to place events
        return []

    import random
    randomtimes = []
    if isinstance(segments[0], (int, long, float)):
        segments = [segments]
    for seg in segments:
        t = seg[0]
        while True:
            wait = random.expovariate(rate)
            t += wait
            if t < seg[-1]:
                randomtimes.append(t)
            else:
                break
    return randomtimes


##@}
