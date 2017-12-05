#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (C) 2013 Reed Essick
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

## \defgroup laldetchar_py_idq_ovl OVL Module
## \ingroup laldetchar_py_idq
## Synopsis
# ~~~
# from laldetchar.idq import ovl
# ~~~
# \author Reed Essick (<reed.essick@ligo.org>)

"""Module that encapsulates all the key OVL scripts as functions .
"""

import os
import sys
import math
import time
import pickle
import shutil
import glob
from laldetchar.idq import event
from laldetchar import git_version

__author__ = 'Reed Essick <reed.essick@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

## \addtogroup laldetchar_py_idq_ovl
# @{

# Encapsulates all the key OVL scripts as functions

# ===================================================================================================
#
#
#                                  define global variables
#
#
# ===================================================================================================
# dictionary for indecies for .stats files
global sD
sD = {
    'livetime': 0,
    '#gwtrg': 1,
    'vchan': 2,
    'vthr': 3,
    'vwin': 4,
    'dsec': 5,
    '#auxtrg': 6,
    'vexp': 7,
    'vact': 8,
    'vsig': 9,
    'useP': 10,
    'eff': 11,
    'dt': 12,
    'eff/dt': 13,
    }

# dictionary for indecies for .vetolist files
global vD
vD = {
    'livetime': 0,
    '#gwtrg': 1,
    'bin': 2,
    'vchan': 3,
    'vthr': 4,
    'vwin': 5,
    'dsec': 6,
    '#auxtrg': 7,
    'vexp': 8,
    'vact': 9,
    'vsig': 10,
    'useP': 11,
    'eff': 12,
    'dt': 13,
    'eff/dt': 14,
    'c_vact': 15,
    'c_dsec': 16,
    'c_auxtrg': 17,
    'metric_exp': 18,
    }


# dictionary for indecies for KW-TRIGGERS 32 second .trg files
global kw_trgD
kw_trgD = {
    'GPSstart': 0,
    'GPSstop': 1,
    'GPScent': 2,
    'fcent': 3,
    'uenergy': 4,
    'nenergy': 5,
    'npix': 6,
    'signif': 7,
    'chan': 8,
    }
kw_trgD['tcent'] = kw_trgD['GPScent']


# ===================================================================================================
#
#
#                          utility functions
#
#
# ===================================================================================================
def effbydt_to_rank(effbydt_exp):
    """ conversion from eff/dt to rank """
    return 1. - math.exp(-effbydt_exp / 100.0)  # a transformend version of eff/dt, mapping it to [0,1]
                                              # normalization by 100 chosen based off of previous experience to make "ranks" reasonable
###
def rank_to_effbydt(rank):
    """ conversion from rank to eff/dt """
    return -100 * math.log(1. - rank)

###
def vsig_to_rank(vsig):
    """ conversion from vsig to rank """
    return 1. - math.exp(-vsig / 100.0)

###
def rank_to_vsig(rank):
    """ conversoin from rank to vsig """
    return -100.0 * math.log(1. - rank)

###
def useP_to_rank(useP):
    """ conversion from useP to rank """
    return 1. - math.exp(- useP / 0.5)

###
def rank_to_useP(rank):
    """ conversion from rank to useP """
    return -0.5 * math.log(1. - rank)

# ================================================
def __check_dirname(dir):
    """ ensures that dir ends in '/' or is an empty string"""
    if dir != '' and dir[-1] != '/':
        dir += '/'
    return dir


#=================================================
def __load_channels(auxdir, analysis_range=False):
    """ loads channels for OVL analysis 
  expect auxdir to be a directory in which auxiliary triggers are stored. Specifically, expect the following directory structure:
    ~/auxdir/GPSstart_GPSstop/channel_name.trg
  """
    auxdir = __check_dirname(auxdir)
    channels = []

    for dir in sorted(os.listdir(auxdir)):
        L = dir.split('_')
        if len(L) == 2:  # must have GPSstart_GPSstop format
            try:
                [start, stop] = [int(i) for i in L]  # GPSstart and GPSstop must be castable as integers
                if not analysis_range or analysis_range[0] <= start \
                    < analysis_range[1] or analysis_range[0] < stop \
                    <= analysis_range[1] or analysis_range[0] <= start \
                    and stop <= analysis_range[1]:  # if analysis_range is given, we select only the directories of interest
                    for l in os.listdir(auxdir + dir):
                        if '.trg' in l:
                            chan = l.split('.trg')[0]
                            if chan not in channels or len(channels) == 0:
                                channels.append(chan)
            except:
                pass

    return channels


#=================================================

def load_trg(
    trgdir,
    analysis_range,
    chans,
    thr,
    trgtype='kw',
    loadtype='full',
    suffix='.trg',
    ):
    """ written to implement trigger discovery and loading in a single function
  should be used to load both veto triggers (from auxdir) and gw triggers (from gwdir)
  """

    if trgtype == 'kw':
        signif_col = event.col_kw['signif']
        tcent_col = event.col_kw['tcent']
    else:
        raise ValueError('unkown trgtype=%s' % trgtype)

  # determine which directories overlap with analysis_range
    dirs = [l for l in os.listdir(trgdir) if '_' in l and '.' not in l
            and '-' not in l]  # and ("H" not in l) and ("L" not in l)]
    target_dirs = []
    for _dir in sorted(dirs):
        [start, stop] = [int(b) for b in _dir.split('_')]  # assume GPSstart_GPSstop
        if start < analysis_range[1] and start >= analysis_range[0] \
            or stop <= analysis_range[1] and stop > analysis_range[0] \
            or start <= analysis_range[0] and stop >= analysis_range[1]:
            target_dirs.append(_dir)

  # iterate over directories and load triggers
    trg = []
    seg = [analysis_range]
    for _dir in sorted(target_dirs):
        _dir = __check_dirname(_dir)
        segfiles = [l for l in os.listdir(trgdir + '/' + _dir) if '.seg'
                     in l and l[0] in [chan[0] for chan in chans]]  # look for segment files
        if len(segfiles):
            seg += event.andsegments([event.load(trgdir + '/' + _dir
                    + '/' + segfiles[0]), [analysis_range]])
        else:
            seg += event.andsegments([[sorted([float(b.strip('/'))
                    for b in _dir.split('_')])], [analysis_range]])
        for chan in chans:
            path = trgdir + '/' + _dir + '/' + chan + suffix
            if os.path.exists(path):
                if loadtype == 'full':
                    trg += event.loadtrg(path, threshold=thr,
                            signif=signif_col)
                elif loadtype == 'veto':
                    trg += event.loadvetotrg(path, threshold=thr,
                            signif=signif_col, tcent=tcent_col)
                else:
                    raise ValueError('unkown loadtype=%s' % loadtype)

    seg = event.fixsegments(seg)
    return (trg, seg)


#=================================================

def round_robin_bin_to_segments(binNo, totBins, __params):
    """ defines segments for a round-robin analysis 
  expect binNo to be an integer (non-negative)
  expect totBins to be an integer > binNo
  expect __params to be a params object
  """

  # slice the analysis range into small segments, each a minute in duration
    min_segs = [[min, min + 60] for min in
                range(__params.analysis_range[0],
                __params.analysis_range[1], 60)]

  # create a list for desired segments
    win_segs = []
    for index in range(len(min_segs) / totBins):
        win_segs.append(min_segs[binNo + index * totBins])
    if binNo <= len(min_segs) % totBins:
        win_segs.append(min_segs[len(min_segs) - 1 - len(min_segs)
                        % totBins + binNo])

    return win_segs

# ===================================================================================================
#
#
#                             params object and handlers
#
#
# ===================================================================================================

class params(object):
    """ an object that represents the params structure. contains all information necessary for ovl to run. """

    required_keys = [
        'analysis_range',
        'auxdir',
        'gwdir',
        'gwchans',
        'gwthr',
        'ifos',
        'gwsets',
        'notused',
        'windows',
        'thresholds',
        'Psigthr',
        'effbydtthr',
        'safety',
        'metric',
        'suffix',
        'trgtype',
        ]

    def __init__(
        self,
        analysis_range,
        auxdir,
        gwdir,
        gwchans,
        gwthr,
        ifos,
        gwsets,
        scisegs=False,
        vetosegs=False,
        channels=False,
        notused=[],
        windows=[0.025, 0.050, 0.100, 0.150, 0.200],
        thresholds=[
            15,
            25,
            30,
            50,
            100,
            200,
            400,
            800,
            1600,
            ],
        Psigthr=10 ** -5,
        effbydtthr=3.0,
        safety='None',
        metric='eff/dt',
        suffix='.trg',
        trgtype='kw',
        ):
        """ constructs a params object. "metric" can be any of the columns defined in ovl.sD, but should "eff/dt", "vsig" (poisson probability), or "useP" (use percentage)"""

    # comments after a variable are examples of what that variable should be
    # these fileds do not have ``standard'' values and must be supplied
        self.analysis_range = analysis_range  # [959126400, 959731200]
        self.auxdir = auxdir  # '/archive/home/lindy/public_html/triggers/s6-merged/'
        self.gwdir = gwdir  # '/archive/home/lindy/public_html/auxmvc/test3/'
        self.gwchans = gwchans  # ['chan', 'chan',...]
        self.gwthr = float(gwthr)  # 35
        self.ifos = ifos  # ['L1', ifo',...]
        self.gwsets = gwsets  # ['kwl1-35', 'gwset', ...]

    # these fields have ``standard'' values and do not have to be supplied
        self.windows = sorted(windows)
        self.thresholds = sorted(thresholds)
        self.Psigthr = Psigthr
        self.effbydtthr = effbydtthr
        if safety == 'None':
            self.safety = None
        else:
            self.safety = safety
        self.notused = notused + gwchans  # notused = 'PEM-EY_MAGY PEM-EY_MAGZ PEM-LVEA_MAGX PEM-RADIO_LVEA ASC-QPDX_DC ASC-QPDX_P ASC-QPDX_Y LSC-DARM_CTRL LSC-DARM_ERR LSC-ETMX_EXC_DAQ OMC-PD_SUM_OUT_DAQ OMC-QPD1_Y_OUT_DAQ OMC-QPD3_SUM_IN1_DAQ OMC-QPD4_SUM_IN1_DAQ OMC-NULLSTREAM_OUT_DAQ SUS-ETMX_SUSPIT_IN SUS-ETMX_SUSPOS_IN'.split()
        self.metric = metric  # "effbydt", "Ppoisson", else?
        self.suffix = suffix
        self.trgtype = trgtype

    # these parameters are only included if they are supplied
        if scisegs:
            self.scisegs = scisegs  # ['/archive/home/lindy/public_html/triggers/s6-segments/l1_science.seg']
        if vetosegs:
            self.vetosegs = vetosegs  # ['/archive/home/lindy/public_html/auxmvc/test3/burstdq/L1-VETOTIME_CAT2-959126400-604800.xml']
        if channels:
            self.channels = channels

    def check(self):
        """ checks to make sure all required fields are supplied """
        keys = sorted(vars(self).keys())
        for key in self.required_keys:
            if key not in keys:
                print key
                return False
            else:
                pass
        else:
            return True


# =================================================

def load_params(params_filename, params_dir=''):
    """ load a params file 
  expect params_filename to be a string corresponding to an ASCII params file
  expect params_dir to be the directory containing params_filename
  """

    params_dir = __check_dirname(params_dir)
    f = open(params_dir + params_filename, 'r')
    dic = dict()
    for line in f:
        if line[0] != '#':
            line = [b.strip() for b in line.strip().split('=')]
            if line[0] == 'analysis_range':
                line[1] = [int(l.strip(',')) for l in line[1].strip('['
                           ).strip(']').split()]
            elif line[0] in ['windows', 'thresholds']:
                line[1] = [float(l.strip(',')) for l in
                           line[1].strip('[').strip(']').split()]
            elif line[0] in [
                'gwchans',
                'notused',
                'gwsets',
                'ifos',
                'scisegs',
                'vetosegs',
                ]:
                line[1] = [str(l.strip(',').strip("\'")) for l in
                           line[1].strip('[').strip(']').split()]
            elif line[0] in ['effbydtthr', 'Psigthr', 'gwthr']:
                line[1] = float(line[1])
            else:
                pass  # line[1] should remain a string

            dic[line[0]] = line[1]

    if not dic.has_key('scisegs'):
        dic['scisegs'] = False
    if not dic.has_key('vetosegs'):
        dic['vetosegs'] = False
    if not dic.has_key('channels'):
        dic['channels'] = False

    new_params = params(
        dic['analysis_range'],
        dic['auxdir'],
        dic['gwdir'],
        dic['gwchans'],
        dic['gwthr'],
        dic['ifos'],
        dic['gwsets'],
        notused=dic['notused'],
        windows=dic['windows'],
        thresholds=dic['thresholds'],
        Psigthr=dic['Psigthr'],
        effbydtthr=dic['effbydtthr'],
        safety=dic['safety'],
        scisegs=dic['scisegs'],
        vetosegs=dic['vetosegs'],
        channels=dic['channels'],
        metric=dic['metric'],
        )

    if new_params.check():
        return new_params
    else:
        raise KeyError('not all required keys present in %s'
                       % params_filename)


# =================================================

def write_params(__params, params_filename, params_dir=''):
    """ write a params file 
  expect __params to be a params object
  expect params_filename to be the name to be assigned to the ASCII version of __params
  expect params_dir to be the directory into which params_filename will be written
  """
    params_dir = __check_dirname(params_dir)
    dic = vars(__params)
    f = open(params_dir + params_filename, 'w')
    for key in dic.keys():
        print >> f, key + ' = ' + str(dic[key])
    f.close()
    return params


# ===================================================================================================
#
#
#                                 MAIN OVL FUNCTIONS
#
#
# ===================================================================================================

def redundancies(__params, output_dir="./", verbose=False, write_channels=False):
    """
    computes the intersection and overlap of vetosegments for each possible configuration. 
    This should contain all information necessary to determine which channels are redundant.
    """
    output_dir = __check_dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # path to auxiliary channel KW triggers
    auxdir = __params.auxdir

    # path to gw channel KW triggers
    gwdir = __params.gwdir

    # load list of auxiliary channels
    if vars(__params).has_key('channels'):
        channels = event.loadlist(__params.channels)
    else:
        channels = __load_channels(auxdir,
                                   analysis_range=__params.analysis_range)

    # we need the same channels list as for recalculate.py so we can associate numbers with channels
    for remove in __params.notused:
        channels = [channel for channel in channels
                    if channel.find(remove) == -1]

    # create a channels file
    if write_channels:
        ch_file = open(output_dir + 'channels.txt', 'w')
        for chan in channels:
            print >> ch_file, chan
        ch_file.close()

    ### segments for this analysis range
    gwseg = [__params.analysis_range]

    # load (and clean) a list of segments for science runs
    scisegs = dict()
    if vars(__params).has_key('scisegs'):
        for (i, ifo) in enumerate(__params.ifos):
            segs = event.load(__params.scisegs[i])

            # segwizard format
            if len(segs) > 0 and len(segs[0]) > 2:
                # assume that the biggest number will be the stop time of the first segment
                index = segs[0].index(max(segs[0]))
                segs = [seg[index - 1:index + 1] for seg in segs]
            gwseg = event.andsegments([gwseg, segs])  # intersect segs with analysis range
    
    # load (and clean) a list of veto segements
    vetosegs = dict()
    if vars(__params).has_key('vetosegs'):
        for (i, ifo) in enumerate(__params.ifos):
            segs = event.load(params.vetosegs[i])
            if len(segs) > 0 and len(segs[0]) > 2:
                index = (segs[0], index(max(segs[0])))
                segs = [seg[index - 1:index + 1] for seg in segs]
            gwseg = event.removesegments( gwseg, segs )

    livetime = event.livetime(gwseg) ### total amount of time in the analysis

    # create a list of all veto triggers to avoid loading multiple times
    allvtrg = [0 for i in channels]
    allvseg = [[[0 for i in __params.thresholds] for j in __params.windows] for k in channels]

    nthr = len(__params.thresholds)
    nwin = len(__params.windows)
    nconfig = nthr*nwin

    ### iterate over channels, creating segments and comparing them
    filename = "redundancies.txt"
    f = open(output_dir + filename, "w")

    print >>f, livetime
    if verbose:
        print livetime

    s = ""
    for chan in channels:
        for vwin in __params.windows:
            for vthr in __params.thresholds:
                s += " (%s,%f,%f)"%(chan, vwin, vthr)
    print >>f, s
    if verbose:
        print s

    event.col = event.col_veto
    for indA, chanA in enumerate(channels):
        if allvtrg[indA] == 0: ### load chanA triggers
            allvtrg[indA], _ = load_trg( auxdir, __params.analysis_range, [chanA], __params.thresholds[0], trgtype=__params.trgtype, loadtype='veto', suffix=__params.suffix)
  
        for indvwinA, vwinA in enumerate(__params.windows):
            for indvthrA, vthrA in enumerate(__params.thresholds):

                if allvseg[indA][indvwinA][indvthrA] == 0:
                    vetosegA = event.vetosegs(allvtrg[indA], vwinA, vthrA) ### generate chanA segments
                    vetosegA = event.andsegments([vetosegA, gwseg])              
                    livetimeA = event.livetime(vetosegA)
                    allvseg[indA][indvwinA][indvthrA] = (vetosegA, livetimeA)
                else:
                    vetosegA, livetimeA = allvseg[indA][indvwinA][indvthrA]

                s = ""

                for indB, chanB in enumerate(channels[indA+1:]):
                    indB += indA+1
                    if allvtrg[indB] == 0: ### load chanB triggers
                        allvtrg[indB], _ = load_trg( auxdir, __params.analysis_range, [chanB], __params.thresholds[0], trgtype=__params.trgtype, loadtype='veto', suffix=__params.suffix)
                    for indvwinB, vwinB in enumerate(__params.windows):
                        for indvthrB, vthrB in enumerate(__params.thresholds):

                            if allvseg[indB][indvwinB][indvthrB] == 0:
                                vetosegB = event.vetosegs(allvtrg[indB], vwinB, vthrB) ### generate chanB segments
                                vetosegB = event.andsegments([vetosegB, gwseg])
                                livetimeB = event.livetime(vetosegB)
                                allvseg[indB][indvwinB][indvthrB] = (vetosegB, livetimeB)
                            else:
                                vetosegB, livetimeB = allvseg[indB][indvwinB][indvthrB]

                            ### compare segments
                            i = event.livetime( event.andsegments([vetosegA, vetosegB]) )
                            u = event.livetime( event.andsegments([vetosegA, vetosegB]) )

                            if livetimeA and livetimeB:
                                stat1 = i*livetime/(livetimeA*livetimeB)
                            else:
                                stat1 = 1
                            if u:
                                stat2 = i/u
                            else:
                                stat2 = 0

                            s += "(%9.3f,%9.3f,%9.3f,%9.3f,%5.1f,%1.4f) "%(livetimeA, livetimeB, i, u, stat1, stat2)

                print >> f, s
                if verbose:
                    print s

                allvseg[indA][indvwinA][indvthrA] = 0 ### forget about segments for this config
                            
        allvtrg[indA] = 0 ### forget about triggers we don't need anymore

    return filename

#=================================================

def recalculate(
    run,
    incremental,
    __params,
    output_dir='./',
    source_dir='./',
    eval=False,
    binNo=False,
    totBins=False,
    track=False,
    verbose=False,
    write_channels=False,
    ):
    """ an encapsulation of recalculate.py and recalculate_eval.py 
  run is the iteration number for the algorithm. expect a non-negative integer
  incremental is a limit on the number of independently-applied vconfigs. expect a non-negative integer
  __params should be a params object
  binNo, totBins set up round-robin analysis. expect both to be non-negative integers with binNo < totBins
  eval is a flag telling the algorithm to always remove data (and to use the analysis data set for round-robin algorithms)
  track is a flag that creates a trackfile.pickle
  verbose is a flag that causes progess statements to be printed (including a duplicate of the .stats file)
  """

  # check for consistent input, specifically round_robin binning
    source_dir = __check_dirname(source_dir)
    output_dir = __check_dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if binNo or totBins:
        if binNo >= totBins:
            sys.exit('binNo not consistent with totBins')
        else:
            round_robin = True
    else:
        round_robin = False
    if verbose:
        print 'round_robin = ' + repr(round_robin)
        print 'track = ' + repr(track)
        print 'eval = ' + repr(eval)

  # load thresholds
    sigthr = -math.log(__params.Psigthr)
    effbydtthr = __params.effbydtthr

  # path to auxiliary channel KW triggers
    auxdir = __params.auxdir

  # path to gw channel KW triggers
    gwdir = __params.gwdir

  # load list of auxiliary channels
    if vars(__params).has_key('channels'):
        channels = event.loadlist(__params.channels)
    else:
        channels = __load_channels(auxdir,
                                   analysis_range=__params.analysis_range)

  # we need the same channels list as for recalculate.py so we can associate numbers with channels
    for remove in __params.notused:
        channels = [channel for channel in channels
                    if channel.find(remove) == -1]

  # create a channels file
    if write_channels:
        ch_file = open(output_dir + 'channels.txt', 'w')
        for chan in channels:
            print >> ch_file, chan
        ch_file.close()

  # create a checkeroard windo array within the analysis range
    if round_robin:
        win_segs = round_robin_bin_to_segments(binNo, totBins, __params)

  # load (and clean) a list of segments for science runs
    scisegs = dict()
    if vars(__params).has_key('scisegs'):
        for (i, ifo) in enumerate(__params.ifos):
            segs = event.load(__params.scisegs[i])

      # segwizard format
            if len(segs) > 0 and len(segs[0]) > 2:

        # assume that the biggest number will be the stop time of the first segment
                index = segs[0].index(max(segs[0]))
                segs = [seg[index - 1:index + 1] for seg in segs]
            scisegs[ifo] = \
                event.andsegments([[__params.analysis_range], segs])  # intersect segs with analysis range
    else:
        for (i, ifo) in enumerate(__params.ifos):
            scisegs[ifo] = [__params.analysis_range]

  # window science segments for round-robin analysis
    if round_robin:
        for ifo in __params.ifos:
            if not eval:  # use all segments except win_segs
                scisegs[ifo] = event.removesegments([scisegs[ifo],
                        event.andsegments([win_segs, scisegs[ifo]])])
            if eval:  # use only win_segs
                scisegs[ifo] = event.andsegments([win_segs,
                        scisegs[ifo]])

  # load (and clean) a list of veto segements
    vetosegs = dict()
    if vars(__params).has_key('vetosegs'):
        for (i, ifo) in enumerate(__params.ifos):
            segs = event.load(params.vetosegs[i])
            if len(segs) > 0 and len(segs[0]) > 2:
                index = (segs[0], index(max(segs[0])))
                segs = [seg[index - 1:index + 1] for seg in segs]
            vetosegs[ifo] = \
                event.andsegments([[__params.analysis_range], segs])

  # create a list of all veto triggers to avoid loading multiple times
    allvtrg = [0 for i in range(len(channels))]

    statsfiles = []  # a list of the completed stats files
    trackfiles = []  # list of completed track files

  # iterate analysi over all sets of GW triggers listed in __params.gwsets
    for gwset in __params.gwsets:
        gwset += '-' + str(__params.gwthr)

    # main output file
        filename = gwset + '.stats.' + repr(run)
        if eval:
            filename += '.eval'
        if round_robin:
            if not eval:
                filename += '.training.bin_' + repr(binNo) + '_out_of_' \
                    + repr(totBins)
            else:
                filename += '.analysis.bin_' + repr(binNo) + '_out_of_' \
                    + repr(totBins)
        f = open(output_dir + filename, 'w')

    # tracking file
        if track:
            trackfilename = gwset + '.track.' + repr(run)
            if eval:
                trackfilename += '.eval'
            if round_robin:
                if not eval:
                    trackfilename += '.training.bin_' + repr(binNo) \
                        + '_out_of_' + repr(totBins)
                else:
                    trackfilename += '.analysis.bin_' + repr(binNo) \
                        + '_out_of_' + repr(totBins)
            trackfilename = output_dir + trackfilename + '.pickle'
            tfp = open(trackfilename, 'w')

        [method, type] = gwset.split('-')

        if __params.trgtype == 'kw':
            event.col = event.col_kw
        else:
            raise ValueError('unkown trgtype=%s' % __params.trgtype)

    # get the appropriate directories (which are bounded by days and do not necessarily line up with params.analysis_range)
        (gwtrg, gwseg) = load_trg(
            gwdir,
            __params.analysis_range,
            __params.gwchans,
            __params.gwthr,
            trgtype=__params.trgtype,
            loadtype='full',
            suffix=__params.suffix,
            )

    # exclude GW triggers which fall outside of analysis segments
        for index in event.ifo[method].values():
            gwtrg = event.include(gwtrg, gwseg, tcent=index)

    # exclude GW triggers which fall outside of science segments
        if vars(__params).has_key('scisegs'):
            for ifo in event.ifo[method].keys():
                gwtrg = event.include(gwtrg, scisegs[ifo],
                        tcent=event.ifo[method][ifo])
                gwseg = event.andsegments([gwseg, scisegs[ifo]])

    # apply prior veto segments to prefilter events
        if vars(__params).has_key('vetosegs'):
            for ifo in event.ifo[method].keys():
                gwtrg = event.exclude(gwtrg, vetosegs[ifo],
                        tcent=event.ifo[method][ifo])
                gwseg = event.removesegments(gwseg, vetosegs[ifo])

        if run == 0:
            if __params.safety != None:
                safethr = dict(event.loadstringtable(params.safety))
            else:
                safethr = dict()

      # create a channel list including all possible combinations of channel, threshold and window
            stats = [[
                0,
                0,
                vchan,
                vthr,
                vwin,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                ] for vchan in range(len(channels)) for vthr in
                    __params.thresholds for vwin in __params.windows]

      # remove veto channels wich aren't at the site that we are interested in
            stats = [line for line in stats if len([key for key in
                     event.ifo[method].keys() if key[0]
                     == channels[line[sD['vchan']]][0]]) > 0]

      # remove veto parameters whcih are determined to be unsafe
            stats = [line for line in stats
                     if not safethr.has_key(channels[line[sD['vchan'
                     ]]]) or float(safethr[channels[line[sD['vchan'
                     ]]]]) <= line[sD['vthr']]]

      # never calculate incremental statistics for the first pass
            incremental = 0
        else:

          # load data from previous run
            o_filename = gwset + '.stats.' + repr(run - 1)
            if round_robin:
                o_filename += '.training.bin_' + repr(binNo) \
                    + '_out_of_' + repr(totBins)
            stats = event.load(source_dir + o_filename)

      # keep only those configurations that perform well enough
            stats = [line for line in stats if line[sD['vsig']]
                     > sigthr and line[sD['eff/dt']] > effbydtthr]

      # reverse sort by threshold so highest thresholds are chosen first
            stats.sort(key=lambda line: line[sD['vthr']], reverse=True)

      # sort by metric : "eff/dt", "vsig", else?
#      stats.sort(key=lambda line: line[sD['eff/dt']], reverse=True)
            stats.sort(key=lambda line: line[sD[__params.metric]],
                       reverse=True)

    # print heading and identifiers into main output file
        print >> f, '# analysis_range = ' \
            + repr(__params.analysis_range)
        print >> f, '# run = ' + repr(run)
        print >> f, '# incremental = ' + repr(incremental)
        if round_robin:
            print >> f, '# binNo = ' + repr(binNo)
            print >> f, '# totBins = ' + repr(totBins)
        print >> f, \
            '# %10s %7s %5s %4s %7s %11s %8s %6s %5s %8s %9s %9s %9s %9s %9s %9s' \
            % (
            'livetime',
            '#gwtrg',
            'vchan',
            'vthr',
            'vwin',
            'dsec',
            '#auxtrg',
            'vexp',
            'vact',
            'vsig',
            'useP',
            'veff',
            'dt',
            'eff/dt',
            'time(usr)',
            'time(real)',
            )
        t0 = time.clock()
        t1 = time.time()

        if track:  # build headers in tfp
            tfparray = [0] * (len(stats) + 1)
            tfparray[0] = [
                'livetime',
                '#gwtrg',
                'vchan',
                'vthr',
                'vwin',
                'dsec',
                '#auxtrg',
                'vact',
                'vsig',
                'gwtrg_vetoed',
                ]

        if verbose:  # print column headings, etc
            print 'analysis_range = ' + repr(__params.analysis_range)
            print 'run = ' + repr(run)
            print 'incremental = ' + repr(incremental)
            if round_robin:
                print 'binNo = ' + repr(binNo)
                print 'totBins = ' + repr(totBins)
            print '%12s %7s %5s %4s %7s %11s %8s %6s %5s %8s %9s %9s %9s %9s %9s %9s' \
                % (
                'livetime',
                '#gwtrg',
                'vchan',
                'vthr',
                'vwin',
                'dsec',
                '#auxtrg',
                'vexp',
                'vact',
                'vsig',
                'useP',
                'eff',
                'dt',
                'eff/dt',
                'time(usr)',
                'time(real)',
                )

    # ###########
    # BEGIN the analysis
    # ###########
        for lineidx in xrange(len(stats)):
            line = stats[lineidx]

            ngwtrg = len(gwtrg)  # '#gwtrg'
            livetime = event.livetime(gwseg)
            if ngwtrg <= 0 or livetime <= 0:  # stop if there are no more gwtrg
                if verbose:
                    print 'ngwtrg = ' + str(ngwtrg)
                    print 'livetime = ' + str(livetime)
                    print 'gwseg = ' + str(gwseg)
                    print 'line = ' + str(line)
                    print 'no more gwtrg or livetime'
                break
            gwrate = float(ngwtrg) / float(livetime)

      # only load the veto channel once and save it for later use
            (vchan, vthr, vwin) = (int(line[sD['vchan']]),
                                   line[sD['vthr']], line[sD['vwin']])

            if allvtrg[vchan] == 0:
                allvtrg[vchan] = []

        # find the bounds so that we pull the correct directories
                (allvtrg[vchan], _) = load_trg(
                    auxdir,
                    __params.analysis_range,
                    [channels[vchan]],
                    __params.thresholds[0],
                    trgtype=__params.trgtype,
                    loadtype='veto',
                    suffix=__params.suffix,
                    )

            vtrg = allvtrg[vchan]  # full set of veto triggers
            nauxtrg = len(event.include([trg for trg in vtrg
                          if trg[event.col_veto['signif']] >= vthr],
                          gwseg, tcent=event.col_veto['tcent']))  # number of auxiliary triggers
            event.col = event.col_veto

            vetoseg = event.vetosegs(vtrg, vwin, vthr)  # builds veto segments based on configuration parameters
            vetoseg = event.andsegments(vetoseg, gwseg)  # only count the intersection

            deadsecs = event.livetime(vetoseg)
            deadfrac = float(deadsecs) / float(livetime)
            vexpected = deadsecs * gwrate

      # Identify coincident gwtrgs
            gwtrg_postveto = gwtrg[:]
            if track:
                for ifo in [key for key in event.ifo[method].keys()
                            if key == (channels[vchan])[:2]
                            or channels[vchan][1] == '0' and key[0]
                            == channels[vchan][0]]:
                    [gwtrg_vetoed, gwtrg_postveto] = \
                        event.includeexclude(gwtrg_postveto, vetoseg,
                            tcent=event.ifo[method][ifo])
            else:
                for ifo in [key for key in event.ifo[method].keys()
                            if key == (channels[vchan])[:2]
                            or channels[vchan][1] == '0' and key[0]
                            == channels[vchan][0]]:
                    gwtrg_postveto = event.exclude(gwtrg_postveto,
                            vetoseg, tcent=event.ifo[method][ifo])

            gwseg_postveto = event.removesegments(gwseg, vetoseg)

            vactual = ngwtrg - len(gwtrg_postveto)

            if nauxtrg > 0:
                useP = 1. * vactual / nauxtrg
            else:
                useP = 0.

            vsig = 0. - _gammpln(vactual, vexpected)
            veff = float(vactual) / float(ngwtrg)
            if deadfrac == 0:
                effbydt = 0.
            else:
                effbydt = veff / deadfrac

      # print results to main output file
            print >> f, \
                '%12.3f %7d %5d %4d %7.3f %11.3f %8d %6.2f %5d %8.2f %9.6f %9.5f %9.5f %9.4f %9.2f %9.2f' \
                % (
                livetime,
                ngwtrg,
                vchan,
                vthr,
                vwin,
                deadsecs,
                nauxtrg,
                vexpected,
                vactual,
                vsig,
                useP,
                veff,
                deadfrac,
                effbydt,
                time.clock() - t0,
                time.time() - t1,
                )

            if track:
                if len(gwtrg_vetoed) == 0:
                    gwtrg_vetoed = ['NONE']
                tfparray[lineidx + 1] = [
                    livetime,
                    ngwtrg,
                    vchan,
                    vthr,
                    vwin,
                    deadsecs,
                    nauxtrg,
                    vactual,
                    vsig,
                    gwtrg_vetoed,
                    ]

            if verbose:
                print '%12.3f %7d %5d %4d %7.3f %11.3f %8d %6.2f %5d %8.2f %9.6f %9.5f %9.5f %9.4f %9.2f %9.2f' \
                    % (
                    livetime,
                    ngwtrg,
                    vchan,
                    vthr,
                    vwin,
                    deadsecs,
                    nauxtrg,
                    vexpected,
                    vactual,
                    vsig,
                    useP,
                    veff,
                    deadfrac,
                    effbydt,
                    time.clock() - t0,
                    time.time() - t1,
                    )

            if incremental > 0:  # check whether we remove data vetoed by this configuration
                if eval or vsig > sigthr and effbydt > effbydtthr:
                    gwtrg = gwtrg_postveto
                    gwseg = gwseg_postveto
                incremental -= 1

            for i in range(lineidx + 1, len(stats)):  # check if we need to remember the triggers for vchan
                if int(stats[i][sD['vchan']]) == vchan:
                    break
            else:
                allvtrg[vchan] = 0

            f.flush()

        f.close()
        statsfiles.append(filename)
        if track:
            pickle.dump(tfparray, tfp)
            tfp.close()
            trackfiles.append(trackfilename)

    return (statsfiles, trackfiles)

#==================================================

def optimize(
    run,
    __params,
    output_dir='./',
    source_dir='./',
    binNo=False,
    totBins=False,
    track=False,
    verbose=False,
    write_channels=False,
    ):
    """ given a set of configurations, this finds the optimal order in which to apply them. 
    Computationally expensive, so recalculate() should be called first to reduce the number of configurations. 
    This should be called between recalculate(eval=False) and vetolist_eval() during training
    """

    eval=True ### this will always be the case for ovl.optimize, and saves me editing work below...

  # check for consistent input, specifically round_robin binning
    source_dir = __check_dirname(source_dir)
    output_dir = __check_dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if binNo or totBins:
        if binNo >= totBins:
            sys.exit('binNo not consistent with totBins')
        else:
            round_robin = True
    else:
        round_robin = False
    if verbose:
        print 'round_robin = ' + repr(round_robin)
        print 'track = ' + repr(track)
        print 'eval = ' + repr(eval)
        print 'optimize = True'

  # load thresholds
    sigthr = -math.log(__params.Psigthr)
    effbydtthr = __params.effbydtthr

  # path to auxiliary channel KW triggers
    auxdir = __params.auxdir

  # path to gw channel KW triggers
    gwdir = __params.gwdir

  # load list of auxiliary channels
    if vars(__params).has_key('channels'):
        channels = event.loadlist(__params.channels)
    else:
        channels = __load_channels(auxdir,
                                   analysis_range=__params.analysis_range)

  # we need the same channels list as for recalculate.py so we can associate numbers with channels
    for remove in __params.notused:
        channels = [channel for channel in channels
                    if channel.find(remove) == -1]

  # create a channels file
    if write_channels:
        ch_file = open(output_dir + 'channels.txt', 'w')
        for chan in channels:
            print >> ch_file, chan
        ch_file.close()

  # create a checkeroard windo array within the analysis range
    if round_robin:
        win_segs = round_robin_bin_to_segments(binNo, totBins, __params)

  # load (and clean) a list of segments for science runs
    scisegs = dict()
    if vars(__params).has_key('scisegs'):
        for (i, ifo) in enumerate(__params.ifos):
            segs = event.load(__params.scisegs[i])

      # segwizard format
            if len(segs) > 0 and len(segs[0]) > 2:

        # assume that the biggest number will be the stop time of the first segment
                index = segs[0].index(max(segs[0]))
                segs = [seg[index - 1:index + 1] for seg in segs]
            scisegs[ifo] = \
                event.andsegments([[__params.analysis_range], segs])  # intersect segs with analysis range
    else:
        for (i, ifo) in enumerate(__params.ifos):
            scisegs[ifo] = [__params.analysis_range]

  # window science segments for round-robin analysis
    if round_robin:
        for ifo in __params.ifos:
            if not eval:  # use all segments except win_segs
                scisegs[ifo] = event.removesegments([scisegs[ifo],
                        event.andsegments([win_segs, scisegs[ifo]])])
            if eval:  # use only win_segs
                scisegs[ifo] = event.andsegments([win_segs,
                        scisegs[ifo]])

  # load (and clean) a list of veto segements
    vetosegs = dict()
    if vars(__params).has_key('vetosegs'):
        for (i, ifo) in enumerate(__params.ifos):
            segs = event.load(params.vetosegs[i])
            if len(segs) > 0 and len(segs[0]) > 2:
                index = (segs[0], index(max(segs[0])))
                segs = [seg[index - 1:index + 1] for seg in segs]
            vetosegs[ifo] = \
                event.andsegments([[__params.analysis_range], segs])

  # create a list of all veto triggers to avoid loading multiple times
    allvtrg = [0 for i in range(len(channels))]

    statsfiles = []  # a list of the completed stats files
    trackfiles = []  # list of completed track files

  # iterate analysi over all sets of GW triggers listed in __params.gwsets
    for gwset in __params.gwsets:
        gwset += '-' + str(__params.gwthr)

    # main output file
        filename = gwset + '.stats.' + repr(run)
        if eval:
            filename += '.eval'
        if round_robin:
            if not eval:
                filename += '.training.bin_' + repr(binNo) + '_out_of_' \
                    + repr(totBins)
            else:
                filename += '.analysis.bin_' + repr(binNo) + '_out_of_' \
                    + repr(totBins)
        f = open(output_dir + filename, 'w')

    # tracking file
        if track:
            trackfilename = gwset + '.track.' + repr(run)
            if eval:
                trackfilename += '.eval'
            if round_robin:
                if not eval:
                    trackfilename += '.training.bin_' + repr(binNo) \
                        + '_out_of_' + repr(totBins)
                else:
                    trackfilename += '.analysis.bin_' + repr(binNo) \
                        + '_out_of_' + repr(totBins)
            trackfilename = output_dir + trackfilename + '.pickle'
            tfp = open(trackfilename, 'w')

        [method, type] = gwset.split('-')

        if __params.trgtype == 'kw':
            event.col = event.col_kw
        else:
            raise ValueError('unkown trgtype=%s' % __params.trgtype)

    # get the appropriate directories (which are bounded by days and do not necessarily line up with params.analysis_range)
        (gwtrg, gwseg) = load_trg(
            gwdir,
            __params.analysis_range,
            __params.gwchans,
            __params.gwthr,
            trgtype=__params.trgtype,
            loadtype='full',
            suffix=__params.suffix,
            )

    # exclude GW triggers which fall outside of analysis segments
        for index in event.ifo[method].values():
            gwtrg = event.include(gwtrg, gwseg, tcent=index)

    # exclude GW triggers which fall outside of science segments
        if vars(__params).has_key('scisegs'):
            for ifo in event.ifo[method].keys():
                gwtrg = event.include(gwtrg, scisegs[ifo],
                        tcent=event.ifo[method][ifo])
                gwseg = event.andsegments([gwseg, scisegs[ifo]])

    # apply prior veto segments to prefilter events
        if vars(__params).has_key('vetosegs'):
            for ifo in event.ifo[method].keys():
                gwtrg = event.exclude(gwtrg, vetosegs[ifo],
                        tcent=event.ifo[method][ifo])
                gwseg = event.removesegments(gwseg, vetosegs[ifo])

        if run == 0:
            if __params.safety != None:
                safethr = dict(event.loadstringtable(params.safety))
            else:
                safethr = dict()

      # create a channel list including all possible combinations of channel, threshold and window
            stats = [[
                0,
                0,
                vchan,
                vthr,
                vwin,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                ] for vchan in range(len(channels)) for vthr in
                    __params.thresholds for vwin in __params.windows]

      # remove veto channels wich aren't at the site that we are interested in
            stats = [line for line in stats if len([key for key in
                     event.ifo[method].keys() if key[0]
                     == channels[line[sD['vchan']]][0]]) > 0]

      # remove veto parameters whcih are determined to be unsafe
            stats = [line for line in stats
                     if not safethr.has_key(channels[line[sD['vchan'
                     ]]]) or float(safethr[channels[line[sD['vchan'
                     ]]]]) <= line[sD['vthr']]]

        else:

          # load data from previous run
            o_filename = gwset + '.stats.' + repr(run - 1)
            if round_robin:
                o_filename += '.training.bin_' + repr(binNo) \
                    + '_out_of_' + repr(totBins)
            stats = event.load(source_dir + o_filename)

      # keep only those configurations that perform well enough
            stats = [line for line in stats if line[sD['vsig']]
                     > sigthr and line[sD['eff/dt']] > effbydtthr]

      # reverse sort by threshold so highest thresholds are chosen first
            stats.sort(key=lambda line: line[sD['vthr']], reverse=True)

      # sort by metric : "eff/dt", "vsig", else?
#      stats.sort(key=lambda line: line[sD['eff/dt']], reverse=True)
            stats.sort(key=lambda line: line[sD[__params.metric]],
                       reverse=True)

    # print heading and identifiers into main output file
        print >> f, '# analysis_range = ' \
            + repr(__params.analysis_range)
        print >> f, '# run = ' + repr(run)
        print >> f, '# incremental = optimize'
        if round_robin:
            print >> f, '# binNo = ' + repr(binNo)
            print >> f, '# totBins = ' + repr(totBins)
        print >> f, \
            '# %10s %7s %5s %4s %7s %11s %8s %6s %5s %8s %9s %9s %9s %9s %9s %9s' \
            % (
            'livetime',
            '#gwtrg',
            'vchan',
            'vthr',
            'vwin',
            'dsec',
            '#auxtrg',
            'vexp',
            'vact',
            'vsig',
            'useP',
            'veff',
            'dt',
            'eff/dt',
            'time(usr)',
            'time(real)',
            )
        t0 = time.clock()
        t1 = time.time()

        if track:  # build headers in tfp
            tfparray = [0] * (len(stats) + 1)
            tfparray[0] = [
                'livetime',
                '#gwtrg',
                'vchan',
                'vthr',
                'vwin',
                'dsec',
                '#auxtrg',
                'vact',
                'vsig',
                'gwtrg_vetoed',
                ]
        if verbose:  # print column headings, etc
            print 'analysis_range = ' + repr(__params.analysis_range)
            print 'run = ' + repr(run)
            print 'incremental = optimize'
            if round_robin:
                print 'binNo = ' + repr(binNo)
                print 'totBins = ' + repr(totBins)
            print '%12s %7s %5s %4s %7s %11s %8s %6s %5s %8s %9s %9s %9s %9s %9s %9s' \
                % (
                'livetime',
                '#gwtrg',
                'vchan',
                'vthr',
                'vwin',
                'dsec',
                '#auxtrg',
                'vexp',
                'vact',
                'vsig',
                'useP',
                'eff',
                'dt',
                'eff/dt',
                'time(usr)',
                'time(real)',
                )

    ############
    # BEGIN the analysis
    ############

        lineidx = 0 ### used for track arrays
        while len(stats): ### iterate until we're done with all configurations
            ngwtrg = len(gwtrg)  # '#gwtrg'
            livetime = event.livetime(gwseg)
            if ngwtrg <= 0 or livetime <= 0:  # stop if there are no more gwtrg
                if verbose:
                    print 'ngwtrg = ' + str(ngwtrg)
                    print 'livetime = ' + str(livetime)
                    print 'gwseg = ' + str(gwseg)
                    print 'line = ' + str(line)
                    print 'no more gwtrg or livetime'
                break
            gwrate = float(ngwtrg) / float(livetime)

            ### iterate over remaining veto configurations and find the best one
            result = None
            track_result = None
            best_metric = -10000
            best_lineid = None
            best_gwtrg = None
            best_gwseg = None
            for lineid, line in enumerate(stats):
                vchan, vthr, vwin = int(line[sD['vchan']]), line[sD['vthr']], line[sD['vwin']]

                if allvtrg[vchan] == 0:
                    allvtrg[vchan] = []

                    # find the bounds so that we pull the correct directories
                    (allvtrg[vchan], _) = load_trg(
                        auxdir,
                        __params.analysis_range,
                        [channels[vchan]],
                        __params.thresholds[0],
                        trgtype=__params.trgtype,
                        loadtype='veto',
                        suffix=__params.suffix,
                        )

                vtrg = allvtrg[vchan]  # full set of veto triggers
                nauxtrg = len(event.include([trg for trg in vtrg
                              if trg[event.col_veto['signif']] >= vthr],
                              gwseg, tcent=event.col_veto['tcent']))  # number of auxiliary triggers
                event.col = event.col_veto

                vetoseg = event.vetosegs(vtrg, vwin, vthr)  # builds veto segments based on configuration parameters
                vetoseg = event.andsegments(vetoseg, gwseg)  # only count the intersection

                deadsecs = event.livetime(vetoseg)
                deadfrac = float(deadsecs) / float(livetime)
                vexpected = deadsecs * gwrate

                # Identify coincident gwtrgs
                gwtrg_postveto = gwtrg[:]
                if track:
                    for ifo in [key for key in event.ifo[method].keys()
                                if key == (channels[vchan])[:2]
                                or channels[vchan][1] == '0' and key[0]
                                == channels[vchan][0]]:
                        [gwtrg_vetoed, gwtrg_postveto] = \
                            event.includeexclude(gwtrg_postveto, vetoseg,
                                tcent=event.ifo[method][ifo])
                else:
                    for ifo in [key for key in event.ifo[method].keys()
                                if key == (channels[vchan])[:2]
                                or channels[vchan][1] == '0' and key[0]
                                == channels[vchan][0]]:
                        gwtrg_postveto = event.exclude(gwtrg_postveto,
                                vetoseg, tcent=event.ifo[method][ifo])

                gwseg_postveto = event.removesegments(gwseg, vetoseg)

                vactual = ngwtrg - len(gwtrg_postveto)

                if nauxtrg > 0:
                    useP = 1. * vactual / nauxtrg
                else:
                    useP = 0.

                vsig = 0. - _gammpln(vactual, vexpected)
                veff = float(vactual) / float(ngwtrg)
                if deadfrac == 0:
                    effbydt = 0.
                else:
                    effbydt = veff / deadfrac

                r = [
                    livetime,
                    ngwtrg,
                    vchan,
                    vthr,
                    vwin,
                    deadsecs,
                    nauxtrg,
                    vexpected,
                    vactual,
                    vsig,
                    useP,
                    veff,
                    deadfrac,
                    effbydt,
                    ]
                if r[sD[__params.metric]] > best_metric:
                    best_metric = r[sD[__params.metric]]
                    best_lineid = lineid
                    best_gwtrg = gwtrg_postveto
                    bets_gwseg = gwseg_postveto

                    result = r
                    if track:
                        if len(gwtrg_vetoed) == 0:
                            gwtrg_vetoed = ['NONE']
                        track_result = [
                            livetime,
                            ngwtrg,
                            vchan,
                            vthr,
                            vwin,
                            deadsecs,
                            nauxtrg,
                            vactual,
                            vsig,
                            gwtrg_vetoed,
                            ]

            ### we now know which is the best
            if best_gwtrg == None:
                raise ValueError("could not find a new \"best\" veto config. Perhaps initial best_metric is too high?")

            stats.pop( best_lineid ) ### remove the line
            gwtrg = best_gwtrg ### remove the triggers
            gwseg = bets_gwseg ### remove the segments

            for line in stats:  # check if we need to remember the triggers for vchan
                if int(line[sD['vchan']]) == vchan:
                    break
            else:
                allvtrg[vchan] = 0

            ### report data
            # print results to main output file
            print >> f, \
                '%12.3f %7d %5d %4d %7.3f %11.3f %8d %6.2f %5d %8.2f %9.6f %9.5f %9.5f %9.4f %9.2f %9.2f' \
                % tuple(result + [time.clock() - t0, time.time() - t1] )

            if track:
                tfparray[lineidx + 1] = track_result
                lineidx += 1

            if verbose:
                print '%12.3f %7d %5d %4d %7.3f %11.3f %8d %6.2f %5d %8.2f %9.6f %9.5f %9.5f %9.4f %9.2f %9.2f' \
                    % tuple( result + [time.clock() - t0, time.time() - t1] )

            f.flush()

        ### we've now exhausted all veto configurations in the list
        f.close()
        statsfiles.append(filename)
        if track:
            pickle.dump(tfparray, tfp)
            tfp.close()
            trackfiles.append(trackfilename)

    return (statsfiles, trackfiles)

# =================================================

def vetolist_safety(
    __params,
    output_dir='./',
    source_dir='./',
    verbose=False,
    ):
    """ an encapsulation of vetolist_eval.py 
  __params is a params object
  if present, expect numBins to be a non-negative integer corresponding to the number of bins in a round-robin analysis.
    delegates job to mish_mash_eval() if numBins is present.
  """

    output_dir = __check_dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    source_dir = __check_dirname(source_dir)

    if vars(__params).has_key('channels'):
        channels = event.loadlist(__params.channels)
    else:
        channels = __load_channels(__params.auxdir,
                                   analysis_range=__params.analysis_range)

  # we need the same channels list as for recalculate.py so we can associate numbers with channels

    for remove in __params.notused:
        channels = [channel for channel in channels
                    if channel.find(remove) == -1]

    vetolists = []  # a list of the finished vetolist files

    for gwset in __params.gwsets:
        gwset += '-' + str(__params.gwthr)
        filename = output_dir + gwset + '.0' \
            + '.vetolist.safety'
        f = open(filename, 'w')
        print >> f, \
            '# %10s %7s %4s %-40s %4s %7s %11s %8s %6s %5s %8s %9s %9s %9s %9s %7s %12s %9s %11s' \
            % (
            'livetime',
            '#gwtrg',
            'bin',
            'vchan',
            'vthr',
            'vwin',
            'dsec',
            '#auxtrg',
            'vexp',
            'vact',
            'vsig',
            'useP',
            'eff',
            'dt',
            'eff/dt',
            'c_vact',
            'c_dsec',
            'c_auxtrg',
            '%s_exp' % __params.metric,
            )

        if verbose:
            print '  %10s %7s %4s %-40s %4s %7s %11s %8s %6s %5s %8s %9s %9s %9s %9s %7s %12s %9s %11s' \
                % (
                'livetime',
                '#gwtrg',
                'bin',
                'vchan',
                'vthr',
                'vwin',
                'dsec',
                '#auxtrg',
                'vexp',
                'vact',
                'vsig',
                'useP',
                'eff',
                'dt',
                'eff/dt',
                'c_vact',
                'c_dsec',
                'c_auxtrg',
                '%s_exp' % __params.metric,
                )

    # load in stats.0 data

        stats = event.load(source_dir + gwset + '.stats.0')

    # sort so they are in the correct order
        stats.sort(key=lambda line: line[sD['vthr']], reverse=True)  # this ordering should be consistent with standard ROC format (rank: low->high)
        stats.sort(key=lambda line: line[sD[__params.metric]], reverse=True)

        bin = 0  # put in place to match format for mish_mash_vetolists
        c_vact = 0.
        c_dsec = 0.
        c_auxtrg = 0.
        for line_idx in range(len(stats)):
            line = stats[line_idx]
            livetime = line[sD['livetime']]
            ngwtrg = line[sD['#gwtrg']]
            vchan = channels[int(line[sD['vchan']])]
            vthr = line[sD['vthr']]
            vwin = line[sD['vwin']]
            dsec = line[sD['dsec']]
            nauxtrg = line[sD['#auxtrg']]
            vexp = line[sD['vexp']]
            vact = line[sD['vact']]
            vsig = line[sD['vsig']]
            useP = line[sD['useP']]
            eff = line[sD['eff']]
            dt = line[sD['dt']]
            effbydt = line[sD['eff/dt']]

            c_vact += vact
            c_dsec += dsec
            c_auxtrg += nauxtrg
            metric_exp = stats[line_idx][sD[__params.metric]]

            print >> f, \
                '%12.3f %7d %4d %-40s %4d %7.3f %11.3f %8d %6.2f %5d %8.2f %9.6f %9.5f %9.5f %9.4f %7d %12.3f %9d %11.4f' \
                % (
                livetime,
                ngwtrg,
                bin,
                vchan,
                vthr,
                vwin,
                dsec,
                nauxtrg,
                vexp,
                vact,
                vsig,
                useP,
                eff,
                dt,
                effbydt,
                c_vact,
                c_dsec,
                c_auxtrg,
                metric_exp,
                )
            if verbose:
                print '%12.3f %7d %4d %-40s %4d %7.3f %11.3f %8d %6.2f %5d %8.2f %9.6f %9.5f %9.5f %9.4f %7d %12.3f %9d %11.4f' \
                    % (
                    livetime,
                    ngwtrg,
                    bin,
                    vchan,
                    vthr,
                    vwin,
                    dsec,
                    nauxtrg,
                    vexp,
                    vact,
                    vsig,
                    useP,
                    eff,
                    dt,
                    effbydt,
                    c_vact,
                    c_dsec,
                    c_auxtrg,
                    metric_exp,
                    )

        f.close()
        vetolists.append(filename)

    return vetolists

###
def vetolist_eval(
    run,
    __params,
    numBins=False,
    output_dir='./',
    source_dir='./',
    verbose=False,
    ):
    """ an encapsulation of vetolist_eval.py 
  run is the iteration number (label for which .stats file to use)
  __params is a params object
  if present, expect numBins to be a non-negative integer corresponding to the number of bins in a round-robin analysis.
    delegates job to mish_mash_eval() if numBins is present.
  """

    if numBins:  # delegate if using round-robin data
        return mish_mash_eval(
            run,
            numBins,
            __params,
            output_dir=output_dir,
            source_dir=source_dir,
            verbose=verbose,
            )

    output_dir = __check_dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    source_dir = __check_dirname(source_dir)

    if vars(__params).has_key('channels'):
        channels = event.loadlist(__params.channels)
    else:
        channels = __load_channels(__params.auxdir,
                                   analysis_range=__params.analysis_range)

  # we need the same channels list as for recalculate.py so we can associate numbers with channels

    for remove in __params.notused:
        channels = [channel for channel in channels
                    if channel.find(remove) == -1]

    vetolists = []  # a list of the finished vetolist files

    for gwset in __params.gwsets:
        gwset += '-' + str(__params.gwthr)
        filename = output_dir + gwset + '.' + repr(run) \
            + '.vetolist.eval'
        f = open(filename, 'w')
        print >> f, \
            '# %10s %7s %4s %-40s %4s %7s %11s %8s %6s %5s %8s %9s %9s %9s %9s %7s %12s %9s %11s' \
            % (
            'livetime',
            '#gwtrg',
            'bin',
            'vchan',
            'vthr',
            'vwin',
            'dsec',
            '#auxtrg',
            'vexp',
            'vact',
            'vsig',
            'useP',
            'eff',
            'dt',
            'eff/dt',
            'c_vact',
            'c_dsec',
            'c_auxtrg',
            '%s_exp' % __params.metric,
            )

        if verbose:
            print '# %10s %7s %4s %-40s %4s %7s %11s %8s %6s %5s %8s %9s %9s %9s %9s %7s %12s %9s %11s' \
                % (
                'livetime',
                '#gwtrg',
                'bin',
                'vchan',
                'vthr',
                'vwin',
                'dsec',
                '#auxtrg',
                'vexp',
                'vact',
                'vsig',
                'useP',
                'eff',
                'dt',
                'eff/dt',
                'c_vact',
                'c_dsec',
                'c_auxtrg',
                '%s_exp' % __params.metric,
                )

    # evaluation data

        stats = event.load(source_dir + gwset + '.stats.' + repr(run)
                           + '.eval')

    # training data prior to evaluation. yields expected performances

        rank_stats = event.load(source_dir + gwset + '.stats.'
                                + repr(run - 1))
        sigthr = -math.log(__params.Psigthr)
        effbydtthr = __params.effbydtthr

    # filter and re-order the rank_stats as done in recalculate()

        rank_stats = [line for line in rank_stats if line[sD['vsig']]
                      > sigthr and line[sD['eff/dt']] > effbydtthr]
        rank_stats.sort(key=lambda line: line[sD['vthr']], reverse=True)  # this ordering should be consistent with standard ROC format (rank: low->high)
        rank_stats.sort(key=lambda line: line[sD[__params.metric]],
                        reverse=True)

#    rank_stats.sort(key=lambda line: line[sD['eff/dt']], reverse=True)

        bin = 0  # put in place to match format for mish_mash_vetolists
        c_vact = 0.
        c_dsec = 0.
        c_auxtrg = 0.
        for line_idx in range(len(stats)):
            line = stats[line_idx]
            livetime = line[sD['livetime']]
            ngwtrg = line[sD['#gwtrg']]
            vchan = channels[int(line[sD['vchan']])]
            vthr = line[sD['vthr']]
            vwin = line[sD['vwin']]
            dsec = line[sD['dsec']]
            nauxtrg = line[sD['#auxtrg']]
            vexp = line[sD['vexp']]
            vact = line[sD['vact']]
            vsig = line[sD['vsig']]
            useP = line[sD['useP']]
            eff = line[sD['eff']]
            dt = line[sD['dt']]
            effbydt = line[sD['eff/dt']]

            c_vact += vact
            c_dsec += dsec
            c_auxtrg += nauxtrg
            metric_exp = rank_stats[line_idx][sD[__params.metric]]

            print >> f, \
                '%12.3f %7d %4d %-40s %4d %7.3f %11.3f %8d %6.2f %5d %8.2f %9.6f %9.5f %9.5f %9.4f %7d %12.3f %9d %11.4f' \
                % (
                livetime,
                ngwtrg,
                bin,
                vchan,
                vthr,
                vwin,
                dsec,
                nauxtrg,
                vexp,
                vact,
                vsig,
                useP,
                eff,
                dt,
                effbydt,
                c_vact,
                c_dsec,
                c_auxtrg,
                metric_exp,
                )
            if verbose:
                print '%12.3f %7d %4d %-40s %4d %7.3f %11.3f %8d %6.2f %5d %8.2f %9.6f %9.5f %9.5f %9.4f %7d %12.3f %9d %11.4f' \
                    % (
                    livetime,
                    ngwtrg,
                    bin,
                    vchan,
                    vthr,
                    vwin,
                    dsec,
                    nauxtrg,
                    vexp,
                    vact,
                    vsig,
                    useP,
                    eff,
                    dt,
                    effbydt,
                    c_vact,
                    c_dsec,
                    c_auxtrg,
                    metric_exp,
                    )

        f.close()
        vetolists.append(filename)

    return vetolists


# =================================================

def mish_mash_eval(
    run,
    numBins,
    __params,
    output_dir='./',
    source_dir='./',
    verbose=True,
    ):
    """ an encapsulation of mish-mash_eval.py, which creates combined vetolists for round-robin analysis
  run is the iteration number for the algorithm (labels which .stats.*.training* files to use)
  numBins is a non-negative integer descrbing the number of bins in the round-robin analysis
  __params is a params object
  """

    output_dir = __check_dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    source_dir = __check_dirname(source_dir)

    if vars(__params).has_key('channels'):
        channels = event.loadlist(__params.channels)
    else:
        channels = __load_channels(__params.auxdir,
                                   analysis_range=__params.analysis_range)

  # we need the same channels list as for recalculate.py so we can associate numbers with channels

    for remove in __params.notused:
        channels = [channel for channel in channels
                    if channel.find(remove) == -1]

    sigthr = -math.log(__params.Psigthr)
    effbydtthr = __params.effbydtthr

    vetolists = []  # a list of the finished vetolist files

    for gwset in __params.gwsets:
        gwset += '-' + str(__params.gwthr)

    # create our cumulative data storage items

        filename = output_dir + gwset + '.' + repr(run) \
            + '.eval.mish-mash.vetolist'
        f = open(filename, 'w')
        Stats = []
        aStats = []
        c_livetime = 0.
        c_numgwtrg = 0.

    # print headers on f

        print >> f, \
            '# %10s %7s %4s %-40s %4s %7s %11s %8s %6s %5s %8s %9s %9s %9s %9s %7s %12s %9s %11s' \
            % (
            'livetime',
            '#gwtrg',
            'bin',
            'vchan',
            'vthr',
            'vwin',
            'dsec',
            '#auxtrg',
            'vexp',
            'vact',
            'vsig',
            'useP',
            'eff',
            'dt',
            'eff/dt',
            'c_vact',
            'c_dsec',
            'c_auxtrg',
            '%s_exp' % __params.metric,
            )

        if verbose:
            print '# %10s %7s %4s %-40s %4s %7s %11s %8s %6s %5s %8s %9s %9s %9s %9s %7s %12s %9s %11s' \
                % (
                'livetime',
                '#gwtrg',
                'bin',
                'vchan',
                'vthr',
                'vwin',
                'dsec',
                '#auxtrg',
                'vexp',
                'vact',
                'vsig',
                'useP',
                'eff',
                'dt',
                'eff/dt',
                'c_vact',
                'c_dsec',
                'c_auxtrg',
                '%s_exp' % __params.metric,
                )

    # create a dictionary to store the analysis data

        lookup = dict()

    # load the individual stats files from all bins

        for bin in range(numBins):

      # We use the final training data to order our list and then insert the analysis data while preserving that order.

            stats = event.load(source_dir + gwset + '.stats.'
                               + repr(run - 1) + '.training.bin_'
                               + repr(bin) + '_out_of_' + repr(numBins))  # training data
            astats = event.load(source_dir + gwset + '.stats.'
                                + repr(run) + '.eval.analysis.bin_'
                                + repr(bin) + '_out_of_'
                                + repr(numBins))  # analysis data

      # pull out global parameters from each bin. This is only valid for the analysis data. Training data overlaps

            c_livetime += float(astats[0][0])
            c_numgwtrg += int(astats[0][1])

      # append each line in stats with the bin number as an identifier

            for line_idx in range(len(stats)):
                stats[line_idx] += [bin]

      # fill in lookup using analysis data. sorted by (vchan-vthr-vwin-bin)

            for line in astats:
                lookup['%d-%d-%.3f-%d' % (line[sD['vchan']],
                       line[sD['vthr']], line[sD['vwin']], bin)] = line \
                    + [bin]

      # add training data to the cumulative list

            Stats += stats

    # sort training data (Stats) by the effbydt
#    Stats.sort(key = lambda line: float(line[sD['eff/dt']]), reverse = True)

        Stats.sort(key=lambda line: float(line[sD[__params.metric]]),
                   reverse=True)

    # fill in analysis data (aStats) using the data stored in lookup while preserving the order of Stats

        aStats = []
        for line in Stats:

      # check to make sure the configuration will appear in the analysis data. (to clarify, Stats stores training data from 'run-1' and the config will only be in the analysis data for 'run' if it passed these performance thresholds)

            if line[sD['eff/dt']] > effbydtthr and line[sD['vsig']] \
                > sigthr:
                aStats += [lookup['%d-%d-%.3f-%d' % (line[sD['vchan']],
                           line[sD['vthr']], line[sD['vwin']],
                           line[-1])] + [line[sD[__params.metric]]]]
            else:
                pass

    # compute the cumulative statistics using aStats. We convert to efficiency and fractional deadtime in the output

        c_dsec = 0.
        c_vact = 0.
        c_auxtrg = 0.
        livetime = c_livetime
        numgwtrg = c_numgwtrg
        for line in aStats:
            livetime = line[sD['livetime']]
            ngwtrg = line[sD['#gwtrg']]
            vchan = channels[int(line[sD['vchan']])]
            vthr = line[sD['vthr']]
            vwin = line[sD['vwin']]
            dsec = line[sD['dsec']]
            nauxtrg = line[sD['#auxtrg']]
            vexp = line[sD['vexp']]
            vact = line[sD['vact']]
            vsig = line[sD['vsig']]
            eff = line[sD['eff']]
            dt = line[sD['dt']]
            effbydt = line[sD['eff/dt']]
            bin = line[-2]
            metric_exp = line[-1]

      # compute cumulative statistics

            c_vact += vact
            c_dsec += dsec
            c_auxtrg += nauxtrg

      # compute the updated ngwtrg and livetime for each line

            livetime -= dsec
            numgwtrg -= int(vact)

      # print to file in a readable way:

            print >> f, \
                '%12.3f %7d %4d %-40s %4d %7.3f %11.3f %8d %6.2f %5d %8.2f %9.6f %9.5f %9.5f %9.4f %7d %12.3f %9d %11.4f' \
                % (
                livetime,
                ngwtrg,
                bin,
                vchan,
                vthr,
                vwin,
                dsec,
                nauxtrg,
                vexp,
                vact,
                vsig,
                useP,
                eff,
                dt,
                effbydt,
                c_vact,
                c_dsec,
                c_auxtrg,
                metric_exp,
                )
            if verbose:
                print '%12.3f %7d %4d %-40s %4d %7.3f %11.3f %8d %6.2f %5d %8.2f %9.6f %9.5f %9.5f %9.4f %7d %12.3f %9d %11.4f' \
                    % (
                    livetime,
                    ngwtrg,
                    bin,
                    vchan,
                    vthr,
                    vwin,
                    dsec,
                    nauxtrg,
                    vexp,
                    vact,
                    vsig,
                    useP,
                    eff,
                    dt,
                    effbydt,
                    c_vact,
                    c_dsec,
                    c_auxtrg,
                    metric_exp,
                    )

        f.close()

        vetolists.append(filename)

    return vetolists


# ===================================================================================================
#
#
#                          ovl output handling functions
#
#
# ===================================================================================================

def stats_to_ROC(
    statsfilename,
    __params=False,
    rank_statsfilename=False,
    cleans_statsfilename=False,
    output_dir='./',
    source_dir='./',
    cleans_source_dir='./',
    ):
    """ generates OVL output in the `standard ROC format' from a stats file 
  expect statsfilename to be the output from recalculate(eval=True)
  expect rank_statsfilename to be the output from recalculate(eval=False) that directly preceeded statsfilename (and from which we can derive rankings)
  expect __params to be a params object. Only required if rank_statsfilename is present
  expect clans_statsfilename to be the output from recalcuate(eval=True) using cleans instead of glitches. implemented for compatibility with AuxMVC pipelines
  """

    output_dir = __check_dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    source_dir = __check_dirname(source_dir)
    cleans_source_dir = __check_dirname(cleans_source_dir)

    if 'eval' not in statsfilename:
        go_on = \
            raw_input("you are not using an 'eval' file. do you wish to continue anyway? (yes/no)"
                      )
        if go_on != 'yes' and go_on != 'y':
            return False

    stats = event.load(source_dir + statsfilename)

    if rank_statsfilename and __params:
        sigthr = -math.log(__params.Psigthr)
        effbydtthr = __params.effbydtthr
        rank_stats = event.load(source_dir + rank_statsfilename)

    # filter and re-order the rank_stats as done in recalculate()

        rank_stats = [line for line in rank_stats if line[sD['vsig']]
                      > sigthr and line[sD['eff/dt']] > effbydtthr]
        rank_stats.sort(key=lambda line: line[sD['vthr']], reverse=True)  # this ordering should be consistent with standard ROC format (rank: low->high)
        rank_stats.sort(key=lambda line: line[sD[__params.metric]],
                        reverse=True)

#    rank_stats.sort(key=lambda line: line[sD['eff/dt']], reverse=True)

    if cleans_statsfilename:
        cleans_stats = event.load(cleans_source_dir
                                  + cleans_statsfilename)

    roc_filename = (output_dir + statsfilename + '.roc', 'w')
    f = open(roc_filename, 'w')

    print >> f, '# line 1 : total number of gw triggers'
    if cleans_statsfilename:
        print >> f, '# line 2 : total number of clean samples'
        print >> f, \
            '# line 3 : rank(%s) cumulative_number_of_cleans_removed cumulative_number_of_GWtrg_removed'%__params.metric
    else:
        print >> f, '# line 2 : total amount of livetime'
        print >> f, \
            '# line 3 : rank(%s) cumulative_deadtime cumulative_number_of_GWtrg_removed'%__params.metric

    print >> f, stats[0][sD['#gwtrg']]  # total number of gwtrg
    if cleans_statsfilename:
        print >> f, cleans_stats[0][sD['#gwtrg']]
    else:
        print >> f, stats[0][sD['livetime']]  # total amount of livetime

    c_vact = 0
    c_dsec = 0

    for index in range(len(stats)):
        line = stats[index]
        c_vact += line[sD['vact']]  # the number of glitches removed

        if cleans_statsfilename:
            c_dsec += cleans_stats[index][sD['vact']]  # the number of cleans removed
        else:
            c_dsec += line[sD['dsec']]  # the number of seconds removed

        if rank_statsfilename and __params:
            if __params.metric == "eff/dt":
                rank = effbydt_to_rank(rank_stats[index][sD['eff/dt']])
            elif __params.metric == "vsig":
                rank = vsig_to_rank(rank_stats[index][sD['vsig']]) 
            elif __params.metric == "useP":
                rank = useP_to_rank(rank_stats[index][sD['useP']])
        else:
            rank = 'na'

        print >> f, rank, c_dsec, c_vact

    f.close()

    return roc_filename


# =================================================

def vetolist_to_ROC(
    vetolist_filename,
    cleans_vetolist_filename=False,
    output_dir='./',
    source_dir='./',
    cleans_source_dir='./',
    metric='eff/dt',
    ):
    """ 
  written only for non-round robin analysis. need to expand the funcionality

  converts a vetolist into the standard ROC format 
  """

    output_dir = __check_dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    source_dir = __check_dirname(source_dir)
    cleans_source_dir = __check_dirname(cleans_source_dir)

    vetolist = [line.strip().split() for line in open(source_dir
                + vetolist_filename, 'r') if line[0] != '#']

    if cleans_vetolist_filename:
        cleans_vetolist = [line.strip().split() for line in
                           open(cleans_source_dir
                           + cleans_vetolist_filename, 'r') if line[0]
                           != '#']

    roc_filename = output_dir + vetolist_filename + '.roc'
    f = open(roc_filename, 'w')
    print >> f, '# line 1 : total number of gw triggers'
    if cleans_vetolist_filename:
        print >> f, '# line 2 : total number of clean samples'
        print >> f, \
            '# line 3 : rank(eff/dt) cumulative_number_of_cleans_removed cumulative_number_of_GWtrg_removed'
    else:
        print >> f, '# line 2 : total amount of livetime'
        print >> f, \
            '# line 3 : rank(%s) cumulative_deadtime cumulative_number_of_GWtrg_removed' \
            % metric

    if len(vetolist) == 0:  # empty list
        print >> f, '''0
0
0 0 0'''
    else:
        print >> f, vetolist[0][vD['#gwtrg']]
        if cleans_vetolist_filename:
            print >> f, cleans_vetolist[0][vD['#gwtrg']]
        else:
            print >> f, vetolist[0][vD['livetime']]
        for index in range(len(vetolist)):
            c_vact = vetolist[index][vD['c_vact']]
            if cleans_vetolist_filename:
                c_dsec = cleans_vetolist[index][vD['c_vact']]
            else:
                c_dsec = vetolist[index][vD['c_dsec']]
            if metric == 'eff/dt':
                rank = \
                    effbydt_to_rank(float(vetolist[index][vD['metric_exp'
                                    ]]))
            elif metric == 'vsig':
                rank = \
                    vsig_to_rank(float(vetolist[index][vD['metric_exp'
                                 ]]))
            elif metric == 'useP':
                rank = \
                    useP_to_rank(float(vetolist[index][vD['metric_exp'
                                 ]]))
            else:
                raise ValueError('unknown metric=', metric)
            print >> f, rank, c_dsec, c_vact

    f.close()

    return roc_filename


# =================================================

def vetolist_to_segments(
    vetolist_filename,
    gpsstart,
    gpsstop,
    fractional_dt,
    trg_dir='./',
    scisegs=False,
    output_filename=False,
    ):
    """ 
  expects file structure to be: trg_dir / basename-(gpsstart/1e5) / *trg
    these will be multi-channel files

  need some sort of scisegs or something to define when is "good" time
  """

  # load vetolist and figure out which channels are present

    v = []
    allvtrg = {}
    for line in open(vetolist_filename, 'r'):
        if line[0] != '#':
            line = line.strip().split()
            v.append(line)
            if line[vD['vchan']] not in allvtrg.keys():
                allvtrg[line[vD['vchan']]] = []

    vchans = allvtrg.keys()

  # fill in allvtrg using triggers in trg_dir

    for dir in sorted(os.listdir(trg_dir)):
        try:
            dir_start = int(dir.split('-')[-1])  # get gpsstart time of directory
            if dir_start <= gpsstart / 1e5 <= dir_start + 1 or dir_stop \
                <= gpsstop / 1e5 <= dir_start + 1:  # directory covers part of [gpsstart, gpsstop]

        # load triggers from .trg files (multichannel files)

                for f in sorted(os.listdir(trg_dir + '/' + dir)):
                    if f[-4:] != '.trg':
                        continue
                    name = f.strip('.trg').split('-')
                    if gpsstart <= int(name[-2]) <= gpsstop or gpsstart \
                        <= int(name[-2]) + int(name[-1]) <= gpsstop:
                        for line in open(trg_dir + '/' + dir + '/' + f,
                                'r'):
                            if line[0] != '#':
                                line = line.strip().split()
                                channel = line[kw_trgD['chan']]
                                if channel in vchans:
                                    allvtrg[channel].append([float(line[kw_trgD['GPScent'
        ]]), float(line[kw_trgD['signif']])])
        except:
            pass

  # filter with scisegs if present

    if scisegs:
        segs = event.load(scisegs)  # load in scisegs
        segs = event.andsegments([[[gpsstart, gpsstop]], segs])
    else:
        segs = [[gpsstart, gpsstop]]

  # total livetime present

    lvtm = event.livetime(segs)

  # current dt

    dt = 0.

  # build vetosegs and measure time removed

    v_ind = 0
    veto_segs = []
    while v_ind < len(v):

    # pull out configuration parameters

        line = v[v_ind]
        vchan = line[vD['vchan']]
        vwin = float(line[vD['vwin']])
        vthr = float(line[vD['vthr']])

    # generate segments

        vtrg = allvtrg[vchan]
        vtrg.sort(key=lambda line: line[0])
        v_seg = event.vetosegs(vtrg, vwin, threshold=vthr, tcent=0,
                               signif=1)
        v_seg = event.andsegments([v_seg, segs])  # take only the overlap with segs

        v_dt = event.livetime(v_seg)  # additional deadtime

    # if we haven't exceeded the requested amount of deadtime or this step puts us closer to the requested amount of deadtime, we include these segments

        if v_dt + dt < lvtm * fractional_dt:  # or ( abs(v_dt+dt - lvtm*fractional_dt) < abs(dt - lvtm*fractional_dt) ):
            segs = event.removesegments(segs, v_seg)  # remove v_seg from segs
            dt += v_dt
            veto_segs = event.orsegments([v_seg, veto_segs])  # add v_seg to total list of veto_segs
        else:
            break

        v_ind += 1

  # write segments to an output file

    if output_filename:
        f = open(output_filename, 'w')
        for s in veto_segs:
            print >> f, s[0], s[1]
        f.close()

    return (veto_segs, float(dt) / lvtm)


# ===================================================================================================
#
#
#                              functions for use with iDQ
#
#
# ===================================================================================================

def patfile_to_GPStimes(auxmvc_pat, skip_lines=1):
    """ reads in auxmvc *.pat files and extracts lists of GPStimes suitable for predict()
  expect auxmvc_pat to be the full path to the *.pat file
  skip_lines is a non-negative integer corresponding to the number of lines we skip when reading the *.pat file

  for multiple *.pat files, call this method repeatedly and concatinate the returned lists
  """

    pat_lines = open(auxmvc_pat).readlines()
    pat_lines = pat_lines[skip_lines:]

    variables = dict([(line, i) for (i, line) in enumerate(pat_lines[0].split())])

    GPStimes = []
    if variables.has_key('i'):
        for line in pat_lines[1:]:
            line = line.strip().split()
            GPStimes.append([float(line[variables['GPS_s']])
                            + float(line[variables['GPS_ms']]) * 1e-3,
                            int(int(line[variables['i']]))])

    elif variables.has_key('unclean'):
        for line in pat_lines[1:]:
            line = line.strip().split()
            GPStimes.append([float(line[variables['GPS_s']])
                            + float(line[variables['GPS_ms']]) * 1e-3,
                            int(float(line[variables['unclean']]))])

    else:
        for line in pat_lines[1:]:
            line = line.strip().split()
            GPStimes.append([float(line[variables['GPS_s']])
                            + float(line[variables['GPS_ms']]) * 1e-3,
                            int(float(line[-1]))]) # class stored in the last column. This may be buggy...

    return GPStimes


# =================================================

def kw_trigfile_to_GPStimes(kw_trgfile, gwchan):
    """ reads in KW_trigfiles and generates lists of GPS times suitable for predict()
  expect kw_trgfile to be the full path to the *.trg file
    expects the kw_trgfile to be the direct output of the KW pipeline, with a name similar to : H-KW_TRIGGERS-1028247424-32.trg
  expect gwchan to be the name of the GW channel

  for multiple *.pat files, call this method repeatedly and concatinate the returned lists
  """

    GPStimes = []
    for line in open(kw_trgfile):
        if line[0] != '#':
            line = line.strip().split()
            if line[kw_trgD['chan']] == gwchan:
                GPStimes.append([float(line[kw_trgD['GPScent']]), 1])

    return GPStimes


# =================================================

def predict(
    vetolist,
    headers,
    GPStimes,
    gps_tcent,
    allvtrg=False,
    kw_trgfiles=False,
    predict_filename=False,
    output_dir='./',
    metric='eff/dt',
    ):
    """ makes predictions about GPS times given a vetolist and KW-TRIGGERS.trg file
  expect vetolist to be an ASCII filename corresponding to the output of vetolist_eval()
  expect kw_trgfile to be an output file directly from the KW pipeline withe a name similar to : H-KW_TRIGGERS-1028247424-32.trg
  expect GPStimes to be a list of tuples (iterables) with : tuple = (GPS-tcent, 0/1)
    0 : clean sample (for compatibility with auxmvc pipelines)
    1 : glitch sample
  expect predict_filename to be a string corresponding to a filename in which we wish to write the prediction results
  expect ouput_dir to be a string corresponding to a directory into which we will write "predict_filename"
  """

    output_dir = __check_dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not allvtrg:

    # load vetolist (as strings) and auxiliary channel names

        v = []  # vetollist
        allvtrg = {}  # dictionary linking channel names to triggers
        for line in open(vetolist, 'r'):
            if line[0] != '#':
                line = line.strip().split()
                v.append(line)
                vchan = line[vD['vchan']]
                if vchan not in allvtrg.keys():
                    allvtrg[vchan] = []

    # load auxilary triggers

        for kw_trgfile in kw_trgfiles:
            for line in open(kw_trgfile, 'r'):
                if line[0] != '#':
                    line = line.strip().split()
                    chan = line[kw_trgD['chan']]
                    if chan in allvtrg.keys():
                        allvtrg[chan].append(line[:-1])
    else:
        v = []
        for line in open(vetolist, 'r'):
            if line[0] != '#':
                line = line.strip().split()
                v.append(line)
                vchan = line[vD['vchan']]
                if vchan not in allvtrg.keys():
                    allvtrg[vchan] = []

  # associate glitches with vconfigs

    gwtrg_predict = []
    c_dsec_act = 0.
    for line in v:

    # pull out configuration parameters

        vchan = line[vD['vchan']]
        vthr = float(line[vD['vthr']])
        vwin = float(line[vD['vwin']])
        vtrg = [(float(l[kw_trgD['GPScent']]), float(l[kw_trgD['signif'
                ]])) for l in allvtrg[vchan]]

        vtrg.sort(key=lambda line: line[0])

    # generate veto segements

        vetoseg = event.vetosegs(vtrg, vwin, vthr, tcent=0, signif=1)

    # winnow GPS times

        [gwtrg_vetoed, GPStimes] = event.includeexclude(GPStimes,
                vetoseg, tcent=gps_tcent)  # tcent comes from the required format of GPStimes

        if metric == 'eff/dt':
            metric_exp = float(line[vD['metric_exp']])
            rank = effbydt_to_rank(metric_exp)
        elif metric == 'vsig':
            metric_exp = float(line[vD['metric_exp']])
            rank = vsig_to_rank(metric_exp)
        elif metric == 'useP':
            metric_exp = float(line[vD['metric_exp']])
            rank = useP_to_rank(metric_exp)
        else:
            raise ValueError('unkown metric=%s' % metric)

    # associate vconfig with each GWtrg vetoed

        for gwtrg in gwtrg_vetoed:

      # need GPStime, glitch/clean, rank, vconfig info (take rank from vetolist!)

            gwtrg_predict.append(gwtrg + [rank, metric_exp, vchan,
                                 vthr, vwin])

  # for remaining GPS times, append gwtrg_predict with rank = 0

    for gwtrg in GPStimes:
        gwtrg_predict.append(gwtrg + [0., 0., 'none', 0., 0.])

  # sort by gps times

    gwtrg_predict.sort(key=lambda line: line[gps_tcent])

    if predict_filename:  # write output file
        f = open(output_dir + predict_filename, 'w')
        s = ''
        for S in headers:
            s += S + ' '
        print >> f, s + 'rank %s vchan vthr vwin' % metric
        for line in gwtrg_predict:
            s = ''
            for l in line:
                if not isinstance(l, str):
                    l = repr(l)
                s += l + ' '
            print >> f, s

        f.close()

    return (gwtrg_predict, predict_filename)


##################################################
# convert veto list to timeseries based on KW triggers provided
# note this has a lot of overlap with ovl.predict(), but most code here will have duplicate functionality for clarity
# may be good to merge these two functions in the future
# vetolist: filename of vetolist file
# allvtrg: kw trigdict of triggers (dictionary of kw tables, key indexed by channel)
# segment: segment over which to construct timeseries
# fs: sampling rate, default = 128 Hz

def vetolist2timeseries(
    vetolist,
    allvtrg,
    segment,
    fs=256,
    metric='eff/dt',
    ):

  # current vetolist format is not a standard structure, idq_realtime reads as
  # vetolist = open(ovl_train_cache, "r").readlines()[-1].strip('\n')

    import numpy as np

  # from ovl.predict(), we have:

    v = []
    for line in open(vetolist, 'r'):
        if line[0] != '#':
            fields = line.strip().split()
            fieldsbykey = dict((key, fields[vD[key]]) for key in
                               vD.keys())
            v.append(fieldsbykey)
            if fieldsbykey['vchan'] not in allvtrg.keys():
                allvtrg[fieldsbykey['vchan']] = []

    livetime = segment[1] - segment[0]
    ts = np.zeros(livetime * fs)
    dt = 1. / fs

    for line in reversed(v):
        triggers = allvtrg[line['vchan']]

    # print line['vchan'], len(triggers), line['vthr'], len([a for a in triggers if a[-1] > float(line['vthr'])])
    # veto segments within our interval

        vetoseg = event.vetosegs(triggers, float(line['vwin']),
                                 float(line['vthr']))
        if metric == 'eff/dt':
            vetorank = effbydt_to_rank(float(line['metric_exp']))
        elif metric == 'vsig':
            vetorank = vsig_to_rank(float(line['metric_exp']))
        elif metric == 'useP':
            vetorank = useP_to_rank(float(line['metric_exp']))
        else:
            raise ValueError('unknown metric=', metric)
        for seg in vetoseg:

      # this calculation assumes timeseries defines start times of veto segments of length dt
      # nearest neighbor to digitize (which is why rounding is done)
      # to be inclusive, round can be replaced by math.floor, math.ceil respectively
      # istart = max(0, int(round((seg[0] - segment[0]) * fs))) # start index
      # istop  = min(len(dt), int(round((seg[1] - segment[0]) * fs))) # 1 + stop index

            istart = max(0, int(math.floor((seg[0] - segment[0]) * fs)))  # start index
            istop = min(len(ts), int(math.ceil((seg[1] - segment[0])
                        * fs)))  # 1 + stop index
            ts[istart:istop] = vetorank

    return ts


# =================================================

def train(
    __params,
    num_runs=9,
    incremental=1000,
    output_dir='./',
    verbose=False,
    write_channels=False,
    optimal=True,
    ):
    """ trains OVL on a specified data set 
  num_runs is a non-negative integer corresponding to the maximum number of iterations for the algorithm
  incremental is the number of vconfigs that are independently applied
  __params is a params object
  output_dir is the directory into which all summary files will be written
  verbose and write_channels are passed to recalculate() as is
  """

  # set up output directory

    output_dir = __check_dirname(output_dir)
    OVL_dir = output_dir + 'ovl/'
    if not os.path.exists(OVL_dir):
        os.makedirs(OVL_dir)

  # write params file into the directory

    write_params(__params, 'params.txt', params_dir=OVL_dir)

  # train OVL

    for run in range(num_runs):
        recalculate(
            run,
            incremental,
            __params,
            source_dir=OVL_dir,
            output_dir=OVL_dir,
            eval=False,
            verbose=verbose,
            )

    if optimal: ### find the best order for the remaining channels. More expensive than just another recalculate, but guaranteed to produce an optimal result
        optimize(
            num_runs,
            __params,
            source_dir=OVL_dir,
            output_dir=OVL_dir,
            track=True,
            verbose=verbose,
            write_channels=write_channels,
            )
    else: ### just another recalculate
        recalculate(  # we only write the channel list on the last iteration
            num_runs,
            incremental,
            __params,
            source_dir=OVL_dir,
            output_dir=OVL_dir,
            track=True,
            eval=True,
            verbose=verbose,
            write_channels=write_channels,
            )

  # generate vetolist
    vetolists = vetolist_eval(num_runs, __params, output_dir=OVL_dir,
                              source_dir=OVL_dir, verbose=verbose)

  # generate standard ROC file
    for v in vetolists:
        v = v.split('/')[-1]
        roc_filename = vetolist_to_ROC(v, source_dir=OVL_dir,
                output_dir=OVL_dir, metric=__params.metric)

    return vetolists

###
def convergent_train(__params, output_dir="./", verbose=False, write_channels=False):
    """
    trains OVL on a specified data set
    increments incremental to keep pace with num_runs and auto-terminates iteration only when incremental==num_configs_in_list
    once we run out of configs to keep in the list, we iterate "num_runs" more times to shake out any loose ends...
        ==> We're still susceptible to Markovian cycles!
    """

    # set up output directory
    output_dir = __check_dirname(output_dir)
    OVL_dir = output_dir + 'ovl/'
    if not os.path.exists(OVL_dir):
        os.makedirs(OVL_dir)

    # write params file into the directory
    write_params(__params, 'params.txt', params_dir=OVL_dir)

    # train OVL
    ### initial run (everything is independent)
    run = 0
    while True:
        statsfiles = recalculate( run, run, __params, source_dir=OVL_dir, output_dir=OVL_dir, eval=False, verbose=verbose, write_channels=False)[0]
        run += 1

        ### see how many surviving configurations there are
        for statsfile in statsfiles:
            nconfigs = len([line for line in open("%s/%s"%(OVL_dir, statsfile), "r").readlines() if line[0] != "#"])
            if nconfigs > run: ### there are configs that are not applied hierarchically
                break
        else: ### only reached if all configs were applied hierarchically for all statsfiles
            break

    incremental = run-1

    # we only write the channel list on the last iteration
    ### find the best possible order for these configurations
    ### more expensive than a simple recalculate run, but guaranteed to reach an optimal result
    optimize( run, __params, source_dir=OVL_dir, output_dir=OVL_dir, track=True, verbose=verbose, write_channels=write_channels)

    # generate vetolist
    vetolists = vetolist_eval(run, __params, output_dir=OVL_dir, source_dir=OVL_dir, verbose=verbose)

    # generate standard ROC file
    for v in vetolists:
        v = v.split('/')[-1]
        roc_filename = vetolist_to_ROC(v, source_dir=OVL_dir,
                output_dir=OVL_dir, metric=__params.metric)

    return vetolists

###
def safety(__params, output_dir="./", verbose=False, write_channels=False):
    """
    perform a safety study, applying all configurations independently and then ranking them by correlation.
    outputs in vetolist format
    """

    # set up output directory
    output_dir = __check_dirname(output_dir)
    OVL_dir = output_dir + "ovl/"
    if not os.path.exists(OVL_dir):
        os.makedirs(OVL_dir)

    # write params file into directory
    write_params(__params, 'params.txt', params_dir=OVL_dir)

    # perform one recalculate run
    statsfiles = recalculate( 0, 0, __params, source_dir=OVL_dir, output_dir=OVL_dir, eval=False, verbose=verbose, write_channels=write_channels, track=True)[0]

    # sort channels and write into vetolist format
    vetolists = vetolist_safety(__params, output_dir=OVL_dir, source_dir=OVL_dir, verbose=verbose)

    # generate standard ROC file
    for v in vetolists:
        v = v.split('/')[-1]
        roc_filename = vetolist_to_ROC(v, source_dir=OVL_dir,
                output_dir=OVL_dir, metric=__params.metric)

    return vetolists

# ===================================================================================================
#
#
#                       special functions for computing statistics
#      Numerical Recipies adaptations by Tom Loredo (loredo@spacenet.tn.cornell.edu)
#       _gserln, _gcfln, _gammpln, _gammqln added by Lindy Blackburn (lindy@mit.edu)
# poisson significance can be calculated: -_gammpln(measured number, expected) = -log(1 - poisscdf(measured - 1, expected) in matlab
# chisq significance can be calculated: -_gammqln(degrees-of-freedom/2, measurement/2) = -log(1 - chi2cdf(measurement, dof)) in matlab
# matlab examples can be calculated more accurtely with gamminc, though the log functions are necessary for very small probabilities
# ===================================================================================================

# Exceptions:

max_iters = 'Too many iterations: '

# ============= Globals ===============

_rt2 = math.sqrt(2.)
_gammln_cof = [
    76.18009173,
    -86.50532033,
    24.01409822,
    -1.231739516e0,
    0.120858003e-2,
    -0.536382e-5,
    ]
_gammln_stp = 2.50662827465


# ============= Gamma, Incomplete Gamma ===========

def _gammln(xx):
    """Logarithm of the gamma function."""

    global _gammln_cof, _gammln_stp
    x = xx - 1.
    tmp = x + 5.5
    tmp = (x + 0.5) * math.log(tmp) - tmp
    ser = 1.
    for j in range(6):
        x = x + 1.
        ser = ser + _gammln_cof[j] / x
    return tmp + math.log(_gammln_stp * ser)


# ========================....

def _gser(
    a,
    x,
    itmax=10000,
    eps=3.e-7,
    ):
    """Series approx'n to the incomplete gamma function."""

    gln = _gammln(a)
    if x < 0.:
        raise bad_arg, x
    if x == 0.:
        return 0.
    ap = a
    sum = 1. / a
    delta = sum
    n = 1
    while n <= itmax:
        ap = ap + 1.
        delta = delta * x / ap
        sum = sum + delta
        if abs(delta) < abs(sum) * eps:
            return (sum * math.exp(-x + a * math.log(x) - gln), gln)
        n = n + 1
    raise StandardError, "%s_%s"%(max_iters, str((abs(delta), abs(sum) * eps)))


# ========================

def _gserln(
    a,
    x,
    itmax=10000,
    eps=3.e-7,
    ):
    """Series approx'n to the incomplete gamma function."""

    gln = _gammln(a)
    if x < 0.:
        raise bad_arg, x
    if x == 0.:
        return 0.
    ap = a
    sum = 1. / a
    delta = sum
    n = 1
    while n <= itmax:
        ap = ap + 1.
        delta = delta * x / ap
        sum = sum + delta
        if abs(delta) < abs(sum) * eps:
            return (math.log(sum) + -x + a * math.log(x) - gln, gln)
        n = n + 1
    raise StandardError, "%s_%s"%(max_iters, str((abs(delta), abs(sum) * eps)))


# ========================

def _gcf(
    a,
    x,
    itmax=1000,
    eps=3.e-7,
    ):
    """Continued fraction approx'n of the incomplete gamma function."""

    gln = _gammln(a)
    gold = 0.
    a0 = 1.
    a1 = x
    b0 = 0.
    b1 = 1.
    fac = 1.
    n = 1
    while n <= itmax:
        an = n
        ana = an - a
        a0 = (a1 + a0 * ana) * fac
        b0 = (b1 + b0 * ana) * fac
        anf = an * fac
        a1 = x * a0 + anf * a1
        b1 = x * b0 + anf * b1
        if a1 != 0.:
            fac = 1. / a1
            g = b1 * fac
            if abs((g - gold) / g) < eps:
                return (g * math.exp(-x + a * math.log(x) - gln), gln)
            gold = g
        n = n + 1
    raise StandardError, "%s_%s"%(max_iters, str((abs(delta), abs(sum) * eps)))


# ========================

def _gcfln(
    a,
    x,
    itmax=1000,
    eps=3.e-7,
    ):
    """Continued fraction approx'n of the incomplete gamma function."""

    gln = _gammln(a)
    gold = 0.
    a0 = 1.
    a1 = x
    b0 = 0.
    b1 = 1.
    fac = 1.
    n = 1
    while n <= itmax:
        an = n
        ana = an - a
        a0 = (a1 + a0 * ana) * fac
        b0 = (b1 + b0 * ana) * fac
        anf = an * fac
        a1 = x * a0 + anf * a1
        b1 = x * b0 + anf * b1
        if a1 != 0.:
            fac = 1. / a1
            g = b1 * fac
            if abs((g - gold) / g) < eps:
                return (math.log(g) + -x + a * math.log(x) - gln, gln)
            gold = g
        n = n + 1
    raise StandardError, "%s_%s"%(max_iters, str((abs(delta), abs(sum) * eps)))


# ========================

def _gammp(a, x):
    """lower Incomplete gamma function."""

    if a == 0:
        return 1.
    if x < 0. or a <= 0.:
        raise ValueError, (a, x)
    if x < a + 1.:
        return _gser(a, x)[0]
    else:
        return 1. - _gcf(a, x)[0]


# ========================

def _gammpln(a, x):
    """lower Incomplete gamma function."""

    if a == 0:
        return math.log(1.)
    if x < 0. or a <= 0.:
        raise ValueError, (a, x)
    if x < a + 1.:
        return _gserln(a, x)[0]
    else:
        return math.log(1. - _gcf(a, x)[0])


# ========================

def _gammq(a, x):
    """upper Incomplete gamma function."""

    if a == 0:
        return 0.
    if x < 0. or a <= 0.:
        raise ValueError, repr((a, x))
    if x < a + 1.:
        return 1. - _gser(a, x)[0]
    else:
        return _gcf(a, x)[0]


# ========================

def _gammqln(a, x):
    """upper Incomplete gamma function."""

    if a == 0:
        return math.log(0.)
    if x < 0. or a <= 0.:
        raise ValueError, repr((a, x))
    if x < a + 1.:
        return math.log(1. - _gser(a, x)[0])
    else:
        return _gcfln(a, x)[0]


# ======== Error function, normal CDF and inverse ================

def _ncdf_inv(p):
    """Inverse of the normal CDF."""

    c0 = 2.515517
    c1 = 0.802853
    c2 = 0.010328
    d1 = 1.432788
    d2 = 0.189269
    d3 = 0.001308

    sign = -1.
    if p > 0.5:
        sign = 1.
        p = 1. - p
    arg = -2. * math.log(p)
    t = math.sqrt(arg)
    g = t - (c0 + t * (c1 + t * c2)) / (1. + t * (d1 + t * (d2 + t
            * d3)))
    return sign * g


# ========================

def _erfcc(x):
    """Complementary error function."""

    z = abs(x)
    t = 1. / (1. + 0.5 * z)
    r = t * math.exp(-z * z - 1.26551223 + t * (1.00002368 + t
                     * (.37409196 + t * (.09678418 + t * (-.18628806
                     + t * (.27886807 + t * (-1.13520398 + t
                     * (1.48851587 + t * (-.82215223 + t
                     * .17087277)))))))))
    if x >= 0.:
        return r
    else:
        return 2. - r


# ========================

def _ncdf(x):
    """Cumulative normal dist'n."""

    global _rt2
    return 1. - 0.5 * _erfcc(x / _rt2)


# ========================

def _ncdf_sig(nsig):
    """Cummulative normal dist'n inside nsig sigmas.
....ncdf_sig = 1 - 2 * (upper tail) = 1 - erfc(sigfac/rt(2))"""

    global _rt2
    return 1. - _erfcc(nsig / _rt2)


# =============== Chi squared dist'n ==============

def _pchisq(chisq, nu):
    """Lower tail area of the chi**2 dist'n with nu dof.
....Note that chisq is *not* the reduced chis**2!"""

    hnu = 0.5 * nu
    hchi = 0.5 * chisq
    return _gammp(hnu, hchi)


# ========================

def _qchisq(chisq, nu):
    """Upper tail area of the chi**2 dist'n with nu dof.
....Note that chisq is *not* the reduced chis**2!"""

    hnu = 0.5 * nu
    hchi = 0.5 * chisq
    return _gammq(hnu, hchi)


# ========================

def _chisq_crit(nu, p, tol=1.e-5):
    """Critical chi**2 with lower tail area of p for nu dof."""

    #  For the first guess, use the assyptotic normal limit of the
    # chi**2 distribution:  chi**2 ~ N(nu,sqrt(2*nu)).

    chi = nu + _ncdf_inv(p) * sqrt(2. * nu)
    pcur = _pchisq(chi, nu)

    # Now do a Newton-Raphson loop...

    while 1:
        dfdc = (_pchisq(1.001 * chi, nu) - pcur) / (1e-3 * chi)
        chi = chi - (pcur - p) / dfdc
        pcur = _pchisq(chi, nu)
        if abs(pcur - p) <= tol:
            return chi


##@}
