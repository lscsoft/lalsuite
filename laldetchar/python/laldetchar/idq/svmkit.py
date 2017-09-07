# Copyright (C) 2013 Yingsheng Ji
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the# Free Software Foundation; either version 3 of the License, or (at your# option
# ) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## \defgroup laldetchar_py_idq_svmkit SVM Utils
## \ingroup laldetchar_py_idq
## Synopsis
# ~~~
# from laldetchar.idq import svmkit
# ~~~
# \author Yingsheng Ji (<yingsheng.ji@ligo.org>)

"""Module with utility functions for SVM algorithm"""

import sys
import numpy as np
from laldetchar.idq import auxmvc_utils
from laldetchar import git_version

__author__ = 'Yingsheng Ji <yingsheng.ji@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date


## \addtogroup laldetchar_py_idq_svmkit
# @{

# Utility functions for SVM algorithm.

def ConvertLIGOtoLibSVM(patfile, svmfile):
    """
....This module is used to convert LIGO data format into libSVM format data file
...."""

    # samples = ReadLIGOfile([patfile])

    samples = auxmvc_utils.ReadMVSCTriggers([patfile], Classified=False)
    features = list(samples.dtype.names)
    for name in [
        'GPS_s',
        'GPS_ms',
        'signif',
        'SNR',
        'unclean',
        'i',
        ]:
        features.remove(name)
    file = open(svmfile, 'w')

    for s in samples:
        line = ''.join([(' ' + str(i + 1) + ':' + str(s[var]) if s[var]
                       != 0 else '') for (i, var) in
                       enumerate(features)])
        file.write((('+1' if s['i'] == 1 else '-1')) + line + '\n')

    file.close()


def ConvertLibSVMtoLIGO(patfile, rstfile, datfile):
    """
....This module is used to convert SVM predict result file into LIGO data format.
...."""

    # features = open(patfile).readlines()[1].strip('\n')+" SVMrank"
    # samples = ReadLIGOfile([patfile])

    samples = auxmvc_utils.ReadMVSCTriggers([patfile], Classified=False)

    # ranks = np.loadtxt(rstfile, usecols=(0,),skiprows=0)

    ranks = np.loadtxt(rstfile, usecols=(1, ), skiprows=1)
    if not ranks.shape:  # single row loaded, requires array reshaping
        ranks = ranks.reshape((1, ))

    features = list(samples.dtype.names)
    for f in [
        'GPS_s',
        'GPS_ms',
        'i',
        'unclean',
        'signif',
        'SNR',
        ]:
        features.remove(f)
    file = open(datfile, 'w')
    head = [
        'GPS',
        'i',
        'unclean',
        'signif',
        'SNR',
        'rank',
        ] + features
    file.write(' '.join([str(h) for h in head]) + '\n')

    for (j, s) in enumerate(samples):
        GPS = float(s['GPS_s']) + float(s['GPS_ms']) / 1000
        i = float(s['i'])
        prefix = ' '.join(str(s[f]) for f in ['unclean', 'signif', 'SNR'
                          ])
        suffix = ' '.join([str(s[f]) for f in features])
        file.write(str('%0.3f' % GPS) + ' ' + str(i) + ' ' + prefix
                   + ' ' + str(ranks[j]) + ' ' + suffix + '\n')
    file.close()


def ReadLIGOfile(files):
    """
    Read LIGO data *.pat file
    """

    for (i, f) in enumerate(files):
        features = open(f).readlines()[1].split()
        features.append('i')
        formats = ['g8' for v in range(len(features))]
        for v in ['GPS_s', 'GPS_ms', 'i']:
            formats[features.index(v)] = 'i'
            if i > 0:
                samples = np.concatenate((samples, np.loadtxt(f,
                        skiprows=2, dtype={'names': features,
                        'formats': formats})), axis=0)
            else:
                samples = np.loadtxt(f, skiprows=2,
                        dtype={'names': features, 'formats': formats})
    return samples


##@}
