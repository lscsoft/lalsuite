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

import sys
import os
from optparse import OptionParser
import ConfigParser
import numpy as np
from laldetchar.idq import svmkit
import subprocess

from laldetchar import git_version

__author__ = 'Yingsheng Ji <yingsheng.ji@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

description = """This program runs SVM training tasks."""

parser = OptionParser(version='Name: %%prog\n%s'
                      % git_version.verbose_msg, usage='%prog [options]'
                      , description=description)
parser.add_option('', '--train', dest='svmtrain',
                  help='svm train cmd line')
parser.add_option('', '--scale', dest='svmscale',
                  help='svm scale cmd line')
parser.add_option('', '--rank', dest='rankmode',
                  help='svm rank output mode')
parser.add_option('', '--train-file', dest='trainfile',
                  help='train data file in original format')
parser.add_option('', '--train-file-svm', dest='trainfilesvm',
                  help='train data file in Libsvm format')
parser.add_option('', '--scale-file', dest='scalefile',
                  help='scale data file')
parser.add_option('', '--range-file', dest='rangefile',
                  help='output range file for svm-scale')
parser.add_option('', '--model-file', dest='modelfile',
                  help='output model file for svm-predict')
parser.add_option('-g', '--gamma', dest='gamma',
                  help='svm parameter for svm 2-class train with RBF')
parser.add_option('-c', '--cost', dest='cost',
                  help='svm parameter for svm 2-class train')

(opts, files) = parser.parse_args()
assert os.path.exists(opts.svmscale), 'svm-scale executable not found'
assert os.path.exists(opts.svmtrain), 'svm-train executable not found'

b = (1 if opts.rankmode == 'prob' else 0)

svmkit.ConvertLIGOtoLibSVM(opts.trainfile, opts.trainfilesvm)

cmd_scale = opts.svmscale + ' -s ' + opts.rangefile + ' ' \
    + opts.trainfilesvm + ' > ' + opts.scalefile
subprocess.Popen(cmd_scale, shell=True, stderr=subprocess.PIPE,
                 stdout=subprocess.PIPE).communicate()

cmd_train = opts.svmtrain + ' -b ' + str(b) + ' -g ' + opts.gamma \
    + ' -c ' + opts.cost + ' ' + opts.scalefile + ' ' + opts.modelfile
subprocess.Popen(cmd_train, shell=True, stderr=subprocess.PIPE,
                 stdout=subprocess.PIPE).communicate()

