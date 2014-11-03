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

description = """This program runs SVM evaluation tasks."""

parser = OptionParser(version='Name: %%prog\n%s'
                      % git_version.verbose_msg, usage='%prog [options]'
                      , description=description)
parser.add_option('-s', '--scale', dest='svmscale',
                  help='svm scale cmd line')
parser.add_option('-p', '--predict', dest='svmpredict',
                  help='svm predict cmd line')
parser.add_option('-b', '--rank', dest='rankmode',
                  help='svm rank output mode')
parser.add_option('-i', '--input-file', dest='testfile',
                  help='input data file')
parser.add_option('-r', '--range-file', dest='rangefile',
                  help='range file for svm scale')
parser.add_option('-m', '--model-file', dest='svmmodel',
                  help='svm model file')
parser.add_option('-o', '--output-file', dest='patfile',
                  help='output data file')

(opts, files) = parser.parse_args()
assert os.path.exists(opts.svmscale), 'svm-scale executable not found'
assert os.path.exists(opts.svmpredict), \
    'svm-predict executable not found'

b = (1 if opts.rankmode == 'prob' else 0)

test_file = './' + os.path.split(opts.testfile)[1] + '.mid'
svmkit.ConvertLIGOtoLibSVM(opts.testfile, test_file)
scale_test_file = './' + os.path.split(test_file)[1] + '.scale'
predict_file = './' + os.path.split(test_file)[1] + '.predict'

cmd_scale = opts.svmscale + ' -r ' + opts.rangefile + ' ' + test_file \
    + ' > ' + scale_test_file
subprocess.Popen(cmd_scale, shell=True, stderr=subprocess.PIPE,
                 stdout=subprocess.PIPE).communicate()
cmd_predict = opts.svmpredict + ' -b ' + str(b) + ' ' + scale_test_file \
    + ' ' + opts.svmmodel + ' ' + predict_file
subprocess.Popen(cmd_predict, shell=True, stderr=subprocess.PIPE,
                 stdout=subprocess.PIPE).communicate()

svmkit.ConvertLibSVMtoLIGO(opts.testfile, predict_file, opts.patfile)
cmd_rm = 'rm ' + ' ' + scale_test_file + ' ' + test_file + ' ' \
    + predict_file
subprocess.Popen(cmd_rm, shell=True)
