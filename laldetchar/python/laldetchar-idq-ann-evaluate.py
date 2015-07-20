# Copyright (C) 2015 Young-Min Kim
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


from pyfann import libfann
from laldetchar.idq import auxmvc_utils
from laldetchar.idq import auxmvc
from optparse import OptionParser
from optparse import *
import glob
import sys
import matplotlib
import pdb
import numpy
from laldetchar import git_version

__author__ = 'Young-Min Kim <young-min.kim@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

description = """This program runs ANN evaluation tasks."""

parser = OptionParser(version='Name: %%prog\n%s'%git_version.verbose_msg, 
                                usage='%prog [options]', 
                                description=description)
parser.add_option("-v","--verbose",action="store_true",help="verbose mode for printing run processing.")
parser.add_option("-e","--evaluation-file",action="store",type="string",help="Trigger file name to be used for evaluation")
parser.add_option("-n","--network",action="store",type="string",help="Network file name to be used for evaluation ")
parser.add_option("-s","--saving-results",action="store",type="string",default=False,help="saving results")
parser.add_option("","--training-machine",action="store",type="string",default='fann',help="saving results")
(opts,files)=parser.parse_args()


print "Reading evaluation data : %s" % opts.evaluation_file
test_ann_data = libfann.training_data()
test_pat_data= auxmvc_utils.ReadMVSCTriggers([opts.evaluation_file],Classified=False)

test_ann_data.read_train_from_file(opts.evaluation_file.replace(".pat",".ann"))
test_input = test_ann_data.get_input()
test_output = test_ann_data.get_output()

network = opts.network

print "Importing network file : %s" % network
ann = libfann.neural_net()
ann.create_from_file(network)

#ann.print_parameters()
ann.reset_MSE()

if opts.verbose:
    print "Start evaluating data"

dat_vars = ['GPS','i','w','unclean','signif','SNR','rank'] + list(test_pat_data.dtype.names)[5:-1]
dat_formats = []
for va in dat_vars:
    if va in ['i','unclean']:
        dat_formats.append('i')
    else:
        dat_formats.append('g8')

rankedTriggers = numpy.empty((len(test_pat_data),), dtype={'names':dat_vars,'formats':dat_formats})

for va in dat_vars:
    if va not in['GPS','w','rank']:
        rankedTriggers[va] = test_pat_data[va]

rankedTriggers['GPS'] = test_pat_data['GPS_s'] + test_pat_data['GPS_ms'] * 10 ** -3
rankedTriggers['w'] = numpy.ones(len(test_pat_data))

for i in range(len(test_input)):
    #calculating glitch-rank
    glitch_rank = ann.run(test_input[i])[0]
    # push the rank to rankedTrigger matrix
    #result_file.write(str(patline[0]).split(".")[0]+"."+str(patline[1]*0.001).split(".")[1]+" "+str(test_output[i][0])+" "+str(w_row[i])+" "+str(patline[4])+" "+str(patline[2])+" "+str(patline[3])+" "+str(results[0])+" "+variable_line+"\n")
    rankedTriggers['rank'][i] = glitch_rank

output_filename = opts.saving_results
auxmvc_utils.WriteMVSCTriggers(rankedTriggers, output_filename, Classified=True)

#print "MSE on evaluation data: %f" % ann.get_MSE()
if opts.verbose:
    print "Evaluation results are saved in %s." % output_filename

