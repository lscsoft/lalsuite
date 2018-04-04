#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#       cbcBayesConvergence.py
#
#       Copyright 2012
#       Alex Mellus <flyers585@gmail.com>
#
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

#===============================================================================
# Preamble
#===============================================================================

#standard library imports
import sys
import os

from math import ceil,floor
import cPickle as pickle

from time import strftime

#related third party imports
from numpy import array,exp,cos,sin,arcsin,arccos,sqrt,size,mean,column_stack,cov,unique,hsplit,correlate,log,dot,power,squeeze,sort
from scipy import stats

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

#local application/library specific imports
from pylal import bayespputils as bppu
from pylal import git_version

__author__="Ben Aylott <benjamin.aylott@ligo.org>, Ben Farr <bfarr@u.northwestern.edu>, Will M. Farr <will.farr@ligo.org>, John Veitch <john.veitch@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

def pickle_to_file(obj,fname):
    """
    Pickle/serialize 'obj' into 'fname'.
    """
    filed=open(fname,'w')
    pickle.dump(obj,filed)
    filed.close()

def oneD_dict_to_file(dict,fname):
    filed=open(fname,'w')
    for key,value in dict.items():
        filed.write("%s %s\n"%(str(key),str(value)) )

def multipleFileCB(opt, opt_str, value, parser):
    args=[]
    for arg in parser.rargs:
        if arg[0] != "-":
            args.append(arg)
        else:
            del parser.rargs[:len(args)]
            break
    #Append new files to list if some already specified
    if getattr(parser.values, opt.dest):
        args.extend(getattr(parser.values, opt.dest))
    setattr(parser.values, opt.dest, args)

def cbcBayesConvergence(
                        outdir,data,
                        #Number of live points used in data
                        ns_Nlive,
                        #Threshold for failure for gelman-rubin test
                        gelmanthresh=1.01
                        #
                    ):
    """
    This is a script which calculates convergence statistics for nested 
    sampling runs.
    """

    if data is None:
        raise RuntimeError('You must specify an input data file')
    #
    if outdir is None:
        raise RuntimeError("You must specify an output directory.")

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    #

    import string
    from numpy import loadtxt
    pos_samples = []
    new_data = []
    for d in reversed(data):
        temp = [d]
        new_data.append(temp)

    peparser=bppu.PEOutputParser('ns')
    posfilename=os.path.join(outdir,'posterior_samples.dat')
    
    for i in range(len(data)):
        # Create a posterior object (Npost=None avoids repeated samples which ruin the KS test)
        commonResultsObj=peparser.parse(new_data[i],Nlive=ns_Nlive,Npost=None)
        pos = bppu.Posterior(commonResultsObj)
        pos.write_to_file(posfilename)
    
        with open(posfilename) as f:
            param_arr = string.split(f.readline())
            loadfile = loadtxt(f)
            pos_samples.append(loadfile)
   
    #==================================================================#
    #Create webpage with nested sampling convergence information
    #==================================================================#
    
    import pylal.nsconvergence as nsc
    runs = len(pos_samples)
    html_nsconvergence=bppu.htmlPage('Convergence Information', css=bppu.__default_css_string)

    convergencedir = os.path.join(outdir, 'convergence')
    if not os.path.isdir(convergencedir):
        os.makedirs(convergencedir)    	

    #Summary Section
    html_nsconvergence_stats=html_nsconvergence.add_section('Summary')
    max_l, l_diff = nsc.compare_maxl(pos_samples, param_arr)
    html_nsconvergence_stats.p('Max difference in loglikelihood: %f'%l_diff)
    summary_table_string = ''
    summary_table_header = '<table border="1"><tr><th>Run</th><th>maxloglikelihood</th>'
    #maxposterior column
    if param_arr.count('prior') > 0:
        max_pos, pos_diff = nsc.compare_maxposterior(pos_samples, param_arr)
        html_nsconvergence_stats.p('Max difference in posterior: %f'%pos_diff)
        summary_table_header += '<th>maxposterior</th>'	

    summary_table_header += '</tr>'
    summary_table_string += summary_table_header
    for i in range(runs):
        max_l_val = max_l[i]
        summary_table_string += '<tr><td>%i</td><td>%f</td>'%(i,max_l_val)
        if param_arr.count('prior') > 0:
            max_pos_val = max_pos[i]
            summary_table_string += '<td>%f</td>'%max_pos_val
        summary_table_string += '</tr>'
        
    summary_table_string += '</table>'
    html_nsconvergence_stats.write(summary_table_string)
    
    #KS Test Section
    html_nsconvergence_ks=html_nsconvergence.add_section('KS Test')
    ks_arr = nsc.kstest(pos_samples, param_arr, convergencedir)
    for index, p in enumerate(param_arr):	
        ks_table_string = '<table><caption>%s</caption><tr><th></th>'%p
        for i in range(runs):
            ks_table_string += '<th>Run%i</th>'%i 
        ks_table_string+='</tr>'
        for i in range(runs):
            ks_table_string += '<tr><th>Run%i</th>'%i
            for j in range(i*runs,(i*runs)+runs):
                pval = ks_arr[index][j]	
                ks_table_string += '<td>%f</td>'%pval
            ks_table_string += '</tr>'	
        ks_table_string += '</table>'
        html_nsconvergence_ks.write(ks_table_string)
    for p in param_arr:
        html_nsconvergence_ks.write('<a href="./convergence/ks/'+p+'_ks.png" target="_blank"><img width="35%" src="./convergence/ks/'+p+'_ks.png"/></a>')

    #Gelman Rubin Section
    html_nsconvergence_gelman=html_nsconvergence.add_section('Gelman Rubin')
    
    gelmandir = os.path.join(convergencedir,'gelmanrubin')
    if not os.path.isdir(gelmandir):
        os.makedirs(gelmandir)
    
    gelmanrubin = nsc.gelman_rubin(pos_samples, param_arr, gelmandir)
    warn = False
    warnparams = []
    for index,g in enumerate(gelmanrubin):
        if g > gelmanthresh:
            warn = True
            warnparams.append(index)
    if warn:
        with open(outdir+'/warning.txt', 'w') as warnfile:
            warnfile.write('Gelman-Rubin threshold set to %f\n'%gelmanthresh)
            for i in warnparams:
                warnfile.write('%s has an R-value of %f\n'%(param_arr[i], gelmanrubin[i]))
                    
    colors = ['b', 'r', 'g', 'c', 'k', 'm', 'y', .25, .5, .75]
    for param_index, param in enumerate(param_arr):
        for i in range(runs):
            data_range = []
            for j in range(len(pos_samples[i])):
                data_range.append(j) 
            col = nsc.get_data_col(pos_samples[i], param_arr, param)
            plt.figure(param_index)
            plt.scatter(data_range, col, c = colors[i], s = 5, edgecolors = 'none')
            plt.title('R = ' + str(gelmanrubin[param_index]))
            plt.xlabel('Sample')
            plt.ylabel(param)	
            plt.xlim(0,len(pos_samples[i]))	
            plt.ylim(min(col),max(col))
            plt.savefig(gelmandir+'/'+param)
        
    for p in param_arr:
        html_nsconvergence_gelman.write('<a href="./convergence/gelmanrubin/'+p+'.png" target="_blank"><img width="35%" src="./convergence/gelmanrubin/'+p+'.png"/></a>')


    #Write convergence page
    nsconvergencepage=open(os.path.join(outdir, 'convergence.html'), 'w')
    nsconvergencepage.write(str(html_nsconvergence))
    nsconvergencepage.close()

USAGE='''%prog [options] datafile.dat [datafile2.dat ...]
Generate a web page displaying results of parameter estimation based on the contents
of one or more datafiles containing samples from one of the bayesian algorithms (MCMC, nested sampling).
Options specify which extra statistics to compute and allow specification of additional information.
'''

if __name__=='__main__':

    from optparse import OptionParser
    parser=OptionParser(USAGE)
    parser.add_option("-o","--outpath", dest="outpath",help="make page and plots in DIR", metavar="DIR")
    parser.add_option("-d","--data",dest="data",action="callback",callback=multipleFileCB,help="datafile")
    parser.add_option("--gelmanthresh", action="store", type="float", dest="gelmanthresh", default=1.01, help="set threshold for failure for gelman-rubin test")
    parser.add_option("--Nlive",action="store",default=None,help="(inspnest) Number of live points used in each parallel nested sampling run.",type="int")
    
    (opts,args)=parser.parse_args()

    datafiles=[]
    if args:
      datafiles=datafiles+args
    if opts.data:
      datafiles=datafiles + opts.data
    
    cbcBayesConvergence(
                        opts.outpath,datafiles,                        
                        opts.Nlive,
                        #Threshold for failure for gelman-rubin test
                        gelmanthresh=opts.gelmanthresh
                    )

