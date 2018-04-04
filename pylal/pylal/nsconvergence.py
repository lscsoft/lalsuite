#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import string
import sys
from pylal import bayespputils as bppu
from lalapps.combine_evidence import combine_evidence
import os

def get_data_col(data, param_arr, param):
    """
    Return a one-dimensional array with the column 
    of data corresponding to the specified parameter.
    """
    ret = []
    for index, p in enumerate(param_arr):
        if p == param:
            for i in range(len(data)):
                ret.append(data[i][index])
                
    return ret

def merge_files(data, param_arr, writefile):
    """
    Merge a list of files and create a 'chain' column
    that labels each file.
    """

    wf = open(writefile, 'w')

    chainFound = False
    for i in param_arr:
        wf.write(i+' ')
        if i == 'chain':
            chainFound = True
    if not chainFound:
        wf.write('chain\n')
    run_num = 0
    for d in data:
        for i in d:
            for j in i:
                wf.write(str(j)+' ')
            wf.write(str(run_num)+'\n')
        run_num+=1
    wf.close()


def kstest(pos_samples, param_arr, outdir):
    """
    Compute p-value for each parameter in param_arr for each combination
    of runs in pos_samples and histogram results. 
    """
    runs = len(pos_samples)
    num_params = len(param_arr)
    ksdir = os.path.join(outdir,'ks')
    if not os.path.isdir(ksdir):
        os.makedirs(ksdir)

    k = 1
    ks_arr = []
    for param in param_arr:
        D_arr = []
        p_arr = []
        p_plot_arr = []
        for i in range(runs):
            data1 = get_data_col(pos_samples[i], param_arr, param)
            for j in range(runs):
                data2 = get_data_col(pos_samples[j], param_arr, param)
                D, p = stats.ks_2samp(data1, data2)
                D_arr.append(D)
                p_arr.append(p)
                if i is not j:
                    p_plot_arr.append(p)    
        
        ks_arr.append(p_arr)
        plt.figure(k)
        
        plt.hist(p_plot_arr, bins = 20, normed = True)
        plt.xlabel('p-value')
        plt.ylabel('probability density')
        plt.title(param)
        plt.savefig(outdir+'/ks/'+param+'_ks')
        k+=1
    return ks_arr
    
def gelman_rubin(pos_samples, param_arr, outdir):
    """
    Compute Gelman-Rubin R-statistic for each parameter in param_arr.
    """
    writefile = os.path.join(outdir,'merged_files.dat')
    runs = len(pos_samples)
    R_arr = []
    merge_files(pos_samples, param_arr, writefile)
    for param in param_arr:
        data=bppu.PEOutputParser('common')
        inputFileObj = open(writefile)
        dataObj0 = data.parse(inputFileObj)
        posterior = bppu.Posterior(dataObj0)
        R = posterior.gelman_rubin(param)
        R_arr.append(R)
    return R_arr

def compare_maxl(pos_samples, param_arr):
    """
    Find maximum value for logl for each run and compute maximum
    difference.
    """
    runs = len(pos_samples)
    maxl_arr = []
    for i in range(runs):
        l = get_data_col(pos_samples[i], param_arr, 'logl')
        if not l:
            l = get_data_col(pos_samples[i], param_arr, 'likelihood')
        maxl = max(l)
        maxl_arr.append(maxl)
    max_diff = max(maxl_arr)-min(maxl_arr)
    return (maxl_arr, max_diff)

def compare_maxposterior(pos_samples, param_arr):
    """
    Compute the maximum posterior value (maximum of logl+logprior)
    for each run. Calculate maximum difference in maxposterior.
    """
    runs = len(pos_samples)
    l = []
    p = []
    for i in range(runs):
        l.append(get_data_col(pos_samples[i], param_arr, 'logl'))
        p.append(get_data_col(pos_samples[i], param_arr, 'prior'))

    max_pos_arr = []
    for i in range(runs):
        pos = []
        for j in range(len(pos_samples[i])):
            pos.append(l[i][j]+p[i][j])
        max_pos_arr.append(max(pos))   
    max_diff = max(max_pos_arr)-min(max_pos_arr)
    return (max_pos_arr, max_diff)

