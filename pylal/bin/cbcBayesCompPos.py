#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcBayesCompPos.py Copyright 2010--2012 Benjamin Aylott
#       <benjamin.aylott@ligo.org>, Will M. Farr <will.farr@ligo.org>
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

#standard library imports
import os
import sys
from time import strftime
import copy
import random
import getpass

#related third party imports
import numpy as np
from numpy import floor

import scipy.stats as ss

import matplotlib as mpl
mpl.use("AGG")
from matplotlib import pyplot as plt
from matplotlib import colors as mpl_colors
from matplotlib import cm as mpl_cm
from matplotlib.ticker import FormatStrFormatter,ScalarFormatter,AutoMinorLocator

from scipy import stats

#local application/library specific imports
import pylal.bayespputils as bppu
from pylal import SimInspiralUtils
from pylal import git_version

__author__="Ben Aylott <benjamin.aylott@ligo.org>, Will M. Farr <will.farr@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

#List of parameters to plot/bin . Need to match (converted) column names.
oneDMenu=['mtotal','m1','m2','mchirp','mc','chirpmass','distance','distMPC','dist','iota','psi','eta','q','asym_massratio','spin1','spin2','a1','a2','phi1','theta1','phi2','theta2','costilt1','costilt2','costhetas','cosbeta','phi_orb', 'lambdat', 'dlambdat', 'lambda1', 'lambda2', 'lam_tilde', 'dlam_tilde','theta_jn','a1z','a2z'] + bppu.snrParams
#List of parameter pairs to bin . Need to match (converted) column names.
twoDGreedyMenu=[['mc','eta'],['mchirp','eta'],['chirpmass','eta'],['mc','q'],['mchirp','q'],['chirpmass','q'],['mc','asym_massratio'],['mchirp','asym_massratio'],['chirpmass','asym_massratio'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['dist','m1'],['ra','dec'],['dist','cos(iota)'],['phi_orb','iota'],['theta_jn','dist'],['spin1','spin2'],['spin1','mchirp'],['spin1','m1'],['a1','a2'],['a1','mchirp'],['a1','m1'],['tilt1','tilt2'],['tilt1','mchirp'],['tilt1','m1'],['a1z','a2z']]
#Bin size/resolution for binning. Need to match (converted) column names.

#Convert parameter names to LaTeX; if not mentioned here, just use parameter name.
paramNameLatexMap = {'m1': 'm_1', 'm2' : 'm_2', 'mtotal' : r'M_{\rm tot}', 'mchirp' : r'\mathcal{M}',
                     'mc': r'\mathcal{M}', 'distance' : 'd', 'distMPC' : 'd', 'dist': 'd',
                     'iota': r'\iota', 'psi': '\psi', 'eta': '\eta', 'asym_massratio': 'q', 'a1': 'a_1',
                     'a2': 'a_2', 'phi1': r'\phi_1', 'phi2': r'\phi_2', 'theta1': r'\theta_1', 'theta2': r'\theta_2',
                     'cos(tilt1)': r'\cos t_1', 'cos(tilt2)': r'\cos t_2', 'cos(thetas)': r'\cos \theta_s',
                     'cosbeta': r'\cos \beta', 'phi_orb': r'\phi_{\rm orb}', 'cos(beta)': r'\cos \beta',
                     'cos(iota)': r'\cos \iota', 'tilt1': r't_1', 'tilt2': r't_2', 'ra': r'\alpha', 'dec': r'\delta',
                     'lambdat' : r'\tilde{\Lambda}', 'dlambdat': r'\delta \tilde{\Lambda}',
                     'lambda1' : r'\lambda_1', 'lambda2': r'\lambda_2',
                     'lam_tilde' : r'\tilde{\Lambda}', 'dlam_tilde': r'\delta \tilde{\Lambda}','dchi0':r'\delta\chi_0','dchi1':r'\delta\chi_1','dchi2':r'\delta\chi_2','dchi3':r'\delta\chi_3','dchi4':r'\delta\chi_4','dchi5':r'\delta\chi_5','dchi5l':r'\delta\chi_{5l}','dchi6':r'\delta\chi_6','dchi6l':r'\delta\chi_{6l}','dchi7':r'\delta\chi_7','dbeta2':r'\delta\beta_2','dbeta3':r'\delta\beta_3','dsigma2':r'\delta\sigma_2','dsigma3':r'\delta\sigma_3','dsigma4':r'\delta\sigma_4','dbeta2':r'\delta\beta_2','dbeta3':r'\delta\beta_3' }

# Only these parameters, in this order appear in confidence level table.
clTableParams = ['mchirp', 'mc', 'chirpmass', 'eta', 'q', 'm1', 'm2', 'distance', 'distMPC', 'dist', 'cos(iota)', 'iota', 'theta_jn', 'psi', 'ra', 'dec', 'time', 'phase', 'a1', 'a2', 'costilt1', 'costilt2','dchi0','dchi1','dchi2','dchi3','dchi4','dchi5','dchi5l','dchi6','dchi6l','dchi7','dbeta2','dbeta3','dsigma2','dsigma3','dsigma4','dbeta2','dbeta3']


greedyBinSizes={'mc':0.001,'m1':0.1,'m2':0.1,'mass1':0.1,'mass2':0.1,'mtotal':0.1,'eta':0.001,'q':0.001,'asym_massratio':0.001,'iota':0.05,'time':1e-4,'distance':5.0,'dist':1.0,'mchirp':0.01,'chirpmass':0.01,'a1':0.02,'a2':0.02,'phi1':0.05,'phi2':0.05,'theta1':0.05,'theta2':0.05,'ra':0.05,'dec':0.005,'psi':0.1,'cos(iota)':0.01, 'cos(tilt1)':0.01, 'cos(tilt2)':0.01, 'tilt1':0.05, 'tilt2':0.05, 'cos(thetas)':0.01, 'cos(beta)':0.01,'phi_orb':0.2,'inclination':0.05,'theta_jn':0.05,'spin1':0.02,'spin2':0.02}
for s in bppu.snrParams:
        greedyBinSizes[s]=0.02
for s in bppu.calParams:
        greedyBinSizes[s]=0.02
for s in bppu.tigerParams:
  greedyBinSizes[s]=0.01
#Confidence levels
OneDconfidenceLevels=[0.9]#[0.68,0.9,0.95,0.997]
TwoDconfidenceLevels=OneDconfidenceLevels

#2D plots list
#twoDplots=[['mc','eta'],['mchirp','eta'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['RA','dec'],['ra','dec'],['m1','dist'],['m2','dist'],['psi','iota'],['psi','distance'],['psi','dist'],['psi','phi0'],['dist','cos(iota)']]
twoDplots=[['m1','m2'],['mass1','mass2'],['RA','dec'],['ra','dec'],['cos(thetas)','cos(beta)'],['distance','iota'],['dist','iota'],['dist','cosiota'],['distance','cosiota'],['psi','iota'],['psi','distance'],['psi','phi0'],['dist','cos(iota)'],['phi_orb','iota'],['distance','inclination'],['dist','inclination'],['theta_jn','dist'],['spin1','spin2'],['spin1','mchirp'],['spin1','m1'],['a1','a2'],['a1','mchirp'],['a1','m1'],['tilt1','tilt2'],['tilt1','mchirp'],['tilt1','m1']]
allowed_params=['mtotal','m1','m2','mchirp','mc','chirpmass','q','asym_massratio','distance','distMPC','dist','iota','psi','eta','ra','dec','a1','a2','spin1','spin2','phi1','theta1','phi2','theta2','cos(iota)','cos(tilt1)','cos(tilt2)','tilt1','tilt2','cos(thetas)','cos(beta)','phi_orb','inclination', 'logl', 'lambdat', 'dlambdat', 'lambda1', 'lambda2', 'lam_tilde', 'dlam_tilde','theta_jn','a1z','a2z']+bppu.snrParams+bppu.calParams

def open_url(url,username,password):

    import urllib
    import urllib2
    import urlparse

    parsed_url=urlparse.urlparse(url)
    url=urlparse.urljoin(parsed_url.geturl(),'posterior_samples.dat')


    opener = urllib2.build_opener(urllib2.HTTPCookieProcessor())
    urllib2.install_opener(opener)

    body={'username':username,'password':password}
    txdata = urllib.urlencode(body) # if we were making a POST type request, we could encode a dictionary of values here - using urllib.urlencode
    txheaders = {'User-agent' : 'Mozilla/4.0 (compatible; MSIE 5.5; Windows NT)'} # fake a user agent, some websites (like google) don't like automated exploration

    req = urllib2.Request(parsed_url[0]+'://'+parsed_url[1], txdata, txheaders)

    resp = opener.open(req) # save a cookie
    dump=resp.read()
    resp.close()
    try:
        req = urllib2.Request(url, txdata, txheaders) # create a request object
        handle = opener.open(req) # and open it to return a handle on the url
        data = handle.read()
        f = file('posterior_samples.dat', 'w')
        f.write(data)
        f.close()

    except IOError, e:
        print 'We failed to open "%s".' % url
        if hasattr(e, 'code'):
            print 'We failed with error code - %s.' % e.code
        elif hasattr(e, 'reason'):
            print "The error object has the following 'reason' attribute :", e.reason
            print "This usually means the server doesn't exist, is down, or we don't have an internet connection."
            sys.exit()
    else:
        print 'Here are the headers of the page :'
        print handle.info() # handle.read() returns the page, handle.geturl() returns the true url of the page fetched (in case urlopen has followed any redirects, which it sometimes does)

    return HTMLSource

def all_pairs(L):
    while L:
        i = L.pop()
        for j in L: yield i, j

def open_url_curl(url,args=[]):
    import subprocess
    import urlparse

    kerberos_args = "--insecure -c /tmp/{0}_cookies -b /tmp/{0}_cookies --negotiate --user : --location-trusted".format(getpass.getuser()).split()

    retcode=subprocess.call(['curl'] + kerberos_args + [url] + args)

    return retcode

def test_and_switch_param(common_output_table_header,test,switch):
    try:
        idx=common_output_table_header.index(test)
        common_output_table_header[idx]=switch
        print "Re-labelled %s -> %s"%(test,switch)
    except:
        pass

    return

def compare_plots_one_param_pdf(list_of_pos_by_name,param,analyicPDF=None):
    """
    Plots a gaussian kernel density estimate for a set
    of Posteriors onto the same axis.

    @param list_of_pos: a list of Posterior class instances.

    @param plot1DParams: a dict; {paramName:Nbins}

    """

    from scipy import seterr as sp_seterr

    #Create common figure
    myfig=plt.figure(figsize=(6,4.5),dpi=150)

    list_of_pos=list_of_pos_by_name.values()
    list_of_pos_names=list_of_pos_by_name.keys()

    allmins=map(lambda a: np.min(a[param].samples), list_of_pos)
    allmaxes=map(lambda a: np.max(a[param].samples), list_of_pos)
    min_pos=np.min(allmins)
    max_pos=np.max(allmaxes)
    print 'Found global min: %f, max: %f'%(min_pos,max_pos)

    gkdes={}
    injvals=[]
    for name,posterior in list_of_pos_by_name.items():

        pos_samps=posterior[param].samples
        if posterior[param].injval is not None:
            injvals.append(posterior[param].injval)

        min_pos_temp=np.min(pos_samps)
        max_pos_temp=np.max(pos_samps)

        if min_pos_temp<min_pos:
            min_pos=min_pos_temp
        if max_pos_temp>max_pos:
            max_pos=max_pos_temp

        injpar=posterior[param].injval

        gkdes[name]=posterior[param].gaussian_kde

    if gkdes:
        ind=np.linspace(min_pos,max_pos,101)

        kdepdfs=[]
        for name,gkde in gkdes.items():
            kdepdf=gkde.evaluate(ind)
            kdepdfs.append(kdepdf)
            plt.plot(ind,np.transpose(kdepdf),label=name)
        plt.grid()
        plt.legend()
        plt.xlabel(bppu.plot_label(param))
        plt.xlim(min_pos,max_pos)
        plt.ylabel('Probability Density')
        try:
          plt.tight_layout()
        except:
          pass
        if injvals:
            print "Injection parameter is %f"%(float(injvals[0]))
            injpar=injvals[0]
            if min(pos_samps)<injpar and max(pos_samps)>injpar:
                plt.plot([injpar,injpar],[0,max(kdepdf)],'r-.',scalex=False,scaley=False)
    if analyticPDF is not None:
	plt.plot(ind,map(analyticPDF,ind),'r')
    #
    return myfig#,rkde
#
def compare_plots_one_param_line_hist(list_of_pos_by_name,param,cl,color_by_name,cl_lines_flag=True,legend='right',analyticPDF=None):


    """
    Plots a gaussian kernel density estimate for a set
    of Posteriors onto the same axis.

    @param list_of_pos: a list of Posterior class instances.

    @param plot1DParams: a dict; {paramName:Nbins}

    """

    from scipy import seterr as sp_seterr

    #Create common figure
    myfig=plt.figure(figsize=(6,4.5),dpi=150)
  #myfig.add_axes([0.1,0.1,0.65,0.85])
  #myfig.add_axes([0.15,0.15,0.6,0.76])
    axes=plt.Axes(myfig,[0.15,0.15,0.6,0.76])
    myfig.add_axes(axes)
    majorFormatterX=ScalarFormatter(useMathText=True)
    majorFormatterX.format_data=lambda data:'%.6g'%(data)
    majorFormatterY=ScalarFormatter(useMathText=True)
    majorFormatterY.format_data=lambda data:'%.6g'%(data)
    majorFormatterX.set_scientific(True)
    majorFormatterY.set_scientific(True)

    list_of_pos=list_of_pos_by_name.values()
    list_of_pos_names=list_of_pos_by_name.keys()

    allmins=map(lambda a: np.min(a[param].samples), list_of_pos)
    allmaxes=map(lambda a: np.max(a[param].samples), list_of_pos)
    min_pos=np.min(allmins)
    max_pos=np.max(allmaxes)

    injvals=[]

    patch_list=[]
    max_y=0

    posbins=np.linspace(min_pos,max_pos,50)

    for name,posterior in list_of_pos_by_name.items():
        colour=color_by_name[name]
        #myfig.gca(autoscale_on=True)
        if posterior[param].injval:
            injvals.append(posterior[param].injval)

        try:
            n,bins=np.histogram(posterior[param].samples,bins=posbins,normed=True,new=True)
        except:
            n,bins=np.histogram(posterior[param].samples,bins=posbins,normed=True)
        if min(bins)==max(bins):
            print 'Skipping '+param
            continue
        locmaxy=max(n)
        if locmaxy>max_y: max_y=locmaxy
#(n, bins, patches)=plt.hist(posterior[param].samples,bins=bins,facecolor='white',label=name,normed=True,hold=True,color=color_by_name[name])#range=(min_pos,max_pos)
        (n, bins, patches)=plt.hist(posterior[param].samples,bins=bins,histtype='step',label=name,normed=True,hold=True,color=color_by_name[name])
        patch_list.append(patches[0])

    Nchars=max(map(lambda d:len(majorFormatterX.format_data(d)),axes.get_xticks()))
    if Nchars>8:
      Nticks=3
    elif Nchars>5:
      Nticks=4
    elif Nchars>4:
      Nticks=6
    else:
      Nticks=6
    locatorX=mpl.ticker.MaxNLocator(nbins=Nticks)
    locatorX.view_limits(bins[0],bins[-1])
    axes.xaxis.set_major_locator(locatorX)

    plt.xlim(min_pos,max_pos)
    top_cl_intervals_list={}
    pos_names=list_of_pos_by_name.keys()


    for name,posterior in list_of_pos_by_name.items():
        #toppoints,injectionconfidence,reses,injection_area,cl_intervals=bppu.greedy_bin_one_param(posterior,{param:greedyBinSizes[param]},[cl])
        cl_intervals=posterior[param].prob_interval([cl])
        colour=color_by_name[name]
        if cl_intervals[0] is not None and cl_lines_flag:
            try:
                plt.plot([cl_intervals[0][0],cl_intervals[0][0]],[0,max_y],color=colour,linestyle='--')
                plt.plot([cl_intervals[0][1],cl_intervals[0][1]],[0,max_y],color=colour,linestyle='--')
            except:
                print "MAX_Y",max_y,[cl_intervals[0][0],cl_intervals[0][0]],[cl_intervals[0][1],cl_intervals[0][1]]
        top_cl_intervals_list[name]=(cl_intervals[0][0],cl_intervals[0][1])

    if cl_lines_flag:
        pos_names.append(str(int(cl*100))+'%')
        patch_list.append(mpl.lines.Line2D(np.array([0.,1.]),np.array([0.,1.]),linestyle='--',color='black'))

    plt.grid()
    plt.xlim(min_pos,max_pos)
    if legend is not None:
      oned_legend=plt.figlegend(patch_list,pos_names,'right')
      for text in oned_legend.get_texts():
        text.set_fontsize('small')
    plt.xlabel(bppu.plot_label(param))
    plt.ylabel('Probability density')
    plt.draw()
    #plt.tight_layout()
    if injvals:
        print "Injection parameter is %f"%(float(injvals[0]))
        injpar=injvals[0]
        #if min(pos_samps)<injpar and max(pos_samps)>injpar:
        plt.plot([injpar,injpar],[0,max_y],'r-.',scalex=False,scaley=False,linewidth=4,label='Injection')

    #
    if analyticPDF is not None:
	plt.plot(posbins,map(analyticPDF,posbins),'r')
    return myfig,top_cl_intervals_list#,rkde

#
def compute_ks_pvalue_matrix(list_of_pos_by_name, param):
    """Returns a matrix of ks p-value tests between pairs of
    posteriors on the 1D marginalized distributions for param."""

    poss=list_of_pos_by_name.values()

    N=len(poss)

    matrix=np.zeros((N,N))
    matrix[:,:]=float('nan')

    for i in range(N):
        pi=poss[i]
        for j in range(i+1, N):
            pj=poss[j]

            d,pvalue=ss.ks_2samp(pi[param].samples.flatten(), pj[param].samples.flatten())

            matrix[i,j]=pvalue
            matrix[j,i]=pvalue

    return matrix
        
def compare_plots_one_param_line_hist_cum(list_of_pos_by_name,param,cl,color_by_name,cl_lines_flag=True,analyticCDF=None,legend='auto'):

    """
    Plots a gaussian kernel density estimate for a set
    of Posteriors onto the same axis.

    @param list_of_pos: a list of Posterior class instances.

    @param plot1DParams: a dict; {paramName:Nbins}

    """

    from scipy import seterr as sp_seterr

    #Create common figure
    myfig=plt.figure(figsize=(6,4.5),dpi=150)
    myfig.add_axes([0.15,0.15,0.6,0.76])
    list_of_pos=list_of_pos_by_name.values()
    list_of_pos_names=list_of_pos_by_name.keys()

    injvals=[]
    allmins=map(lambda a: np.min(a[param].samples), list_of_pos)
    allmaxes=map(lambda a: np.max(a[param].samples), list_of_pos)
    min_pos=np.min(allmins)
    max_pos=np.max(allmaxes)
 
    patch_list=[]
    max_y=1.

    posbins=np.linspace(min_pos,max_pos,50)

    for name,posterior in list_of_pos_by_name.items():
        colour=color_by_name[name]
        #myfig.gca(autoscale_on=True)
        if posterior[param].injval:
            injvals.append(posterior[param].injval)

        try:
            n,bins=np.histogram(posterior[param].samples,bins=posbins,normed=True,new=True)
        except:
            n,bins=np.histogram(posterior[param].samples,bins=posbins,normed=True)

        if min(bins)==max(bins):
            print 'Skipping '+param
            continue

        (n, bins, patches)=plt.hist(posterior[param].samples,bins=bins,histtype='step',label=name,normed=True,hold=True,color=color_by_name[name],cumulative='True')#range=(min_pos,max_pos)

        patch_list.append(patches[0])

    top_cl_intervals_list={}
    pos_names=list_of_pos_by_name.keys()


    for name,posterior in list_of_pos_by_name.items():
        #toppoints,injectionconfidence,reses,injection_area,cl_intervals=bppu.greedy_bin_one_param(posterior,{param:greedyBinSizes[param]},[cl])
        cl_intervals=posterior[param].prob_interval([cl])
        colour=color_by_name[name]
        if cl_intervals[0] is not None and cl_lines_flag:
            try:
                plt.plot([cl_intervals[0][0],cl_intervals[0][0]],[0,max_y],color=colour,linestyle='--')
                plt.plot([cl_intervals[0][1],cl_intervals[0][1]],[0,max_y],color=colour,linestyle='--')
            except:
                print "MAX_Y",max_y,[cl_intervals[0][0],cl_intervals[0][0]],[cl_intervals[0][1],cl_intervals[0][1]]
        top_cl_intervals_list[name]=(cl_intervals[0][0],cl_intervals[0][1])

    if cl_lines_flag:
        pos_names.append(str(int(cl*100))+'%')
        patch_list.append(mpl.lines.Line2D(np.array([0.,1.]),np.array([0.,1.]),linestyle='--',color='black'))

    plt.grid()
    plt.xlim(min_pos,max_pos)
    plt.ylim(0,1)
    if legend:
      oned_legend=plt.figlegend(patch_list,pos_names,'right')
      for text in oned_legend.get_texts():
          text.set_fontsize('small')
    plt.xlabel(bppu.plot_label(param))
    plt.ylabel('Cumulative Probability')
    plt.draw()
    #plt.tight_layout()
    if injvals:
        print "Injection parameter is %f"%(float(injvals[0]))
        injpar=injvals[0]
        #if min(pos_samps)<injpar and max(pos_samps)>injpar:
        plt.plot([injpar,injpar],[0,max_y],'r-.',scalex=False,scaley=False,linewidth=4,label='Injection')
    if analyticCDF is not None:
	plt.plot(posbins,map(analyticCDF,posbins),'r')
    return myfig,top_cl_intervals_list#,rkde


def compare_bayes(outdir,names_and_pos_folders,injection_path,eventnum,username,password,reload_flag,clf,ldg_flag,contour_figsize=(4.5,4.5),contour_dpi=250,contour_figposition=[0.15,0.15,0.5,0.75],fail_on_file_err=True,covarianceMatrices=None,meanVectors=None,Npixels2D=50):

    injection=None

    if injection_path is not None and os.path.exists(injection_path) and eventnum is not None:
        eventnum=int(eventnum)
        import itertools
        injections = SimInspiralUtils.ReadSimInspiralFromFiles([injection_path])
        if eventnum is not None:
            if(len(injections)<eventnum):
                print "Error: You asked for event %d, but %s contains only %d injections" %(eventnum,injection_path,len(injections))
                sys.exit(1)
            else:
                injection=injections[eventnum]
    
    #Create analytic likelihood functions if covariance matrices and mean vectors were given
    analyticLikelihood = None
    if covarianceMatrices and meanVectors:
	analyticLikelihood = bppu.AnalyticLikelihood(covarianceMatrices, meanVectors)
    peparser=bppu.PEOutputParser('common')
    pos_list={}
    tp_list={}
    common_params=None
    working_folder=os.getcwd()
    for name,pos_folder in names_and_pos_folders:
        import urlparse

        pos_folder_url=urlparse.urlparse(pos_folder)
        pfu_scheme,pfu_netloc,pfu_path,pfu_params,pfu_query,pfu_fragment=pos_folder_url

        if 'http' in pfu_scheme:

            """
            Retrieve a file over http(s).
            """
            downloads_folder=os.path.join(os.getcwd(),"downloads")
            pos_folder_parse=urlparse.urlparse(pos_folder)
            pfp_scheme,pfp_netloc,pfp_path,pfp_params,pfp_query,pfp_fragment=pos_folder_parse
            head,tail=os.path.split(pfp_path)
            if tail is 'posplots.html' or tail:
                pos_file_part=head
            else:
                pos_file_part=pfp_path
            pos_file_url=urlparse.urlunsplit((pfp_scheme,pfp_netloc,os.path.join(pos_file_part,'posterior_samples.dat'),'',''))
            print pos_file_url
            pos_file=os.path.join(os.getcwd(),downloads_folder,"%s.dat"%name)

            if not os.path.exists(pos_file):
                reload_flag=True

            if reload_flag:
                if os.path.exists(pos_file):
                    os.remove(pos_file)
                if not os.path.exists(downloads_folder):
                    os.makedirs(downloads_folder)
                open_url_curl(pos_file_url,args=["-o","%s"%pos_file])

        elif pfu_scheme is '' or pfu_scheme is 'file':
            pos_file=os.path.join(pos_folder,'%s.dat'%name)
            # Try looking for posterior_samples.dat if name.dat doesn't exist
            if not os.path.exists(pos_file):
                print '%s does not exist, trying posterior_samples.dat'%(pos_file)
                pos_file=os.path.join(pos_folder,'posterior_samples.dat')
        else:
            print "Unknown scheme for input data url: %s\nFull URL: %s"%(pfu_scheme,str(pos_folder_url))
            exit(0)

        print "Reading posterior samples from %s ..."%pos_file

        try:
            common_output_table_header,common_output_table_raw=peparser.parse(open(pos_file,'r'))
        except:
            print 'Unable to read file '+pos_file
            continue

        test_and_switch_param(common_output_table_header,'distance','dist')
        test_and_switch_param(common_output_table_header,'chirpmass','mchirp')
        test_and_switch_param(common_output_table_header,'mc','mchirp')
        test_and_switch_param(common_output_table_header,'asym_massratio','q')
        test_and_switch_param(common_output_table_header,'massratio', 'eta')
        test_and_switch_param(common_output_table_header,'RA','ra')
        test_and_switch_param(common_output_table_header,'rightascension','ra')
        test_and_switch_param(common_output_table_header,'declination','dec')
        test_and_switch_param(common_output_table_header,'tilt_spin1','tilt1')
        test_and_switch_param(common_output_table_header,'tilt_spin2','tilt2')

        if 'LI_MCMC' in name or 'FU_MCMC' in name:

            try:

                idx=common_output_table_header.index('iota')
                print "Inverting iota!"

                common_output_table_raw[:,idx]= np.pi*np.ones(len(common_output_table_raw[:,0])) - common_output_table_raw[:,idx]

            except:
                pass


        # try:
        #     print "Converting phi_orb-> 2phi_orb"
        #     idx=common_output_table_header.index('phi_orb')
        #     common_output_table_header[idx]='2phi_orb'
        #     common_output_table_raw[:,idx]= 2*common_output_table_raw[:,idx]
        # except:
        #     pass

        try:
            print "Converting iota-> cos(iota)"
            idx=common_output_table_header.index('iota')
            common_output_table_header[idx]='cos(iota)'
            common_output_table_raw[:,idx]=np.cos(common_output_table_raw[:,idx])
        except:
            pass

        #try:
        #    print "Converting tilt1 -> cos(tilt1)"
        #    idx=common_output_table_header.index('tilt1')
        #    common_output_table_header[idx]='cos(tilt1)'
        #    common_output_table_raw[:,idx]=np.cos(common_output_table_raw[:,idx])
        #except:
        #    pass

        #try:
        #    print "Converting tilt2 -> cos(tilt2)"
        #    idx=common_output_table_header.index('tilt2')
        #    common_output_table_header[idx]='cos(tilt2)'
        #    common_output_table_raw[:,idx]=np.cos(common_output_table_raw[:,idx])
        #except:
        #    pass

        try:
            print "Converting thetas -> cos(thetas)"
            idx=common_output_table_header.index('thetas')
            common_output_table_header[idx]='cos(thetas)'
            common_output_table_raw[:,idx]=np.cos(common_output_table_raw[:,idx])
        except:
            pass

        try:
            print "Converting beta -> cos(beta)"
            idx=common_output_table_header.index('beta')
            common_output_table_header[idx]='cos(beta)'
            common_output_table_raw[:,idx]=np.cos(common_output_table_raw[:,idx])
        except:
            pass

        try:
            idx=common_output_table_header.index('f_ref')
            injFrefs=np.unique(common_output_table_raw[:,idx])
            if len(injFrefs) == 1:
              injFref = injFrefs[0]
              print "Using f_ref in results as injected value"
        except:
            injFref = None
            pass

        pos_temp=bppu.Posterior((common_output_table_header,common_output_table_raw),SimInspiralTableEntry=injection, injFref=injFref)

        if 'a1' in pos_temp.names and min(pos_temp['a1'].samples)[0] < 0:
          pos_temp.append_mapping('spin1', lambda a:a, 'a1')
          pos_temp.pop('a1')
          pos_temp.append_mapping('a1', lambda a:np.abs(a), 'spin1')
        if 'a2' in pos_temp.names and min(pos_temp['a2'].samples)[0] < 0:
          pos_temp.append_mapping('spin2', lambda a:a, 'a2')
          pos_temp.pop('a2')
          pos_temp.append_mapping('a2', lambda a:np.abs(a), 'spin2')


        if 'm1' in pos_temp.names and 'm2' in pos_temp.names:
          print "Calculating total mass"
          pos_temp.append_mapping('mtotal', lambda m1,m2: m1+m2, ['m1','m2'])
        if 'mass1' in pos_temp.names and 'mass2' in pos_temp.names:
          print "Calculating total mass"
          pos_temp.append_mapping('mtotal', lambda m1,m2: m1+m2, ['mass1','mass2'])

        try:
            idx=common_output_table_header.index('m1')

            idx2=common_output_table_header.index('m2')

            if pos_temp['m1'].mean<pos_temp['m2'].mean:
                print "SWAPPING MASS PARAMS!"
                common_output_table_header[idx]='x'
                common_output_table_header[idx2]='m1'
                common_output_table_header[idx]='m2'
                pos_temp=bppu.Posterior((common_output_table_header,common_output_table_raw),SimInspiralTableEntry=injection)
        except:
            pass

        pos_list[name]=pos_temp

        if common_params is None:
            common_params=pos_temp.names
        else:
            set_of_pars = set(pos_temp.names)
            common_params=list(set_of_pars.intersection(common_params))

    print "Common parameters are %s"%str(common_params)

    if injection is None and injection_path is not None:
        import itertools
        injections = SimInspiralUtils.ReadSimInspiralFromFiles([injection_path])
        injection=bppu.get_inj_by_time(injections,pos_temp.means['time'])
    if injection is not None:
        for pos in pos_list.values():
            pos.set_injection(injection)

    set_of_pars = set(allowed_params)
    common_params=list(set_of_pars.intersection(common_params))

    print "Using parameters %s"%str(common_params)

    if not os.path.exists(os.path.join(os.getcwd(),'results')):
        os.makedirs('results')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    pdfdir=os.path.join(outdir,'pdfs')
    if not os.path.exists(pdfdir):
        os.makedirs(pdfdir)

    greedy2savepaths=[]

    if common_params is not [] and common_params is not None: #If there are common parameters....
        colorlst=bppu.__default_color_lst

        if len(common_params)>1: #If there is more than one parameter...
            temp=copy.copy(common_params)
            #Plot two param contour plots

            #Assign some colours to each different analysis result
            color_by_name={}
            hatches_by_name={}
            my_cm=mpl_cm.Dark2
            cmap_size=my_cm.N
            color_idx=0
            color_idx_max=len(names_and_pos_folders)
            cmap_array=my_cm(np.array(range(cmap_size)))
            #cmap_array=['r','g','b','c','m','k','0.5','#ffff00']
            hatches=['/','\\','|','-','+','x','o','O','.','*']
            ldg='auto'
            if not ldg_flag:
              ldg=None
            for name,infolder in names_and_pos_folders:
                #color_by_name=cmap_array[color_idx]
                color_by_name[name]=cmap_array[int(floor(color_idx*cmap_size/color_idx_max)),:]
                color_idx=(color_idx+1) % color_idx_max
                hatches_by_name[name]=hatches[color_idx]

            for i,j in all_pairs(temp):#Iterate over all unique pairs in the set of common parameters
                pplst=[i,j]
                rpplst=pplst[:]
                rpplst.reverse()

                pplst_cond=(pplst in twoDplots)
                rpplst_cond=(rpplst in twoDplots)
                if pplst_cond or rpplst_cond:#If this pair of parameters is in the plotting list...

                    try:
                        print '2d plots: building ',i,j
                        greedy2Params={i:greedyBinSizes[i],j:greedyBinSizes[j]}
                    except KeyError:
                        continue

                    name_list=[]
                    cs_list=[]

                    slinestyles=['solid', 'dashed', 'dashdot', 'dotted']

                    fig=bppu.plot_two_param_kde_greedy_levels(pos_list,greedy2Params,TwoDconfidenceLevels,color_by_name,figsize=contour_figsize,dpi=contour_dpi,figposition=contour_figposition,legend=ldg,line_styles=slinestyles,hatches_by_name=hatches_by_name,Npixels=Npixels2D)
                    if fig is None: continue
                    #fig=bppu.plot_two_param_greedy_bins_contour(pos_list,greedy2Params,TwoDconfidenceLevels,color_by_name,figsize=contour_figsize,dpi=contour_dpi,figposition=contour_figposition)
                    greedy2savepaths.append('%s-%s.png'%(pplst[0],pplst[1]))
                    fig.savefig(os.path.join(outdir,'%s-%s.png'%(pplst[0],pplst[1])),bbox_inches='tight')
                    fig.savefig(os.path.join(pdfdir,'%s-%s.pdf'%(pplst[0],pplst[1])),bbox_inches='tight')


            plt.clf()
        oned_data={}
        #confidence_levels={}
        confidence_levels=[{},{},{},{}]
        confidence_uncertainty={}
        for param in common_params:
            print "Plotting comparison for '%s'"%param

            cl_table_header='<table><th>Run</th>'
            cl_table={}
            save_paths=[]
            cl_table_min_max_str='<tr><td> Min | Max </td>'
            level_index=0
            for confidence_level in OneDconfidenceLevels:
		if analyticLikelihood:
		  pdf=analyticLikelihood.pdf(param)
		  cdf=analyticLikelihood.cdf(param)
		else:
		  pdf=None
		  cdf=None
		  
                cl_table_header+='<th colspan="2">%i%% (Lower|Upper)</th>'%(int(100*confidence_level))
                hist_fig,cl_intervals=compare_plots_one_param_line_hist(pos_list,param,confidence_level,color_by_name,cl_lines_flag=clf,legend=ldg,analyticPDF=pdf)
                hist_fig2,cl_intervals=compare_plots_one_param_line_hist_cum(pos_list,param,confidence_level,color_by_name,cl_lines_flag=clf,analyticCDF=cdf,legend=ldg)

                # Save confidence levels and uncertainty
                #confidence_levels[param]=[]
                confidence_levels[level_index][param]=[]
                
                for name,pos in pos_list.items():
                    median=pos[param].median
                    low,high=cl_intervals[name]
                    #confidence_levels[param].append((name,low,median,high))
                    confidence_levels[level_index][param].append((name,low,median,high))
                    
                level_index=level_index+1
                cl_bounds=[]
                poses=[]
                for name,pos in pos_list.items():
                    cl_bounds.append(cl_intervals[name])
                    poses.append(pos[param])
                confidence_uncertainty[param]=bppu.confidence_interval_uncertainty(confidence_level, cl_bounds, poses)

                save_path=''
                if hist_fig is not None:
                    save_path=os.path.join(outdir,'%s_%i.png'%(param,int(100*confidence_level)))
                    save_path_pdf=os.path.join(pdfdir,'%s_%i.pdf'%(param,int(100*confidence_level)))
                    try:
                      plt.tight_layout(hist_fig)
                      plt.tight_layout(hist_fig2)
                    except:
                      pass
                    hist_fig.savefig(save_path,bbox_inches='tight')
                    hist_fig.savefig(save_path_pdf,bbox_inches='tight')
                    save_paths.append(save_path)
                    save_path=os.path.join(outdir,'%s_%i_cum.png'%(param,int(100*confidence_level)))
                    save_path_pdf=os.path.join(pdfdir,'%s_%i_cum.pdf'%(param,int(100*confidence_level)))
                    hist_fig2.savefig(save_path,bbox_inches='tight')
                    hist_fig2.savefig(save_path_pdf,bbox_inches='tight')
                    save_paths.append(save_path)
                min_low,max_high=cl_intervals.values()[0]
                for name,interval in cl_intervals.items():
                    low,high=interval
                    if low<min_low:
                        min_low=low
                    if high>max_high:
                        max_high=high
                    try:
                        cl_table[name]+='<td>%s</td><td>%s</td>'%(low,high)
                    except:
                        cl_table[name]='<td>%s</td><td>%s</td>'%(low,high)
                cl_table_min_max_str+='<td>%s</td><td>%s</td>'%(min_low,max_high)
            cl_table_str=cl_table_header
            for name,row_contents in cl_table.items():
                cl_table_str+='<tr><td>%s<font color="%s"></font></td>'%(name,str(mpl_colors.rgb2hex(color_by_name[name][0:3])))#,'&#183;'.encode('utf-8'))

                cl_table_str+=row_contents+'</tr>'
            cl_table_str+=cl_table_min_max_str+'</tr>'
            cl_table_str+='</table>'

            cl_uncer_str='<table> <th>Confidence Relative Uncertainty</th> <th>Confidence Fractional Uncertainty</th> <th>Confidence Percentile Uncertainty</th>\n'
            cl_uncer_str+='<tr> <td> %g </td> <td> %g </td> <td> %g </td> </tr> </table>'%(confidence_uncertainty[param][0], confidence_uncertainty[param][1], confidence_uncertainty[param][2])

            ks_matrix=compute_ks_pvalue_matrix(pos_list, param)

            N=ks_matrix.shape[0]+1

            # Make up KS-test table
            ks_table_str='<table><th colspan="%d"> K-S test p-value matrix </th>'%N

            # Column headers
            ks_table_str+='<tr> <td> -- </td> '
            for name,pos in pos_list.items():
                ks_table_str+='<td> %s </td>'%name
            ks_table_str+='</tr>'

            # Now plot rows of matrix
            for i in range(len(pos_list)):
                ks_table_str+='<tr> <td> %s </td>'%(pos_list.keys()[i])
                for j in range(len(pos_list)):
                    if i == j:
                        ks_table_str+='<td> -- </td>'
                    elif ks_matrix[i,j] < 0.05:
                        # Failing at suspiciously low p-value
                        ks_table_str+='<td> <b> %g </b> </td>'%ks_matrix[i,j]
                    else:
                        ks_table_str+='<td> %g </td>'%ks_matrix[i,j]

                ks_table_str+='</tr>'

            ks_table_str+='</table>'

            oned_data[param]=(save_paths,cl_table_str,ks_table_str,cl_uncer_str)
            
    # Watch out---using private variable _logL
    max_logls = [[name,max(pos._logL)] for name,pos in pos_list.items()]
    dics = [pos.DIC for name, pos in pos_list.items()]

    return greedy2savepaths,oned_data,confidence_uncertainty,confidence_levels,max_logls,dics

def output_confidence_levels_tex(clevels,outpath):
    """Outputs a LaTeX table of parameter and run medians and confidence levels."""
    outfile=open(os.path.join(outpath,'confidence_table.tex'), 'w')
    for level_index in range(len(OneDconfidenceLevels)):
        params=clevels[level_index].keys()

        clevels_by_name={}
        for param in clTableParams:
            if param in params:
                for name,low,med,high in clevels[level_index][param]:
                    if name in clevels_by_name:
                        clevels_by_name[name].append((param,low,med,high))
                    else:
                        clevels_by_name[name] = [(param,low,med,high)]

        try:
            outfile.write('confidence level %1.3g\n'%OneDconfidenceLevels[level_index])
            outfile.write(r'\begin{tabular}{|l||')
            for param in clTableParams:
                if param in params:
                    outfile.write('c|')
            outfile.write('}\n')

            outfile.write(r'\hline ')
            for param in clTableParams:
                if param in params:
                    tparam=paramNameLatexMap.get(param,param)
                    outfile.write(r'& $%s$ '%tparam)
            outfile.write('\\\\ \n \\hline \\hline ')

            for name,levels in clevels_by_name.items():
                outfile.write(name)
                for param,low,med,high in levels:
                    outfile.write(r' & $%0.5g^{%0.5g}_{%0.5g}$ '%(med,high,low))
                outfile.write('\\\\ \n')

            outfile.write('\\hline \n \\end{tabular}')
        finally:
            outfile.write('\n\n')

    outfile.close()

def output_confidence_levels_dat(clevels,outpath):
    """Outputs a LaTeX table of parameter and run medians and confidence levels."""
    outfile=open(os.path.join(outpath,'confidence_table.dat'), 'w')
    for level_index in range(len(OneDconfidenceLevels)):
        params=clevels[level_index].keys()

        clevels_by_name={}
        for param in clTableParams:
            if param in params:
                for name,low,med,high in clevels[level_index][param]:
                    if name in clevels_by_name:
                        clevels_by_name[name].append((param,low,med,high))
                    else:
                        clevels_by_name[name] = [(param,low,med,high)]

        try:
            outfile.write('%1.3g\t'%OneDconfidenceLevels[level_index])
            for param in clTableParams:
                if param in params:
                    tparam=paramNameLatexMap.get(param,param)
                    outfile.write('%s\t'%param)
            outfile.write('\n') 

            for name,levels in clevels_by_name.items():
                outfile.write(name)
                for param,low,med,high in levels:
                    outfile.write('\t%6.6g - %6.6g'%(low,high))
                outfile.write('\n')
        finally:
            outfile.write('\n')

    outfile.close()

def output_confidence_uncertainty(cluncertainty, outpath):
    outfile=open(os.path.join(outpath, 'confidence_uncertainty.dat'), 'w')
    try:
        params=cluncertainty.keys()
        uncer=cluncertainty.values()

        outfile.write('# Uncertainty in confidence levels.\n')
        outfile.write('# First row is relative uncertainty (wrt to parameter mean).\n')
        outfile.write('# Second row is fractional uncertainty (wrt to combined conf interval).\n')
        outfile.write('# Third row is percentile uncertainty (wrt combined samples).\n')
        outfile.write('# ')
        for param in params:
            outfile.write(str(param) + ' ')
        outfile.write('\n')

        rel = np.array([d[0] for d in uncer])
        fracs = np.array([d[1] for d in uncer])
        quants = np.array([d[2] for d in uncer])

        np.savetxt(outfile, np.reshape(rel, (1, -1)))
        np.savetxt(outfile, np.reshape(fracs, (1, -1)))
        np.savetxt(outfile, np.reshape(quants, (1,-1)))
    finally:
        outfile.close()
            
if __name__ == '__main__':
    from optparse import OptionParser
    parser=OptionParser()
    parser.add_option("-o","--outpath", dest="outpath",help="Make page and plots in DIR.", metavar="DIR")
    parser.add_option("-p","--pos",dest="pos_list",action="append",help="Path to folders containing output of cbcBayesPostProc.")
    parser.add_option("-n","--name",dest="names",action="append",help="Name of posterior result e.g. followupMCMC 2.5PN (optional)")
    parser.add_option("-i","--inj",dest="inj",help="Path of injection XML containing SimInspiralTable (optional).")
    parser.add_option("-e","--eventnum",dest="eventnum",help="Sim ID of injection described in injection XML (optional).")
    parser.add_option("-u",dest="username",help="User name for https authenticated content (optional).")
    parser.add_option("-x",dest="password",help="Password for https authenticated content (optional).")
    parser.add_option("--reload",dest="reload_flag",action="store_true",help="Re-download all pos files (optional).")
    parser.add_option("--hide-cl-lines",dest="clf",action="store_false",default=True,help="Hide confidence level lines on 1D plots for clarity (optional).")
    parser.add_option("--contour-dpi",dest="cdpi",default=250,help="DPI for contour plot (optional).")
    parser.add_option("--contour-width",dest="cw",default=7,help="Width (in inches) of contour plots (optional).")
    parser.add_option("--contour-height",dest="ch",default=6,help="Height (in inches) of contour plots (optional).")
    parser.add_option("--contour-plot-width",dest="cpw",default=0.5,help="Relative width of plot element 0.15<width<1 (optional).")
    parser.add_option("--contour-plot-height",dest="cph",default=0.76,help="Relative height of plot element 0.15<width<1 (optional).")
    parser.add_option("--no-legend",dest="ldg_flag",action="store_false",default=True,help="Hide legend (optional).")
    parser.add_option("--ignore-missing-files",dest="readFileErr",default=False,action="store_true",help="Do not fail when files are missing (optional).")
    parser.add_option("-c","--covarianceMatrix",dest="covarianceMatrices",action="append",default=None,help="CSV file containing covariance (must give accompanying mean vector CSV. Can add more than one matrix.")
    parser.add_option("-m","--meanVectors",dest="meanVectors",action="append",default=None,help="Comma separated list of locations of the multivariate gaussian described by the correlation matrix.  First line must be list of params in the order used for the covariance matrix.  Provide one list per covariance matrix.")
    parser.add_option("--no2D",dest="no2d",action="store_true",default=False,help="Disable 2D plots")
    parser.add_option("--npixels-2d",dest="npixels_2d",action="store",type="int",default=50,help="Number of pixels on a side of the 2D plots (default 50)",metavar="N")
    
    (opts,args)=parser.parse_args()

    if opts.outpath is None:
        print "No output directory specified. Output will be saved to PWD : %s"%os.getcwd()
        outpath=os.getcwd()
    else:
        outpath=opts.outpath

    if opts.pos_list is None:
        print "No input paths given!"
        exit(1)

    if opts.names is None:
        print "No names given, making some up!"
        names=[]
        for i in range(len(opts.pos_list)):
            names.append(str(i))
    else:
        names=opts.names

    if len(opts.pos_list)!=len(names):
        print "Either add names for all posteriors or dont put any at all!"

    # Sort inputs alphabetically
    names,pos_list = zip(*sorted(zip(names,opts.pos_list)))
    
    if opts.no2d:
        twoDplots=[]


    greedy2savepaths,oned_data,confidence_uncertainty,confidence_levels,max_logls,dics=compare_bayes(outpath,zip(names,pos_list),opts.inj,opts.eventnum,opts.username,opts.password,opts.reload_flag,opts.clf,opts.ldg_flag,contour_figsize=(float(opts.cw),float(opts.ch)),contour_dpi=int(opts.cdpi),contour_figposition=[0.15,0.15,float(opts.cpw),float(opts.cph)],fail_on_file_err=not opts.readFileErr,covarianceMatrices=opts.covarianceMatrices,meanVectors=opts.meanVectors,Npixels2D=int(opts.npixels_2d))

    ####Print Confidence Levels######
    output_confidence_levels_tex(confidence_levels,outpath)    
    output_confidence_levels_dat(confidence_levels,outpath)
    
    ####Save confidence uncertainty#####
    output_confidence_uncertainty(confidence_uncertainty,outpath)

    ####Print HTML!#######

    compare_page=bppu.htmlPage('Compare PDFs (single event)',css=bppu.__default_css_string)

    param_section=compare_page.add_section('Meta')

    param_section_write='<div><p>This comparison was created from the following analyses</p>'
    param_section_write+='<table border="1">'
    param_section_write+='<th>Analysis</th> <th> max(log(L)) </th> <th> DIC </th>'
    for (name,logl_max), dic in zip(max_logls, dics):
        param_section_write+='<tr><td><a href="%s">%s</a></td> <td>%g</td> <td>%.1f</td></tr>'%(dict(zip(names,pos_list))[name],name,logl_max,dic)
    param_section_write+='</table></div>'

    param_section.write(param_section_write)
    param_section.write('<div><p><a href="confidence_table.tex">LaTeX table</a> of medians and confidence levels.</p></div>')
    if oned_data:

        param_section=compare_page.add_section('1D marginal posteriors')

        for param_name,data in oned_data.items():
            param_section.h3(param_name)
            save_paths,cl_table_str,ks_table_str,cl_uncer_str=data
            clf_toggle=False
            for save_path in save_paths:
                head,plotfile=os.path.split(save_path)
                param_section.write('<img src="%s"/>'%str(plotfile))

            param_section.write(cl_table_str)
            param_section.write(cl_uncer_str)
            param_section.write(ks_table_str)

    if greedy2savepaths:

        param_section=compare_page.add_section('2D greedy bin histograms')
        for plot_path in greedy2savepaths:
            temp,param_name=os.path.split(plot_path)
            param_name=param_name.split('.')[0]
            head,plotfile=os.path.split(plot_path)
            param_section.write('<img src="%s"/>'%str(plotfile))#str(os.path.relpath(plot_path,outpath)))



    compare_page_footer=compare_page.add_section('')
    compare_page_footer.p('Produced using cbcBayesCompPos.py at '+strftime("%Y-%m-%d %H:%M:%S")+' .')

    cc_args=copy.deepcopy(sys.argv)

    #remove username and password
    try:
        user_idx=cc_args.index('-u')
        cc_args[user_idx+1]='<LIGO username>'
    except:
        pass

    try:
        pass_idx=cc_args.index('-x')
        cc_args[pass_idx+1]='<LIGO password>'
    except:
        pass

    cc_args_str=''
    for cc_arg in cc_args:
        cc_args_str+=cc_arg+' '

    compare_page_footer.p('Command line: %s'%cc_args_str)
    compare_page_footer.p(git_version.verbose_msg)


    resultspage=open(os.path.join(outpath,'index.html'),'w')
    resultspage.write(str(compare_page))
    resultspage.close()
