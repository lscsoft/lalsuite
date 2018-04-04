#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcBayesPostProc.py
#       Copyright 2010 Benjamin Aylott <benjamin.aylott@ligo.org>, John Veitch <john.veitch@ligo.org>,
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

import sys, os
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from pylal import SimInspiralUtils
from pylal import bayespputils as bppu
from pylal import git_version

from time import strftime

__author__="John Veitch <john.veitch@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

def addm1m2(pos):
    injection = pos.injection()
    # Produce m1,m2 if requested and not available
    if ('mc' in pos.names or 'mchirp' in pos.names) and \
        'eta' in pos.names and \
        ('mass1' not in pos.names or 'm1' not in pos.names) and\
        ('m2' not in pos.names or 'm2' not in pos.names):

        if 'mc' in pos.names:
            mchirp_name='mc'
        else:
            mchirp_name='mchirp'
        inj_mass1=None
        inj_mass2=None
        if injection:
            inj_mass1,inj_mass2=bppu.mc2ms(injection.mchirp,injection.eta)
        mass1_samps,mass2_samps=bppu.mc2ms(pos[mchirp_name].samples,pos['eta'].samples)
        mass1_pos=bppu.PosteriorOneDPDF('m1',mass1_samps,injected_value=inj_mass1)
        mass2_pos=bppu.PosteriorOneDPDF('m2',mass2_samps,injected_value=inj_mass2)
        pos.append(mass1_pos)
        pos.append(mass2_pos)
    return


def getInjByTime(time,injections):
    import itertools
    injections = itertools.ifilter(lambda a: abs(float(a.get_end()) - time) < 0.1, injections)
    if len(injections)!=1:
        print 'ERROR: Found more than one injection with end time %f\n Please specify input files in order and right number'%(time)
        os.exit(1)
    return injections.next()

def makeInjObjs(injfile,posfiles):
    """
    Make a list of results objects from the posfiles and injections
    """
    getByTime=False
    resObjs=[]
    if injfile:
        import itertools
        injections = SimInspiralUtils.ReadSimInspiralFromFiles([injfile])
        if(len(injections)!=len(posfiles)):
            print 'Different numbers of injections and posterior files, attempting to recognise by time'
            getByTime=True
    i=0
    for posfile in posfiles:
        peparser=bppu.PEOutputParser('common')
        resObj=peparser.parse(open(posfile,'r'))
        pos=bppu.Posterior(resObj)
        if not getByTime:
            injection=injections[i]
        else:
            time=pos['time'].mean
            injection=getInjByTime(time,injections)
        if injection is None:
            continue
        pos.set_injection(injection)
        resObjs.append(pos)
        i=i+1
    return resObjs

def makeOutputPage(objs, params, outdir, confidencelevels):
    """
    Make a summary page with table of results and plots for each param in params
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    #Bin size/resolution for binning. Need to match (converted) column names.
    GreedyRes={'mc':0.025,'m1':0.1,'m2':0.1,'mass1':0.1,'mass2':0.1,'mtotal':0.1,'eta':0.001,'iota':0.01,'time':1e-3,'distance':1.0,'dist':1.0,'mchirp':0.025,'a1':0.02,'a2':0.02,'phi1':0.05,'phi2':0.05,'theta1':0.05,'theta2':0.05,'ra':0.05,'dec':0.05}


    html=bppu.htmlPage('Injection Summary',css=bppu.__default_css_string)
    html_meta=html.add_section('Summary')
    html_meta.p('Analysed %i injections.'%(len(objs)))
    
    # Make a directory for stdacc and errorbar plots
    accdir=os.path.join(outdir,'stdacc')
    if not os.path.isdir(accdir):
        os.makedirs(accdir)
    errdir=os.path.join(outdir,'errbar')
    if not os.path.isdir(errdir):
        os.makedirs(errdir)
    boxdir=os.path.join(outdir,'boxplot')
    if not os.path.isdir(boxdir):
        os.makedirs(boxdir)
    # Calculate confidence intervals for each injection and parameter
    
    #for par in params:
    #    par=par.lower()
    #    print 'Binning %s to determine confidence intervals'%(par)
    #    
    #    for obj in objs:
    #        try:
    #            obj[par]
    #        except KeyError:
    #            print 'No input chain for %s, skipping'%(par)
    #            continue
    #        try:
    #            GreedyRes[par]
    #        except KeyError:
    #            print 'No bin size for %s, skipping'%(par)
    #            continue
    #        binParams={par: GreedyRes[par]}

     #       toppoints,injectionconfidence,reses,injection_area=bppu.greedy_bin_one_param(obj,binParams, confidencelevels)
     #       oneDContCL, oneDContInj = bppu.contigious_interval_one_param(obj, binParams, confidencelevels)
    
    print 'Calculating std accuracies'
    # Calculate Std Accuracies
    stdacc={}
    mean={}
    std={}
    for par in params:
        if not reduce( lambda a,b: a and b, map(lambda z: par in z.names, objs)):
            print '%s not found in all objects, skipping'%(par)
            continue
        stdacc[par]=[]
        mean[par]=[]
        std[par]=[]
        for obj in objs:
            oneDpos=obj[par]
            stdacc[par].append(oneDpos.stacc)
            mean[par].append(oneDpos.mean)
            std[par].append(oneDpos.stdev)
    
    html_scatter=html.add_section('Standard Accuracy')
    html_scatter_content='<table>'
    print 'Making std accuracy plots'
    # Make a scatter plot for std accuracy against each injection parameter
    for p1 in params: # X axis, injection values
        try:
            injval=map(lambda o: o._getinjpar(p1), objs)
        except KeyError:
            print 'Error! Unable to find parameter %s in injection!'%(p1)
            continue
        html_scatter_content+='<tr><td>%s</td>'%(p1)
        for p2 in params: # Y axis
            try:
                stdacc[p2]
            except KeyError:
                print 'No stdacc data for %s, skipping'%(p2)
                html_scatter_content+='<td>%s not found</td>'%(p2)
                continue
            fig = scatterAcc(injval,stdacc[p2],p1,p2)
            figname=p1+'-'+p2+'.png'
            figpath=os.path.join(accdir,figname)
            fig.savefig(figpath)
            html_scatter_content+='<td><img src="%s" alt="%s-%s" /></td>'%(figpath,p1,p2)
        html_scatter_content+='</tr>'
    html_scatter_content+='</table>'
    html_scatter.write(html_scatter_content)
    print 'Making errorbar plots'
    # Make an errorbar plot for each parameter against each injection parameter
    html_err=html.add_section('Parameter Estimates')
    html_err_content='<table>'
    for p1 in params: # X axis, injection values
        try:
            injval=map(lambda o: o._getinjpar(p1),objs)
        except KeyError:
            print 'Error! Unable to find parameter %s in injection!'%(p1)
            continue
        html_err_content+='<tr><td>%s</td>'%(p1)
        for p2 in params:
            try:
                mean[p2]
            except KeyError:
                print 'No mean for %s, skipping'%(p2)
                html_err_content+='<td>%s not found</td>'%(p2)
                continue
            yinjval=map(lambda o: o._getinjpar(p2), objs)
            fig = plt.figure()
            plt.errorbar(injval,mean[p2],std[p2],linestyle='None')
            plt.plot(injval,yinjval,'gx')
            plt.xlabel(p1)
            plt.ylabel(p2)
            figname=p1+'-'+p2+'.png'
            figpath=os.path.join(errdir,figname)
            fig.savefig(figpath)
            html_err_content+='<td><img src="%s" alt="%s-%s" /></td>'%(figpath,p1,p2)
        html_err_content+='</tr>'
    html_err_content+='</table>'
    html_err.write(html_err_content)

    # Box and whiskers plot for each parameter pair
    html_box=html.add_section('Box plots')
    html_box_content='<table>'
    print 'Making box plots'
    for p1 in params:
        # Error checking to be put here
        injval=map(lambda o: o._getinjpar(p1),objs)
        html_box_content+='<tr><td>%s</td>'%(p1)
        for p2 in params:
            try:
                mean[p2]
            except KeyError:
                print 'No mean for %s, skipping'%(p2)
                html_box_content+='<td>%s not found</td>'%(p2)
                continue
            posteriors=map(lambda o: o[p2].samples, objs)
            yinjval=map(lambda o: o._getinjpar(p2),objs)
            fig=plt.figure()
            upper=max(injval)
            lower=min(injval)
            boxwidth=0.75*(upper-lower)/len(injval)
            plt.boxplot(posteriors,positions=injval,widths=boxwidth )
            plt.plot(injval,yinjval,'gx')
            plt.xlabel(p1)
            plt.ylabel(p2)
            plt.xlim(lower-0.5*boxwidth,upper+0.5*boxwidth)
            figname=p1+'-'+p2+'.png'
            figpath=os.path.join(boxdir,figname)
            fig.savefig(figpath)
            html_box_content+='<td><img src="%s" alt="%s-%s" /></td>'%(figpath,p1,p2)
        html_box_content+='</tr>'
    html_box_content+='</table>'
    html_box.write(html_box_content)

    html_footer=html.add_section('')
    html_footer.p('Produced using cbcBayesInjProc.py at '+strftime("%Y-%m-%d %H:%M:%S")+' .')
    html_footer.p(git_version.verbose_msg)

    # Save page
    resultspage=open(os.path.join(outdir,'injectionsummary.html'),'w')
    resultspage.write(str(html))
    resultspage.close()

def scatterAcc(x,y,xname,yname):
    print 'Plotting %s - %s scatter plot'%(xname,yname)
    myfig = plt.figure()
    plt.scatter(x,y,label='Std Acc of %s'%(yname))
    plt.xlabel(xname)
    plt.ylabel(yname)
    return myfig

if __name__=='__main__':
    from optparse import OptionParser
    usage ="""%prog -i injections.xml -o outdir -p param1 [-p param2 ... ] posterior_1.dat [posterior_2.dat ... ]
    Parse a list of hardware injections and compile summary statistics
    page for all parameter specified with the -p option. Output is stored in outdir.
    
    If you specify the same number of posterior.dat files as there are injections,
    the code will assume the same ordering as in the injection.dat file and associate each
    posterior with that injection. Otherwise, it will attempt to match posteriors to injections
    by the end time of the injection.
    """
    parser=OptionParser(usage)
    parser.add_option("-i","--inj",help="Injection XML",metavar="INJ.XML",default=None,type="string")
    parser.add_option("-p","--parameter",help="Name of parameter to analyse, can be specified multiple times",metavar="mchirp",action="append",default=[],type="string")
    parser.add_option("-o","--outdir",help="output directory",metavar="PATH",type="string",default="./")

    (opts,args)=parser.parse_args()
    if not os.path.isdir(opts.outdir):
        os.makedirs(opts.outdir)

    
    confidencelevels=[0.65,0.9,0.95]
    results_objects=makeInjObjs(opts.inj,args)
    makeOutputPage(results_objects,opts.parameter,opts.outdir,confidencelevels)

