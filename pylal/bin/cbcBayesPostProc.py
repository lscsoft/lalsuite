#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#       cbcBayesPostProc.py
#
#       Copyright 2010
#       Benjamin Aylott <benjamin.aylott@ligo.org>,
#       Benjamin Farr <bfarr@u.northwestern.edu>,
#       Will M. Farr <will.farr@ligo.org>,
#       John Veitch <john.veitch@ligo.org>
#       Vivien Raymond <vivien.raymond@ligo.org>
#       Salvatore Vitale <salvatore.vitale@ligo.org>
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
import numpy as np
from numpy import array,exp,cos,sin,arcsin,arccos,sqrt,size,mean,column_stack,cov,unique,hsplit,correlate,log,dot,power,squeeze,sort
from scipy import stats

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cbcBayesPlotSpinDisk as cbcdiskspin

# Default font properties
fig_width_pt = 246  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
matplotlib.rcParams.update(
        {'axes.labelsize': 16,
        'font.size':       16,
        'legend.fontsize': 16,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'text.usetex': False,
        'figure.figsize': fig_size,
        'font.family': "serif",
        'font.serif': ['Times','Palatino','New Century Schoolbook','Bookman','Computer Modern Roman','Times New Roman','Liberation Serif'],
        'font.weight':'normal',
        'font.size':16,
        'savefig.dpi': 120
        })

#local application/library specific imports
from pylal import SimInspiralUtils
from pylal import bayespputils as bppu
from pylal import git_version

from glue.ligolw import table
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils

__author__="Ben Aylott <benjamin.aylott@ligo.org>, Ben Farr <bfarr@u.northwestern.edu>, Will M. Farr <will.farr@ligo.org>, John Veitch <john.veitch@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

def email_notify(address,path):
    import smtplib
    import subprocess
    import socket
    import os
    address=address.split(',')
    SERVER="localhost"
    USER=os.environ['USER']
    HOST=socket.getfqdn()
    FROM=USER+'@'+HOST
    SUBJECT="LALInference result is ready at "+HOST+"!"
    # Guess the web space path for the clusters
    fslocation=os.path.abspath(path)
    webpath='posplots.html'
    if 'public_html' in fslocation:
        k='public_html/'
    elif 'WWW' in fslocation:
        k='WWW/'
    elif 'www_html' in fslocation:
        k='www_html/'
    else:
        k=None
    if k is not None:
        (a,b)=fslocation.split(k)
        webpath=os.path.join('~%s'%(USER),b,webpath)
        onweb=True
    else:
        (c,d)=fslocation.split(USER)
        for k in ['public_html/','WWW/','www_html/']:
            trypath=os.path.normpath(c+USER+'/'+k+d)
            #Follow symlinks
            if os.path.realpath(trypath)==os.path.normpath(fslocation):
                (a,b)=trypath.split(k)
                print USER
                print b
                print webpath
                webpath=os.path.join('~%s'%(USER),b,webpath)
                print webpath
                onweb=True
                break
            else:
                webpath=os.path.join(fslocation,'posplots.html')
                onweb=False
    if 'atlas' in HOST:
        url="https://atlas1.atlas.aei.uni-hannover.de/"
    elif 'cit' in HOST or 'caltech' in HOST:
        url="https://ldas-jobs.ligo.caltech.edu/"
    elif 'ligo-wa' in HOST:
        url="https://ldas-jobs.ligo-wa.caltech.edu/"
    elif 'ligo-la' in HOST:
        url="https://ldas-jobs.ligo-la.caltech.edu/"
    elif 'uwm' in HOST or 'nemo' in HOST:
        url="https://ldas-jobs.phys.uwm.edu/"
    elif 'phy.syr.edu' in HOST:
        url="https://sugar-jobs.phy.syr.edu/"
    elif 'arcca.cf.ac.uk' in HOST:
        url="https://geo2.arcca.cf.ac.uk/"
    elif 'vulcan' in HOST:
        url="https://galahad.aei.mpg.de/"
    else:
        if onweb:
          url="http://%s/"%(HOST)
        else:
          url=HOST+':'
    url=url+webpath

    TEXT="Hi "+USER+",\nYou have a new parameter estimation result on "+HOST+".\nYou can view the result at "+url+"\n"
    cmd='echo "%s" | mail -s "%s" "%s"'%(TEXT,SUBJECT,', '.join(address))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.STDOUT, shell=True)
    (out, err) = proc.communicate()
    #print "program output %s error %s:"%(out,err)

#Import content handler
from pylal.SimInspiralUtils import ExtractSimInspiralTableLIGOLWContentHandler
lsctables.use_in(ExtractSimInspiralTableLIGOLWContentHandler)


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

    def floatable(str):
      try:
        float(str)
        return True
      except ValueError:
        return False

    for arg in parser.rargs:
      # stop on --foo like options
      if arg[:2] == "--" and len(arg) > 2:
        break
      # stop on -a, but not on -3 or -3.0
      if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
        break
      args.append(arg)

    del parser.rargs[:len(args)]
    #Append new files to list if some already specified
    if getattr(parser.values, opt.dest):
        oldargs = getattr(parser.values, opt.dest)
        oldargs.extend(args)
        args = oldargs
    setattr(parser.values, opt.dest, args)

def dict2html(d,parent=None):
    if not d: return ""
    out=bppu.htmlChunk('div',parent=parent)
    tab=out.tab()
    row=out.insert_row(tab)
    for key in d.keys():
        out.insert_td(row,str(key))
    row2=out.insert_row(tab)
    for val in d.values():
        out.insert_td(row2,str(val))
    return out

def extract_hdf5_metadata(h5grp,parent=None):
    import h5py
    #out=bppu.htmlChunk('div',parent=parent)
    sec=bppu.htmlSection(h5grp.name,htmlElement=parent)
    dict2html(h5grp.attrs,parent=sec)
    for group in h5grp:
        if(isinstance(h5grp[group],h5py.Group)):
            extract_hdf5_metadata(h5grp[group],sec)
    return h5grp


def cbcBayesPostProc(
                        outdir,data,oneDMenus,twoDGreedyMenu,GreedyRes,
                        confidence_levels,twoDplots,
                        #misc. optional
                        injfile=None,eventnum=None,
                        trigfile=None,trignum=None,
                        skyres=None,
                        #direct integration evidence
                        dievidence=False,boxing=64,difactor=1.0,
                        #elliptical evidence
                        ellevidence=False,
                        #manual input of bayes factors optional.
                        bayesfactornoise=None,bayesfactorcoherent=None,
                        #manual input for SNR in the IFOs, optional.
                        snrfactor=None,
                        #nested sampling options
                        ns_flag=False,ns_Nlive=None,
                        #spinspiral/mcmc options
                        ss_flag=False,ss_spin_flag=False,
                        #lalinferenceMCMC options
                        li_flag=False,deltaLogP=None,fixedBurnins=None,nDownsample=None,oldMassConvention=False,
                        #followupMCMC options
                        fm_flag=False,
                        #spin frame for the injection
                        inj_spin_frame='OrbitalL',
                        # on ACF?
                        noacf=False,
                        #Turn on 2D kdes
                        twodkdeplots=False,
                        #Turn on R convergence tests
                        RconvergenceTests=False,
                        # Save PDF figures?
                        savepdfs=True,
                        #List of covariance matrix csv files used as analytic likelihood
                        covarianceMatrices=None,
                        #List of meanVector csv files used, one csv file for each covariance matrix
                        meanVectors=None,
                        #header file
                        header=None,
                        psd_files=None,
                        greedy=True ## If true will use greedy bin for 1-d credible regions. Otherwise use 2-steps KDE
                    ):
    """
    This is a demonstration script for using the functionality/data structures
    contained in pylal.bayespputils . It will produce a webpage from a file containing
    posterior samples generated by the parameter estimation codes with 1D/2D plots
    and stats from the marginal posteriors for each parameter/set of parameters.
    """
    votfile=None
    if eventnum is not None and injfile is None:
        print "You specified an event number but no injection file. Ignoring!"

    if trignum is not None and trigfile is None:
        print "You specified a trigger number but no trigger file. Ignoring!"

    if trignum is None and trigfile is not None:
        print "You specified a trigger file but no trigger number. Taking first entry (the case for GraceDB events)."
        trignum=0

    if data is None:
        raise RuntimeError('You must specify an input data file')
    #
    if outdir is None:
        raise RuntimeError("You must specify an output directory.")

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    #
    if fm_flag:
        peparser=bppu.PEOutputParser('fm')
        commonResultsObj=peparser.parse(data)

    elif ns_flag and not ss_flag:
        peparser=bppu.PEOutputParser('ns')
        commonResultsObj=peparser.parse(data,Nlive=ns_Nlive)

    elif ss_flag and not ns_flag:
        peparser=bppu.PEOutputParser('mcmc_burnin')
        commonResultsObj=peparser.parse(data,spin=ss_spin_flag,deltaLogP=deltaLogP)

    elif li_flag:
        peparser=bppu.PEOutputParser('inf_mcmc')
        commonResultsObj=peparser.parse(data,outdir=outdir,deltaLogP=deltaLogP,fixedBurnins=fixedBurnins,nDownsample=nDownsample,oldMassConvention=oldMassConvention)

    elif ss_flag and ns_flag:
        raise RuntimeError("Undefined input format. Choose only one of:")

    elif '.xml' in data[0]:
        peparser=bppu.PEOutputParser('xml')
        commonResultsObj=peparser.parse(data[0])
        thefile=open(data[0],'r')
        votfile=thefile.read()
    elif '.hdf' in data[0] or '.h5' in data[0]:
        if len(data) > 1:
            peparser = bppu.PEOutputParser('hdf5s')
            commonResultsObj=peparser.parse(data,deltaLogP=deltaLogP,fixedBurnins=fixedBurnins,nDownsample=nDownsample)
        else:
            fixedBurnins = fixedBurnins if fixedBurnins is not None else None
            peparser = bppu.PEOutputParser('hdf5')
            commonResultsObj=peparser.parse(data[0],deltaLogP=deltaLogP,fixedBurnin=fixedBurnins,nDownsample=nDownsample)
    else:
        peparser=bppu.PEOutputParser('common')
        commonResultsObj=peparser.parse(open(data[0],'r'),info=[header,None])
        # check if Nest (through nest2post) has produced an header file with git and CL info. If yes copy in outdir
        if os.path.isfile(data[0]+"_header.txt"):
          import shutil
          shutil.copy2(data[0]+"_header.txt", os.path.join(outdir,'nest_headers.txt'))

    #Extract f_ref from CRO if present.  This is needed to calculate orbital angular momentum
    #  when converting spin parameters.  Ideally this info will be provided in the
    #  SimInspiralTable in the near future.
    ps,samps = commonResultsObj
    try:
        f_refIdx = ps.index('f_ref')
        injFref = unique(samps[:,f_refIdx])
        if len(injFref) > 1:
            print "ERROR: Expected f_ref to be constant for all samples.  Can't tell which value was injected!"
            print injFref
            injFref = None
        else:
            injFref = injFref[0]
    except ValueError:
        injFref = None

    #Select injections using tc +/- 0.1s if it exists or eventnum from the injection file
    injection=None
    if injfile and eventnum is not None:
        print 'Looking for event %i in %s\n'%(eventnum,injfile)
        xmldoc = utils.load_filename(injfile,contenthandler=ExtractSimInspiralTableLIGOLWContentHandler)
        siminspiraltable=lsctables.table.get_table(xmldoc,lsctables.SimInspiralTable.tableName)
        injection=siminspiraltable[eventnum]
	#injections = SimInspiralUtils.ReadSimInspiralFromFiles([injfile])
	#if(len(injections)!=1): raise RuntimeError('Error: something unexpected happened while loading the injection file!\n')
        #injection=injections[0]

    #Get trigger
    triggers = None
    if trigfile is not None and trignum is not None:
        triggers = bppu.readCoincXML(trigfile, trignum)

    ## Load Bayes factors ##
    # Add Bayes factor information to summary file #
    if bayesfactornoise is not None:
        bfile=open(bayesfactornoise,'r')
        BSN=bfile.read()
        bfile.close()
        if(len(BSN.split())!=1):
          BSN=BSN.split()[0]
        print 'BSN: %s'%BSN
    if bayesfactorcoherent is not None:
        bfile=open(bayesfactorcoherent,'r')
        BCI=bfile.read()
        bfile.close()
        print 'BCI: %s'%BCI

    if snrfactor is not None:
        if not os.path.isfile(snrfactor):
            print "No snr file provided or wrong path to snr file\n"
            snrfactor=None
        else:
            snrstring=""
            snrfile=open(snrfactor,'r')
            snrs=snrfile.readlines()
            snrfile.close()
            for snr in snrs:
                if snr=="\n":
                    continue
                snrstring=snrstring +" "+str(snr[0:-1])+" ,"
            snrstring=snrstring[0:-1]

    #Create an instance of the posterior class using the posterior values loaded
    #from the file and any injection information (if given).
    pos = bppu.Posterior(commonResultsObj,SimInspiralTableEntry=injection,inj_spin_frame=inj_spin_frame, injFref=injFref,SnglInpiralList=triggers,votfile=votfile)

    #Create analytic likelihood functions if covariance matrices and mean vectors were given
    analyticLikelihood = None
    if covarianceMatrices and meanVectors:
        analyticLikelihood = bppu.AnalyticLikelihood(covarianceMatrices, meanVectors)

        #Plot only analytic parameters
        oneDMenu = analyticLikelihood.names
        if twoDGreedyMenu:
            twoDGreedyMenu = []
            for i in range(len(oneDMenu)):
                for j in range(i+1,len(oneDMenu)):
                    twoDGreedyMenu.append([oneDMenu[i],oneDMenu[j]])
        twoDplots = twoDGreedyMenu

    if eventnum is None and injfile is not None:
        import itertools
        injections = SimInspiralUtils.ReadSimInspiralFromFiles([injfile])

        if(len(injections)<1):
            try:
                print 'Warning: Cannot find injection with end time %f' %(pos['time'].mean)
            except KeyError:
                print "Warning: No 'time' column!"

        else:
            try:
                injection = itertools.ifilter(lambda a: abs(float(a.get_end()) - pos['time'].mean) < 0.1, injections).next()
                pos.set_injection(injection)
            except KeyError:
                print "Warning: No 'time' column!"

    pos.extend_posterior()
    # create list of names in oneDMenus dic
    oneDMenu=[]
    for i in oneDMenus.keys():
      oneDMenu+=oneDMenus[i]
    #Perform necessary mappings
    functions = {'cos':cos,'sin':sin,'exp':exp,'log':log}
    for pos_name in oneDMenu:
        if pos_name not in pos.names:
            for func in functions.keys():
                old_pos_name = pos_name.replace(func,'')
                if pos_name.find(func)==0 and old_pos_name in pos.names:
                    print "Taking %s of %s ..."% (func,old_pos_name)
                    pos.append_mapping(pos_name,functions[func],old_pos_name)

    #Remove samples with NaNs in requested params
    requested_params = set(pos.names).intersection(set(oneDMenu))
    pos.delete_NaN_entries(requested_params)

    #Remove non-analytic parameters if analytic likelihood is given:
    if analyticLikelihood:
        dievidence_names = ['post','posterior','logl','prior','likelihood','cycle','chain']
        [pos.pop(param) for param in pos.names if param not in analyticLikelihood.names and param not in dievidence_names]

    ##Print some summary stats for the user...##
    #Number of samples
    print "Number of posterior samples: %i"%len(pos)
    # Means
    print 'Means:'
    print str(pos.means)
    #Median
    print 'Median:'
    print str(pos.medians)
    #maxL
    print 'maxL:'
    max_pos,max_pos_co=pos.maxL
    print max_pos_co

    # Save posterior samples
    posfilename=os.path.join(outdir,'posterior_samples.dat')
    pos.write_to_file(posfilename)

    #==================================================================#
    #Create web page
    #==================================================================#

    html=bppu.htmlPage('Posterior PDFs',css=bppu.__default_css_string,javascript=bppu.__default_javascript_string)

    #Create a section for meta-data/run information
    html_meta=html.add_section('Summary')
    table=html_meta.tab()
    row=html_meta.insert_row(table,label='thisrow')
    td=html_meta.insert_td(row,'',label='Samples')
    SampsStats=html.add_section_to_element('Samples',td)
    SampsStats.p('Produced from '+str(len(pos))+' posterior samples.')
    if 'chain' in pos.names:
        acceptedChains = unique(pos['chain'].samples)
        acceptedChainText = '%i of %i chains accepted: %i'%(len(acceptedChains),len(data),acceptedChains[0])
        if len(acceptedChains) > 1:
            for chain in acceptedChains[1:]:
                acceptedChainText += ', %i'%(chain)
        SampsStats.p(acceptedChainText)
    if 'cycle' in pos.names:
        SampsStats.p('Longest chain has '+str(pos.longest_chain_cycles())+' cycles.')
    filenames='Samples read from %s'%(data[0])
    if len(data) > 1:
        for fname in data[1:]:
            filenames+=', '+str(fname)
    SampsStats.p(filenames)
    td=html_meta.insert_td(row,'',label='SummaryLinks')
    legend=html.add_section_to_element('Sections',td)

    # Create a section for HDF5 metadata if available
    if '.h5' in data[0] or '.hdf' in data[0]:
        html_hdf=html.add_section('Metadata',legend=legend)
        import h5py
        with h5py.File(data[0],'r') as h5grp:
            extract_hdf5_metadata(h5grp,parent=html_hdf)

    #Create a section for model selection results (if they exist)
    if bayesfactorcoherent is not None or bayesfactornoise is not None:
        html_model=html.add_section('Model selection',legend=legend)
        if bayesfactornoise is not None:
            html_model.p('log Bayes factor ( coherent vs gaussian noise) = %s, Bayes factor=%f'%(BSN,exp(float(BSN))))
        if bayesfactorcoherent is not None:
            html_model.p('log Bayes factor ( coherent vs incoherent OR noise ) = %s, Bayes factor=%f'%(BCI,exp(float(BCI))))

    if dievidence:
        html_model=html.add_section('Direct Integration Evidence',legend=legend)
        log_ev = log(difactor) + pos.di_evidence(boxing=boxing)
        ev=exp(log_ev)
        evfilename=os.path.join(outdir,"evidence.dat")
        evout=open(evfilename,"w")
        evout.write(str(ev))
        evout.write(" ")
        evout.write(str(log_ev))
        evout.close()
        print "Computing direct integration evidence = %g (log(Evidence) = %g)"%(ev, log_ev)
        html_model.p('Direct integration evidence is %g, or log(Evidence) = %g.  (Boxing parameter = %d.)'%(ev,log_ev,boxing))
        if 'logl' in pos.names:
            log_ev=pos.harmonic_mean_evidence()
            html_model.p('Compare to harmonic mean evidence of %g (log(Evidence) = %g).'%(exp(log_ev),log_ev))

    if ellevidence:
        try:
            html_model=html.add_section('Elliptical Evidence',legend=legend)
            log_ev = pos.elliptical_subregion_evidence()
            ev = exp(log_ev)
            evfilename=os.path.join(outdir, 'ellevidence.dat')
            evout=open(evfilename, 'w')
            evout.write(str(ev) + ' ' + str(log_ev))
            evout.close()
            print 'Computing elliptical region evidence = %g (log(ev) = %g)'%(ev, log_ev)
            html_model.p('Elliptical region evidence is %g, or log(Evidence) = %g.'%(ev, log_ev))

            if 'logl' in pos.names:
                log_ev=pos.harmonic_mean_evidence()
                html_model.p('Compare to harmonic mean evidence of %g (log(Evidence = %g))'%(exp(log_ev), log_ev))
        except IndexError:
            print 'Warning: Sample size too small to compute elliptical evidence!'

    #Create a section for SNR, if a file is provided
    if snrfactor is not None:
        html_snr=html.add_section('Signal to noise ratio(s)',legend=legend)
        html_snr.p('%s'%snrstring)

    # Create a section for the DIC
    html_dic = html.add_section('Deviance Information Criterion', legend=legend)
    html_dic.p('DIC = %.1f'%pos.DIC)
    dicout = open(os.path.join(outdir, 'dic.dat'), 'w')
    try:
        dicout.write('%.1f\n'%pos.DIC)
    finally:
        dicout.close()

    #Create a section for summary statistics
    # By default collapse section are collapsed when the page is opened or reloaded. Use start_closed=False option as here below to change this
    tabid='statstable'
    html_stats=html.add_collapse_section('Summary statistics',legend=legend,innertable_id=tabid,start_closed=False)
    html_stats.write(str(pos))
    statfilename=os.path.join(outdir,"summary_statistics.dat")
    statout=open(statfilename,"w")
    statout.write("\tmaP\tmaxL\tKL\tstdev\tmean\tmedian\tstacc\tinjection\tvalue\n")

    warned_about_kl = False
    for statname,statoned_pos in pos:

      statmax_pos,max_i=pos._posMaxL()
      statmaxL=statoned_pos.samples[max_i][0]
      try:
          statKL = statoned_pos.KL()
      except ValueError:
          if not warned_about_kl:
              print("Unable to compute KL divergence")
              warned_about_kl = True
          statKL = None

      statmax_pos,max_j=pos._posMap()
      statmaxP=statoned_pos.samples[max_j][0]
      statmean=str(statoned_pos.mean)
      statstdev=str(statoned_pos.stdev)
      statmedian=str(squeeze(statoned_pos.median))
      statstacc=str(statoned_pos.stacc)
      statinjval=str(statoned_pos.injval)

      statarray=[str(i) for i in [statname,statmaxP,statmaxL,statKL,statstdev,statmean,statmedian,statstacc,statinjval]]
      statout.write("\t".join(statarray))
      statout.write("\n")

    statout.close()

    #==================================================================#
    #Generate sky map, WF, and PSDs
    #==================================================================#

    skyreses=None
    sky_injection_cl=None
    inj_position=None
    tabid='skywftable'
    html_wf=html.add_collapse_section('Sky Localization and Waveform',innertable_id=tabid)

    table=html_wf.tab(idtable=tabid)
    row=html_wf.insert_row(table,label='SkyandWF')
    skytd=html_wf.insert_td(row,'',label='SkyMap',legend=legend)
    html_sky=html.add_section_to_element('SkyMap',skytd)
    #If sky resolution parameter has been specified try and create sky map...
    if skyres is not None and \
       ('ra' in pos.names and 'dec' in pos.names):

        if pos['dec'].injval is not None and pos['ra'].injval is not None:
            inj_position=[pos['ra'].injval,pos['dec'].injval]
        else:
            inj_position=None

        hpmap = pos.healpix_map(float(skyres), nest=True)
        bppu.plot_sky_map(hpmap, outdir, inj=inj_position, nest=True)

        if inj_position is not None:
            html_sky.p('Injection found at p = %g'%bppu.skymap_inj_pvalue(hpmap, inj_position, nest=True))

        html_sky.write('<a href="skymap.png" target="_blank"><img src="skymap.png"/></a>')

        html_sky_write='<table border="1" id="statstable"><tr><th>Confidence region</th><th>size (sq. deg)</th></tr>'

        areas = bppu.skymap_confidence_areas(hpmap, confidence_levels)
        for cl, area in zip(confidence_levels, areas):
            html_sky_write+='<tr><td>%g</td><td>%g</td></tr>'%(cl, area)
        html_sky_write+=('</table>')

        html_sky.write(html_sky_write)
    else:
        html_sky.write('<b>No skymap generated!</b>')
    wfdir=os.path.join(outdir,'Waveform')
    if not os.path.isdir(wfdir):
        os.makedirs(wfdir)
    try:
        wfpointer= bppu.plot_waveform(pos=pos,siminspiral=injfile,event=eventnum,path=wfdir)
    except  Exception,e:
        wfpointer = None
        print "Could not create WF plot. The error was: %s\n"%str(e)
    wftd=html_wf.insert_td(row,'',label='Waveform',legend=legend)
    wfsection=html.add_section_to_element('Waveforms',wftd)
    if wfpointer is not None:
      wfsection.write('<a href="Waveform/WF_DetFrame.png" target="_blank"><img src="Waveform/WF_DetFrame.png"/></a>')
    else:
      print "Could not create WF plot.\n"
      wfsection.write("<b>No Waveform generated!</b>")

    wftd=html_wf.insert_td(row,'',label='PSDs',legend=legend)
    wfsection=html.add_section_to_element('PSDs',wftd)
    if psd_files is not None:
      psd_files=list(psd_files.split(','))
      psddir=os.path.join(outdir,'PSDs')
      if not os.path.isdir(psddir):
        os.makedirs(psddir)
      try:
        if 'flow' in pos.names:
          f_low = pos['flow'].samples.min()
        else:
          f_low = 30.
        bppu.plot_psd(psd_files,outpath=psddir, f_min=f_low)
        wfsection.write('<a href="PSDs/PSD.png" target="_blank"><img src="PSDs/PSD.png"/></a>')
      except  Exception,e:
        print "Could not create PSD plot. The error was: %s\n"%str(e)
        wfsection.write("<b>PSD plotting failed</b>")
    else:
        wfsection.write("<b>No PSD files provided</b>")

    # Add plots for calibration estimates
    if np.any(['spcal_amp' in param for param in pos.names]) or np.any(['spcal_phase' in param for param in pos.names]):
      wftd=html_wf.insert_td(row,'',label='Calibration',legend=legend)
      wfsection=html.add_section_to_element('Calibration',wftd)
      bppu.plot_calibration_pos(pos, outpath=outdir)
      wfsection.write('<a href="calibration.png" target="_blank"><img src="calibration.png"/></a>')
     # if precessing spins do spin disk
    allin=1.0
    for i in ['a1','tilt_spin1','a2','tilt_spin2']:
      if not i in pos.names:
        allin*=0.0
    if allin ==0.0:
      pass
    else:
      wftd=html_wf.insert_td(row,'',label='DiskPlot',legend=legend)
      wfsection=html.add_section_to_element('DiskPlot',wftd)
      cbcdiskspin.make_disk_plot(pos, outpath=outdir)
      wfsection.write('<a href="comp_spin_pos.png" target="_blank"><img src="comp_spin_pos.png"/></a>')

    #==================================================================#
    #1D posteriors
    #==================================================================#
    onepdfdir=os.path.join(outdir,'1Dpdf')
    if not os.path.isdir(onepdfdir):
        os.makedirs(onepdfdir)

    sampsdir=os.path.join(outdir,'1Dsamps')
    if not os.path.isdir(sampsdir):
        os.makedirs(sampsdir)
    reses={}

    for i in oneDMenus.keys():
      rss=bppu.make_1d_table(html,legend,i,pos,oneDMenus[i],noacf,GreedyRes,onepdfdir,sampsdir,savepdfs,greedy,analyticLikelihood,nDownsample)
      reses.update(rss)




    tabid='onedconftable'
    html_ogci=html.add_collapse_section('1D confidence intervals (greedy binning)',legend=legend,innertable_id=tabid)
    html_ogci_write='<table id="%s" border="1"><tr><th/>'%tabid
    clasciiout="#parameter \t"
    confidence_levels.sort()
    for cl in confidence_levels:
        html_ogci_write+='<th>%f</th>'%cl
        clasciiout+="%s\t"%('%.02f'%cl)
    if injection:
        html_ogci_write+='<th>Injection Confidence Level</th>'
        html_ogci_write+='<th>Injection Confidence Interval</th>'
        clasciiout+="Injection_Confidence_Level\t"
        clasciiout+="Injection_Confidence_Interval"
    clasciiout+='\n'
    html_ogci_write+='</tr>'
    #Generate new BCI html table row
    printed=0
    for par_name in oneDMenu:
      par_name=par_name.lower()
      try:
              pos[par_name.lower()]
      except KeyError:
      #print "No input chain for %s, skipping binning."%par_name
             continue
      try:
              par_bin=GreedyRes[par_name]
      except KeyError:
              print "Bin size is not set for %s, skipping binning."%par_name
              continue
      binParams={par_name:par_bin}
      injection_area=None
      if greedy:
                    if printed==0:
                        print "Using greedy 1-d binning credible regions\n"
                        printed=1
                    toppoints,injectionconfidence,reses,injection_area,cl_intervals=bppu.greedy_bin_one_param(pos,binParams,confidence_levels)
      else:
                    if printed==0:
                        print "Using 2-step KDE 1-d credible regions\n"
                        printed=1
                    if pos[par_name].injval is None:
                        injCoords=None
                    else:
                        injCoords=[pos[par_name].injval]
                    _,reses,injstats=bppu.kdtree_bin2Step(pos,[par_name],confidence_levels,injCoords=injCoords)
                    if injstats is not None:
                        injectionconfidence=injstats[3]
                        injection_area=injstats[4]

      BCItableline='<tr><td>%s</td>'%(par_name)
      clasciiout+="%s\t"%par_name
      cls=reses.keys()
      cls.sort()

      for cl in cls:
          BCItableline+='<td>%f</td>'%reses[cl]
          clasciiout+="%f\t"%reses[cl]
      if injection is not None:
          if injectionconfidence is not None and injection_area is not None:

              BCItableline+='<td>%f</td>'%injectionconfidence
              BCItableline+='<td>%f</td>'%injection_area
              clasciiout+="%f\t"%injectionconfidence
              clasciiout+="%f"%injection_area

          else:
              BCItableline+='<td/>'
              BCItableline+='<td/>'
              clasciiout+="nan\t"
              clasciiout+="nan"
      BCItableline+='</tr>'
      clasciiout+="\n"
      #Append new table line to section html
      html_ogci_write+=BCItableline

    html_ogci_write+='</table>'
    html_ogci.write(html_ogci_write)

    #===============================#
    # Corner plots
    #===============================#
    cornerdir=os.path.join(outdir,'corner')
    if not os.path.isdir(cornerdir):
      os.makedirs(cornerdir)
    massParams=['mtotal','m1','m2','mc']
    distParams=['distance','distMPC','dist']
    incParams=['iota','inclination','theta_jn']
    polParams=['psi','polarisation','polarization']
    skyParams=['ra','rightascension','declination','dec']
    timeParams=['time']
    spinParams=['spin1','spin2','a1','a2','a1z','a2z','phi1','theta1','phi2','theta2','chi','effectivespin','chi_eff','chi_tot','chi_p','beta','tilt1','tilt2','phi_jl','theta_jn','phi12']
    sourceParams=['m1_source','m2_source','mtotal_source','mc_source','redshift']
    intrinsicParams=massParams+spinParams
    extrinsicParams=incParams+distParams+polParams+skyParams
    sourceFrameParams=sourceParams+distParams
    try:
      myfig=bppu.plot_corner(pos,[0.05,0.5,0.95],parnames=intrinsicParams)
    except:
      myfig=None
    tabid='CornerTable'
    html_corner=''
    got_any=0
    if myfig:
      html_corner+='<tr><td width="100%"><a href="corner/intrinsic.png" target="_blank"><img width="70%" src="corner/intrinsic.png"/></a></td></tr>'
      myfig.savefig(os.path.join(cornerdir,'intrinsic.png'))
      myfig.savefig(os.path.join(cornerdir,'intrinsic.pdf'))
      got_any+=1
    try:
      myfig=bppu.plot_corner(pos,[0.05,0.5,0.95],parnames=extrinsicParams)
    except:
      myfig=None

    if myfig:
      myfig.savefig(os.path.join(cornerdir,'extrinsic.png'))
      myfig.savefig(os.path.join(cornerdir,'extrinsic.pdf'))
      html_corner+='<tr><td width="100%"><a href="corner/extrinsic.png" target="_blank"><img width="70%" src="corner/extrinsic.png"/></a></td></tr>'
      got_any+=1
    try:
      myfig=bppu.plot_corner(pos,[0.05,0.5,0.95],parnames=sourceFrameParams)
    except:
      myfig=None

    if myfig:
      myfig.savefig(os.path.join(cornerdir,'sourceFrame.png'))
      myfig.savefig(os.path.join(cornerdir,'sourceFrame.pdf'))
      html_corner+='<tr><td width="100%"><a href="corner/sourceFrame.png" target="_blank"><img width="70%" src="corner/sourceFrame.png"/></a></td></tr>'
      got_any+=1

    if got_any>0:
      html_corner='<table id="%s" border="1">'%tabid+html_corner
      html_corner+='</table>'
    if html_corner!='':
      html_co=html.add_collapse_section('Corner plots',legend=legend,innertable_id=tabid)
      html_co.write(html_corner)
    if clasciiout:
      fout=open(os.path.join(outdir,'confidence_levels.txt'),'w')
      fout.write(clasciiout)
      fout.close()
    #==================================================================#
    #2D posteriors
    #==================================================================#

    #Loop over parameter pairs in twoDGreedyMenu and bin the sample pairs
    #using a greedy algorithm . The ranked pixels (toppoints) are used
    #to plot 2D histograms and evaluate Bayesian confidence intervals.

    #Make a folder for the 2D kde plots
    margdir=os.path.join(outdir,'2Dkde')
    if not os.path.isdir(margdir):
        os.makedirs(margdir)

    twobinsdir=os.path.join(outdir,'2Dbins')
    if not os.path.isdir(twobinsdir):
        os.makedirs(twobinsdir)

    greedytwobinsdir=os.path.join(outdir,'greedy2Dbins')
    if not os.path.isdir(greedytwobinsdir):
        os.makedirs(greedytwobinsdir)

    #Add a section to the webpage for a table of the confidence interval
    #results.
    tabid='2dconftable'
    html_tcig=html.add_collapse_section('2D confidence intervals (greedy binning)',legend=legend,innertable_id=tabid)
    #Generate the top part of the table
    html_tcig_write='<table id="%s" border="1"><tr><th/>'%tabid
    confidence_levels.sort()
    for cl in confidence_levels:
        html_tcig_write+='<th>%f</th>'%cl
    if injection:
        html_tcig_write+='<th>Injection Confidence Level</th>'
        html_tcig_write+='<th>Injection Confidence Interval</th>'
    html_tcig_write+='</tr>'


    #=  Add a section for a table of 2D marginal PDFs (kde)
    twodkdeplots_flag=twodkdeplots
    if twodkdeplots_flag:
        tabid='2dmargtable'
        html_tcmp=html.add_collapse_section('2D Marginal PDFs',legend=legend,innertable_id=tabid)
        #Table matter
        html_tcmp_write='<table border="1" id="%s">'%tabid

    tabid='2dgreedytable'
    html_tgbh=html.add_collapse_section('2D Greedy Bin Histograms',legend=legend,innertable_id=tabid)
    html_tgbh_write='<table border="1" id="%s">'%tabid

    row_count=0
    row_count_gb=0

    for par1_name,par2_name in twoDGreedyMenu:
        par1_name=par1_name.lower()
        par2_name=par2_name.lower()
        # Don't plot a parameter against itself!
        if par1_name == par2_name: continue
        try:
            pos[par1_name.lower()]
        except KeyError:
            #print "No input chain for %s, skipping binning."%par1_name
            continue
        try:
            pos[par2_name.lower()]
        except KeyError:
            #print "No input chain for %s, skipping binning."%par2_name
            continue
        #Bin sizes
        try:
            par1_bin=GreedyRes[par1_name]
        except KeyError:
            print "Bin size is not set for %s, skipping %s/%s binning."%(par1_name,par1_name,par2_name)
            continue
        try:
            par2_bin=GreedyRes[par2_name]
        except KeyError:
            print "Bin size is not set for %s, skipping %s/%s binning."%(par2_name,par1_name,par2_name)
            continue

        # Skip any fixed parameters to avoid matrix inversion problems
        par1_pos=pos[par1_name].samples
        par2_pos=pos[par2_name].samples
        if (size(unique(par1_pos))<2 or size(unique(par2_pos))<2):
            continue

        #print "Binning %s-%s to determine confidence levels ..."%(par1_name,par2_name)
        #Form greedy binning input structure
        greedy2Params={par1_name:par1_bin,par2_name:par2_bin}

        #Greedy bin the posterior samples
        toppoints,injection_cl,reses,injection_area=\
        bppu.greedy_bin_two_param(pos,greedy2Params,confidence_levels)

        print "BCI %s-%s:"%(par1_name,par2_name)
        print reses

        #Generate new BCI html table row
        BCItableline='<tr><td>%s-%s</td>'%(par1_name,par2_name)
        cls=reses.keys()
        cls.sort()

        for cl in cls:
            BCItableline+='<td>%f</td>'%reses[cl]

        if injection is not None:
            if injection_cl is not None:
                BCItableline+='<td>%f</td>'%injection_cl
                BCItableline+='<td>'+str(injection_area)+'</td>'

            else:
                BCItableline+='<td/>'
                BCItableline+='<td/>'

        BCItableline+='</tr>'

        #Append new table line to section html
        html_tcig_write+=BCItableline


        #= Plot 2D histograms of greedily binned points =#

        #greedy2ContourPlot=bppu.plot_two_param_greedy_bins_contour({'Result':pos},greedy2Params,[0.67,0.9,0.95],{'Result':'k'})
        greedy2ContourPlot=bppu.plot_two_param_kde_greedy_levels({'Result':pos},greedy2Params,[0.67,0.9,0.95],{'Result':'k'})
        greedy2contourpath=os.path.join(greedytwobinsdir,'%s-%s_greedy2contour.png'%(par1_name,par2_name))
        if greedy2ContourPlot is not None:
          greedy2ContourPlot.savefig(greedy2contourpath)
          if(savepdfs): greedy2ContourPlot.savefig(greedy2contourpath.replace('.png','.pdf'))
          plt.close(greedy2ContourPlot)

        greedy2HistFig=bppu.plot_two_param_greedy_bins_hist(pos,greedy2Params,confidence_levels)
        greedy2histpath=os.path.join(greedytwobinsdir,'%s-%s_greedy2.png'%(par1_name,par2_name))
        greedy2HistFig.savefig(greedy2histpath)
        if(savepdfs): greedy2HistFig.savefig(greedy2histpath.replace('.png','.pdf'))
        plt.close(greedy2HistFig)

        greedyFile = open(os.path.join(twobinsdir,'%s_%s_greedy_stats.txt'%(par1_name,par2_name)),'w')

        #= Write out statistics for greedy bins
        for cl in cls:
            greedyFile.write("%lf %lf\n"%(cl,reses[cl]))
        greedyFile.close()

        if [par1_name,par2_name] in twoDplots or [par2_name,par1_name] in twoDplots :
            print 'Generating %s-%s greedy hist plot'%(par1_name,par2_name)

            par1_pos=pos[par1_name].samples
            par2_pos=pos[par2_name].samples

            if (size(unique(par1_pos))<2 or size(unique(par2_pos))<2):
                continue
            head,figname=os.path.split(greedy2histpath)
            head,figname_c=os.path.split(greedy2contourpath)
            if row_count_gb==0:
                html_tgbh_write+='<tr>'
            html_tgbh_write+='<td width="30%"><img width="100%" src="greedy2Dbins/'+figname+'"/>[<a href="greedy2Dbins/'+figname_c+'">contour</a>]</td>'
            row_count_gb+=1
            if row_count_gb==3:
                html_tgbh_write+='</tr>'
                row_count_gb=0

        #= Generate 2D kde plots =#

        if twodkdeplots_flag is True:
            if [par1_name,par2_name] in twoDplots or [par2_name,par1_name] in twoDplots :
                print 'Generating %s-%s plot'%(par1_name,par2_name)

                par1_pos=pos[par1_name].samples
                par2_pos=pos[par2_name].samples

                if (size(unique(par1_pos))<2 or size(unique(par2_pos))<2):
                    continue

                plot2DkdeParams={par1_name:50,par2_name:50}
                myfig=bppu.plot_two_param_kde(pos,plot2DkdeParams)

                figname=par1_name+'-'+par2_name+'_2Dkernel.png'
                twoDKdePath=os.path.join(margdir,figname)

                if row_count==0:
                    html_tcmp_write+='<tr>'
                html_tcmp_write+='<td width="30%"><img width="100%" src="2Dkde/'+figname+'"/></td>'
                row_count+=1
                if row_count==3:
                    html_tcmp_write+='</tr>'
                    row_count=0

                if myfig:
                  myfig.savefig(twoDKdePath)
                  if(savepdfs): myfig.savefig(twoDKdePath.replace('.png','.pdf'))
                  plt.close(myfig)
                else:
                  print 'Unable to generate 2D kde levels for %s-%s'%(par1_name,par2_name)


    #Finish off the BCI table and write it into the etree
    html_tcig_write+='</table>'
    html_tcig.write(html_tcig_write)

    if twodkdeplots_flag is True:
    #Finish off the 2D kde plot table
        while row_count!=0:
            html_tcmp_write+='<td/>'
            row_count+=1
            if row_count==3:
                row_count=0
                html_tcmp_write+='</tr>'
        html_tcmp_write+='</table>'
        html_tcmp.write(html_tcmp_write)
        #Add a link to all plots
        html_tcmp.a("2Dkde/",'All 2D marginal PDFs (kde)')

    #Finish off the 2D greedy histogram plot table
    while row_count_gb!=0:
        html_tgbh_write+='<td/>'
        row_count_gb+=1
        if row_count_gb==3:
            row_count_gb=0
            html_tgbh_write+='</tr>'
    html_tgbh_write+='</table>'
    html_tgbh.write(html_tgbh_write)
    #Add a link to all plots
    html_tgbh.a("greedy2Dbins/",'All 2D Greedy Bin Histograms')

    if RconvergenceTests is True:
        convergenceResults=bppu.convergenceTests(pos,gelman=False)

        if convergenceResults is not None:
            tabid='convtable'
            html_conv_test=html.add_collapse_section('Convergence tests',legend=legend,innertable_id=tabid)
            data_found=False
            for test,test_data in convergenceResults.items():

                if test_data:
                    data_found=True
                    html_conv_test.h3(test)

                    html_conv_table_rows={}
                    html_conv_table_header=''
                    for chain,chain_data in test_data.items():
                        html_conv_table_header+='<th>%s</th>'%chain


                        for data in chain_data:
                            if len(data)==2:
                                try:
                                    html_conv_table_rows[data[0]]+='<td>'+data[1]+'</td>'
                                except KeyError:
                                    html_conv_table_rows[data[0]]='<td>'+data[1]+'</td>'

                    html_conv_table='<table id="%s"><tr><th>Chain</th>'%tabid+html_conv_table_header+'</tr>'
                    for row_name,row in html_conv_table_rows.items():
                        html_conv_table+='<tr><td>%s</td>%s</tr>'%(row_name,row)
                    html_conv_table+='</table>'
                    html_conv_test.write(html_conv_table)
            if data_found is False:
                html_conv_test.p('No convergence diagnostics generated!')

    #Create a section for the covariance matrix
    tabid='covtable'
    html_stats_cov=html.add_collapse_section('Covariance matrix',legend=legend,innertable_id=tabid)
    pos_samples,table_header_string=pos.samples()

    #calculate cov matrix
    try:
        cov_matrix=cov(pos_samples,rowvar=0,bias=1)

    #Create html table
        table_header_list=table_header_string.split()
        cov_table_string='<table border="1" id="%s"><tr><th/>'%tabid
        for header in table_header_list:
            cov_table_string+='<th>%s</th>'%header
        cov_table_string+='</tr>'

        cov_column_list=hsplit(cov_matrix,cov_matrix.shape[1])

        for cov_column,cov_column_name in zip(cov_column_list,table_header_list):
            cov_table_string+='<tr><th>%s</th>'%cov_column_name
            for cov_column_element in cov_column:
                cov_table_string+='<td>%.3e</td>'%(cov_column_element[0])
            cov_table_string+='</tr>'
        cov_table_string+='</table>'
        html_stats_cov.write(cov_table_string)

    except:
        print 'Unable to compute the covariance matrix.'

    #Create a section for run configuration information if it exists
    if pos._votfile is not None:
        html_vot=html.add_section('Run information',legend=legend)
        html_vot.write(pos.write_vot_info())

    html_footer=html.add_section('')
    html_footer.p('Produced using cbcBayesPostProc.py at '+strftime("%Y-%m-%d %H:%M:%S")+' .')

    cc_args=''
    for arg in sys.argv:
        cc_args+=arg+' '

    html_footer.p('Command line: %s'%cc_args)
    html_footer.p(git_version.verbose_msg)

    #Save results page
    resultspage=open(os.path.join(outdir,'posplots.html'),'w')
    resultspage.write(str(html))


    #Close files
    resultspage.close()

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
    #Optional (all)
    parser.add_option("-i","--inj",dest="injfile",help="SimInsipral injection file",metavar="INJ.XML",default=None)
    parser.add_option("-t","--trig",dest="trigfile",help="Coinc XML file",metavar="COINC.XML",default=None)
    parser.add_option("--skyres",dest="skyres",help="Sky resolution to use to calculate sky box size",default=None)
    parser.add_option("--eventnum",dest="eventnum",action="store",default=None,help="event number in SimInspiral file of this signal",type="int",metavar="NUM")
    parser.add_option("--trignum",dest="trignum",action="store",default=None,help="trigger number in CoincTable",type="int",metavar="NUM")
    parser.add_option("--bsn",action="store",default=None,help="Optional file containing the bayes factor signal against noise",type="string")
    parser.add_option("--bci",action="store",default=None,help="Optional file containing the bayes factor coherent signal model against incoherent signal model.",type="string")
    parser.add_option("--snr",action="store",default=None,help="Optional file containing the SNRs of the signal in each IFO",type="string")
    parser.add_option("--dievidence",action="store_true",default=False,help="Calculate the direct integration evidence for the posterior samples")
    parser.add_option("--boxing",action="store",default=64,help="Boxing parameter for the direct integration evidence calculation",type="int",dest="boxing")
    parser.add_option("--evidenceFactor",action="store",default=1.0,help="Overall factor (normalization) to apply to evidence",type="float",dest="difactor",metavar="FACTOR")

    parser.add_option('--ellipticEvidence', action='store_true', default=False,help='Estimate the evidence by fitting ellipse to highest-posterior points.', dest='ellevidence')

    parser.add_option("--plot-2d", action="store_true", default=False,help="Make individual 2-D plots.")
    parser.add_option("--header",action="store",default=None,help="Optional file containing the header line for posterior samples",type="string")
    #NS
    parser.add_option("--ns",action="store_true",default=False,help="(inspnest) Parse input as if it was output from parallel nested sampling runs.")
    parser.add_option("--Nlive",action="store",default=None,help="(inspnest) Number of live points used in each parallel nested sampling run.",type="int")
    parser.add_option("--xflag",action="store_true",default=False,help="(inspnest) Convert x to iota.")
    #SS
    parser.add_option("--ss",action="store_true",default=False,help="(SPINspiral) Parse input as if it was output from SPINspiral.")
    parser.add_option("--spin",action="store_true",default=False,help="(SPINspiral) Specify spin run (15 parameters). ")
    #LALInf
    parser.add_option("--lalinfmcmc",action="store_true",default=False,help="(LALInferenceMCMC) Parse input from LALInferenceMCMC.")
    parser.add_option("--inj-spin-frame",default='OrbitalL', help="The reference frame used for the injection (default: OrbitalL)")
    parser.add_option("--downsample",action="store",default=None,help="(LALInferenceMCMC) approximate number of samples to record in the posterior",type="int")
    parser.add_option("--deltaLogL",action="store",default=None,help="(LALInferenceMCMC) Difference in logL to use for convergence test. (DEPRECATED)",type="float")
    parser.add_option("--deltaLogP",action="store",default=None,help="(LALInferenceMCMC) Difference in logpost to use for burnin criteria.",type="float")
    parser.add_option("--fixedBurnin",dest="fixedBurnin",action="callback",callback=multipleFileCB,help="(LALInferenceMCMC) Fixed number of iteration for burnin.")
    parser.add_option("--oldMassConvention",action="store_true",default=False,help="(LALInferenceMCMC) if activated, m2 > m1; otherwise m1 > m2 in PTMCMC.output.*.00")
    #FM
    parser.add_option("--fm",action="store_true",default=False,help="(followupMCMC) Parse input as if it was output from followupMCMC.")
    # ACF plots off?
    parser.add_option("--no-acf", action="store_true", default=False, dest="noacf")
    # Turn on 2D kdes
    parser.add_option("--twodkdeplots", action="store_true", default=False, dest="twodkdeplots")
    # Turn on R convergence tests
    parser.add_option("--RconvergenceTests", action="store_true", default=False, dest="RconvergenceTests")
    parser.add_option("--nopdfs",action="store_false",default=True,dest="nopdfs")
    parser.add_option("-c","--covarianceMatrix",dest="covarianceMatrices",action="append",default=None,help="CSV file containing covariance (must give accompanying mean vector CSV. Can add more than one matrix.")
    parser.add_option("-m","--meanVectors",dest="meanVectors",action="append",default=None,help="Comma separated list of locations of the multivariate gaussian described by the correlation matrix.  First line must be list of params in the order used for the covariance matrix.  Provide one list per covariance matrix.")
    parser.add_option("--email",action="store",default=None,type="string",metavar="user@ligo.org",help="Send an e-mail to the given address with a link to the finished page.")
    parser.add_option("--archive",action="store",default=None,type="string",metavar="results.tar.gz",help="Create the given tarball with all results")
    parser.add_option("--psdfiles",action="store",default=None,type="string",metavar="H1,L1,V1",help="comma separater list of ASCII files with PSDs, one per IFO")
    parser.add_option("--kdecredibleregions",action="store_true",default=False,help="If given, will use 2-step KDE trees to estimate 1-d credible regions [default false: use greedy binning]")
    parser.add_option("--noplot-source-frame", action="store_true", default=False,help="Don't make 1D plots of source-frame masses")
    (opts,args)=parser.parse_args()

    datafiles=[]
    if args:
      datafiles=datafiles+args
    if opts.data:
      datafiles=datafiles + opts.data

    if opts.fixedBurnin:
      # If only one value for multiple chains, assume it's to be applied to all chains
      if len(opts.fixedBurnin) == 1:
        fixedBurnins = [int(opts.fixedBurnin[0]) for df in datafiles]
      else:
        fixedBurnins = [int(fixedBurnin) for fixedBurnin in opts.fixedBurnin]
    else:
      fixedBurnins = None

    import pylal
    from pylal.bayespputils import massParams,spinParams,cosmoParam,strongFieldParams,distParams,incParams,polParams,skyParams,phaseParams,timeParams,endTimeParams,statsParams,calibParams,snrParams,tidalParams

    oneDMenus={'Masses':None,'SourceFrame':None,'Timing':None,'Extrinsic':None,'Spins':None,'StrongField':None,'Others':None}

    oneDMenus['Masses']= massParams
    oneDMenus['Extrinsic']=incParams+distParams+polParams+skyParams+phaseParams
    oneDMenus['Spins']= spinParams
    oneDMenus['Timing']=timeParams+endTimeParams
    oneDMenus['StrongField']= strongFieldParams
    oneDMenus['Others']=snrParams+statsParams+calibParams
    oneDMenus['SourceFrame']= cosmoParam

    if opts.noplot_source_frame:
      oneDMenus['SourceFrame']= []

    ifos_menu=['h1','l1','v1']
    from itertools import combinations
    for ifo1,ifo2 in combinations(ifos_menu,2):
      oneDMenus['Timing'].append(ifo1+ifo2+'_delay')
    #oneDMenu=[]
    twoDGreedyMenu=[]
    if opts.plot_2d:
        for mp1,mp2 in combinations(massParams,2):
          twoDGreedyMenu.append([mp1, mp2])
        for mp in massParams:
            for d in distParams:
                twoDGreedyMenu.append([mp,d])
        for mp in massParams:
            for sp in spinParams:
                twoDGreedyMenu.append([mp,sp])
        for mp in massParams:
            for dchi in bppu.tigerParams:
                twoDGreedyMenu.append([mp,dchi])
        for dp in distParams:
            for sp in snrParams:
                twoDGreedyMenu.append([dp,sp])
        for dp in distParams:
            for ip in incParams:
                twoDGreedyMenu.append([dp,ip])
        for dp in distParams:
            for sp in skyParams:
                twoDGreedyMenu.append([dp,sp])
        for dp in distParams:
            for sp in spinParams:
                twoDGreedyMenu.append([dp,sp])
        for ip in incParams:
            for sp in skyParams:
                twoDGreedyMenu.append([ip,sp])
        for ip in incParams:
            for sp in spinParams:
                twoDGreedyMenu.append([ip,sp])
            for phip in phaseParams:
                twoDGreedyMenu.append([ip,phip])
            for psip in polParams:
                twoDGreedyMenu.append([ip,psip])
        for sp1 in skyParams:
            for sp2 in skyParams:
                if not (sp1 == sp2):
                    twoDGreedyMenu.append([sp1, sp2])
        for sp1,sp2 in combinations(spinParams,2):
          twoDGreedyMenu.append([sp1, sp2])
        for dc1,dc2 in combinations(bppu.tigerParams,2):
            twoDGreedyMenu.append([dc1,dc2])
        for mp in massParams:
             for tp in tidalParams:
                 if not (mp == tp):
                     twoDGreedyMenu.append([mp, tp])
        for sp1,sp2 in combinations(snrParams,2):
                twoDGreedyMenu.append([sp1,sp2])
        twoDGreedyMenu.append(['lambda1','lambda2'])
        twoDGreedyMenu.append(['lam_tilde','dlam_tilde'])
        twoDGreedyMenu.append(['lambdat','dlambdat'])
        for psip in polParams:
            for phip in phaseParams:
                twoDGreedyMenu.append([psip,phip])
            for sp in skyParams:
                twoDGreedyMenu.append([psip,sp])
            for sp in spinParams:
                twoDGreedyMenu.append([psip,sp])

        for i in calibParams[3:]:
          twoDGreedyMenu.append([i,'distance'])

    #twoDGreedyMenu=[['mc','eta'],['mchirp','eta'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['dist','m1'],['ra','dec']]
    #Bin size/resolution for binning. Need to match (converted) column names.
    greedyBinSizes=bppu.greedyBinSizes


    if opts.plot_2d:
        for dt1,dt2 in combinations(['h1_end_time','l1_end_time','v1_end_time'],2):
          twoDGreedyMenu.append([dt1,dt2])
        for dt1,dt2 in combinations( ['h1l1_delay','l1v1_delay','h1v1_delay'],2):
          twoDGreedyMenu.append([dt1,dt2])

    if opts.deltaLogL and not opts.deltaLogP:
        print("DEPRECATION WARNING: --deltaLogL has been replaced by --deltaLogP.  Using the posterior to define burnin criteria")
        deltaLogP = opts.deltaLogL
    else:
        deltaLogP = opts.deltaLogP

    confidenceLevels=bppu.confidenceLevels
    #2D plots list
    #twoDplots=[['mc','eta'],['mchirp','eta'],['mc', 'time'],['mchirp', 'time'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['RA','dec'],['ra', 'dec'],['m1','dist'],['m2','dist'],['mc', 'dist'],['psi','iota'],['psi','distance'],['psi','dist'],['psi','phi0'], ['a1', 'a2'], ['a1', 'iota'], ['a2', 'iota'],['eta','time'],['ra','iota'],['dec','iota'],['chi','iota'],['chi','mchirp'],['chi','eta'],['chi','distance'],['chi','ra'],['chi','dec'],['chi','psi']]
    twoDplots=twoDGreedyMenu
    cbcBayesPostProc(
                        opts.outpath,datafiles,oneDMenus,twoDGreedyMenu,
                        greedyBinSizes,confidenceLevels,twoDplots,
                        #optional
                        injfile=opts.injfile,eventnum=opts.eventnum,
                        trigfile=opts.trigfile,trignum=opts.trignum,
                        skyres=opts.skyres,
                        # direct integration evidence
                        dievidence=opts.dievidence,boxing=opts.boxing,difactor=opts.difactor,
                        # Ellipitical evidence
                        ellevidence=opts.ellevidence,
                        #manual bayes factor entry
                        bayesfactornoise=opts.bsn,bayesfactorcoherent=opts.bci,
                        #manual input for SNR in the IFOs, optional.
                        snrfactor=opts.snr,
                        #nested sampling options
                        ns_flag=opts.ns,ns_Nlive=opts.Nlive,
                        #spinspiral/mcmc options
                        ss_flag=opts.ss,ss_spin_flag=opts.spin,
                        #LALInferenceMCMC options
                        li_flag=opts.lalinfmcmc,deltaLogP=deltaLogP,fixedBurnins=fixedBurnins,nDownsample=opts.downsample,oldMassConvention=opts.oldMassConvention,
                        #followupMCMC options
                        fm_flag=opts.fm,
                        #injected spin frame
                        inj_spin_frame=opts.inj_spin_frame,
                        # Turn of ACF?
                        noacf=opts.noacf,
                        #Turn on 2D kdes
                        twodkdeplots=opts.twodkdeplots,
                        #Turn on R convergence tests
                        RconvergenceTests=opts.RconvergenceTests,
                        # Also save PDFs?
                        savepdfs=opts.nopdfs,
                        #List of covariance matrix csv files used as analytic likelihood
                        covarianceMatrices=opts.covarianceMatrices,
                        #List of meanVector csv files used, one csv file for each covariance matrix
                        meanVectors=opts.meanVectors,
                        #header file for parameter names in posterior samples
                        header=opts.header,
                        # ascii files (one per IFO) containing  freq - PSD columns
                        psd_files=opts.psdfiles,
                        greedy=not(opts.kdecredibleregions)
                    )

    if opts.archive is not None:
        import subprocess
        subprocess.call(["tar","cvzf",opts.archive,opts.outpath])

    # Send an email, useful for keeping track of dozens of jobs!
    # Will only work if the local host runs a mail daemon
    # that can send mail to the internet
    if opts.email:
        try:
            email_notify(opts.email,opts.outpath)
        except Exception,e:
            print 'Unable to send notification email'
            print "The error was %s\n"%str(e)
#
