#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcBayesBurstPostProc.py
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
from numpy import array,exp,cos,sin,arcsin,arccos,sqrt,size,mean,column_stack,cov,unique,hsplit,correlate,log,dot,power,squeeze,sort
from scipy import stats

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

#local application/library specific imports
from pylal import SimInspiralUtils
from pylal import bayespputils as bppu
from pylal import SimBurstUtils

from pylal import git_version

from glue.ligolw import table
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils
try:
  os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
except:
  os.environ['PATH'] =':/usr/texbin'
print os.environ['PATH']
__author__="Ben Aylott <benjamin.aylott@ligo.org>, Ben Farr <bfarr@u.northwestern.edu>, Will M. Farr <will.farr@ligo.org>, John Veitch <john.veitch@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

def email_notify(address,path):
    import smtplib
    import subprocess
    address=address.split(',')
    SERVER="localhost"
    USER=os.environ['USER']
    HOST=subprocess.check_output(["hostname","-f"]).strip()
    FROM=USER+'@'+HOST
    SUBJECT="LALInference result is ready at "+HOST+"!"
    # Guess the web space path for the clusters
    fslocation=os.path.abspath(path)
    webpath='posplots.html'
    if 'public_html' in fslocation:
        k='public_html'
    elif 'WWW' in fslocation:
        k='WWW'
    else:
        k=None
    if k is not None:
        (a,b)=(fslocation,'')
        while a!=k:
            (a,b)=fslocation.split(a)
            webpath=os.path.join(b,webpath)
    else: webpath=os.path.join(fslocation,'posplots.html')

    if 'atlas.aei.uni-hannover.de' in HOST:
        url="https://atlas1.atlas.aei.uni-hannover.de/"
    elif 'ligo.caltech.edu' in HOST:
        url="https://ldas-jobs.ligo.caltech.edu/"
    elif 'ligo-wa.caltech.edu' in HOST:
        url="https://ldas-jobs.ligo-wa.caltech.edu/"
    elif 'ligo-la.caltech.edu' in HOST:
        url="https://ldas-jobs.ligo-la.caltech.edu/"
    elif 'phys.uwm.edu' in HOST:
        url="https://ldas-jobs.phys.uwm.edu/"
    elif 'phy.syr.edu' in HOST:
        url="https://sugar-jobs.phy.syr.edu/"
    else:
        url=HOST+':'
    url=url+webpath

    TEXT="Hi "+USER+",\nYou have a new parameter estimation result on "+HOST+".\nYou can view the result at "+url
    message="From: %s\nTo: %s\nSubject: %s\n\n%s"%(FROM,', '.join(address),SUBJECT,TEXT)
    server=smtplib.SMTP(SERVER)
    server.sendmail(FROM,address,message)
    server.quit()


class LIGOLWContentHandlerExtractSimInspiralTable(ligolw.LIGOLWContentHandler):
    def __init__(self,document):
      ligolw.LIGOLWContentHandler.__init__(self,document)
      self.tabname=lsctables.SimInspiralTable.tableName
      self.intable=False
      self.tableElementName=''
    def startElement(self,name,attrs):
      if attrs.has_key('Name') and attrs['Name']==self.tabname:
        self.tableElementName=name
        # Got the right table, let's see if it's the right event
        ligolw.LIGOLWContentHandler.startElement(self,name,attrs)
        self.intable=True
      elif self.intable: # We are in the correct table
        ligolw.LIGOLWContentHandler.startElement(self,name,attrs)
    def endElement(self,name):
      if self.intable: ligolw.LIGOLWContentHandler.endElement(self,name)
      if self.intable and name==self.tableElementName: self.intable=False

lsctables.use_in(LIGOLWContentHandlerExtractSimInspiralTable)

class LIGOLWContentHandlerExtractSimBurstTable(ligolw.LIGOLWContentHandler):
    def __init__(self,document):
      ligolw.LIGOLWContentHandler.__init__(self,document)
      self.tabname=lsctables.SimBurstTable.tableName
      self.intable=False
      self.tableElementName=''
    def startElement(self,name,attrs):
      if attrs.has_key('Name') and attrs['Name']==self.tabname:
        self.tableElementName=name
        # Got the right table, let's see if it's the right event
        ligolw.LIGOLWContentHandler.startElement(self,name,attrs)
        self.intable=True
      elif self.intable: # We are in the correct table
        ligolw.LIGOLWContentHandler.startElement(self,name,attrs)
    def endElement(self,name):
      if self.intable: ligolw.LIGOLWContentHandler.endElement(self,name)
      if self.intable and name==self.tableElementName: self.intable=False

lsctables.use_in(LIGOLWContentHandlerExtractSimBurstTable)

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

def cbcBayesBurstPostProc(
                        outdir,data,oneDMenu,twoDGreedyMenu,GreedyRes,
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
                        li_flag=False,deltaLogL=None,fixedBurnins=None,nDownsample=None,oldMassConvention=False,
                        #followupMCMC options
                        fm_flag=False,
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
                        #only save stats about PDF and exit #
                        statsonly=False
                    ):
    """
    This is a demonstration script for using the functionality/data structures
    contained in pylal.bayespputils . It will produce a webpage from a file containing
    posterior samples generated by the parameter estimation codes with 1D/2D plots
    and stats from the marginal posteriors for each parameter/set of parameters.
    """
    got_inspiral_table=0
    got_burst_table=0
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
        commonResultsObj=peparser.parse(data,spin=ss_spin_flag,deltaLogL=deltaLogL)

    elif li_flag:
        peparser=bppu.PEOutputParser('inf_mcmc')
        commonResultsObj=peparser.parse(data,outdir=outdir,deltaLogL=deltaLogL,fixedBurnins=fixedBurnins,nDownsample=nDownsample,oldMassConvention=oldMassConvention)

    elif ss_flag and ns_flag:
        raise RuntimeError("Undefined input format. Choose only one of:")

    elif '.xml' in data[0]:
        peparser=bppu.PEOutputParser('xml')
        commonResultsObj=peparser.parse(data[0])
        thefile=open(data[0],'r')
        votfile=thefile.read()
    elif '.hdf' in data[0] or '.h5' in data[0]:
        peparser = bppu.PEOutputParser('hdf5')
        commonResultsObj = peparser.parse(data[0])
    else:
        peparser=bppu.PEOutputParser('common')
        commonResultsObj=peparser.parse(open(data[0],'r'),info=[header,None])
    
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
        xmldoc = utils.load_filename(injfile,contenthandler=LIGOLWContentHandlerExtractSimBurstTable)
        from glue.ligolw import table
        got_burst_table=1
        try:
            lsctables.use_in(LIGOLWContentHandlerExtractSimBurstTable)
            simtable=table.get_table(xmldoc,lsctables.SimBurstTable.tableName)
        except ValueError:
            lsctables.use_in(LIGOLWContentHandlerExtractSimInspiralTable)
            simtable=table.get_table(xmldoc,lsctables.SimInspiralTable.tableName)
            got_inspiral_table=1
            got_burst_table=0
        
        injection=simtable[eventnum]
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
    if got_inspiral_table==1:
        pos = bppu.Posterior(commonResultsObj,SimInspiralTableEntry=injection,injFref=injFref,SnglInpiralList=triggers,votfile=votfile)
    else:
        pos = bppu.BurstPosterior(commonResultsObj,SimBurstTableEntry=injection,injFref=injFref,SnglBurstList=triggers,votfile=votfile)
    #Create analytic likelihood functions if covariance matrices and mean vectors were given
    analyticLikelihood = None
    if covarianceMatrices and meanVectors:
        analyticLikelihood = bppu.AnalyticLikelihood(covarianceMatrices, meanVectors)

        #Plot only analytic parameters
        oneDMenu = analyticLikelihood.names
        twoDGreedyMenu = []
        for i in range(len(oneDMenu)):
            for j in range(i+1,len(oneDMenu)):
                twoDGreedyMenu.append([oneDMenu[i],oneDMenu[j]])
        twoDplots = twoDGreedyMenu

    if eventnum is None and injfile is not None:
        import itertools
        if got_burst_table==1:
            injections = SimBurstUtils.ReadSimBurstFromFiles([injfile])
        elif got_inspiral_table==1:
            injections = SimInspiralUtils.ReadSimInpiralFromFiles([injfile])

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

        # Compute time delays from sky position
    if ('ra' in pos.names or 'rightascension' in pos.names) \
    and ('declination' in pos.names or 'dec' in pos.names) \
    and 'time' in pos.names:
        from pylal import antenna
        from pylal import xlal,inject
        from pylal.xlal import tools,datatypes
        from pylal import date
        from pylal.date import XLALTimeDelayFromEarthCenter
        from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
        import itertools
        detMap = {'H1': 'LHO_4k', 'H2': 'LHO_2k', 'L1': 'LLO_4k',
                'G1': 'GEO_600', 'V1': 'VIRGO', 'T1': 'TAMA_300'}
        if 'ra' in pos.names:
            ra_name='ra'
        else: ra_name='rightascension'
        if 'dec' in pos.names:
            dec_name='dec'
        else: dec_name='declination'
        ifo_times={}
        my_ifos=['H1','L1','V1']
        for ifo in my_ifos:
            inj_time=None
            if injection:
                inj_time=float(injection.get_end(ifo[0]))
            location=inject.cached_detector[detMap[ifo]].location
            ifo_times[ifo]=array(map(lambda ra,dec,time: array([time[0]+XLALTimeDelayFromEarthCenter(location,ra[0],dec[0],LIGOTimeGPS(float(time[0])))]), pos[ra_name].samples,pos[dec_name].samples,pos['time'].samples))
            loc_end_time=bppu.PosteriorOneDPDF(ifo.lower()+'_end_time',ifo_times[ifo],injected_value=inj_time)
            pos.append(loc_end_time)
        for ifo1 in my_ifos:
            for ifo2 in my_ifos:
                if ifo1==ifo2: continue
                delay_time=ifo_times[ifo2]-ifo_times[ifo1]
                if injection:
                    inj_delay=float(injection.get_end(ifo2[0])-injection.get_end(ifo1[0]))
                else:
                    inj_delay=None
                time_delay=bppu.PosteriorOneDPDF(ifo1.lower()+ifo2.lower()+'_delay',delay_time,inj_delay)
                pos.append(time_delay)

    
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
    
    if bayesfactornoise is not None:
        html_model=html.add_section('Model selection',legend=legend)
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
        
    #Create a section for summary statistics
    tabid='statstable'
    html_stats=html.add_collapse_section('Summary statistics',legend=legend,innertable_id=tabid)    
    html_stats.write(str(pos))
    statfilename=os.path.join(outdir,"summary_statistics.dat")
    statout=open(statfilename,"w")
    statout.write("\tmaxP\tmaxL\tstdev\tmean\tmedian\tstacc\tinjection\tvalue\n")
    
    for statname,statoned_pos in pos:

      statmax_pos,max_i=pos._posMaxL()
      statmaxL=statoned_pos.samples[max_i][0]
      statmax_pos,max_j=pos._posMap()
      statmaxP=statoned_pos.samples[max_j][0]
      statmean=str(statoned_pos.mean)
      statstdev=str(statoned_pos.stdev)
      statmedian=str(squeeze(statoned_pos.median))
      statstacc=str(statoned_pos.stacc)
      statinjval=str(statoned_pos.injval)
      
      statarray=[str(i) for i in [statname,statmaxP,statmaxL,statstdev,statmean,statmedian,statstacc,statinjval]]
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
        wfpointer= bppu.plot_burst_waveform(pos=pos,simburst=injfile,event=eventnum,path=wfdir)
    except:
        wfpointer = None
    wftd=html_wf.insert_td(row,'',label='Waveform',legend=legend)
    wfsection=html.add_section_to_element('Waveforms',wftd)
    if wfpointer is not None:
      wfsection.write('<a href="Waveform/WF_DetFrame.png" target="_blank"><img src="Waveform/WF_DetFrame.png"/></a>')
    else:
      wfsection.write("<b>No Waveform generated!</b>")
      
    wftd=html_wf.insert_td(row,'',label='PSDs',legend=legend)
    wfsection=html.add_section_to_element('PSDs',wftd)
    psd_pointer=None
    if psd_files is not None:
      psd_files=list(psd_files.split(','))
      psddir=os.path.join(outdir,'PSDs')
      if not os.path.isdir(psddir):
        os.makedirs(psddir)
      try:
        psd_pointer=bppu.plot_psd(psd_files,outpath=psddir)    
      except:
        psd_pointer=None
    if psd_pointer:
      wfsection.write('<a href="PSDs/PSD.png" target="_blank"><img src="PSDs/PSD.png"/></a>')
    else:
      wfsection.write("<b>No PSD file found!</b>")
    if statsonly:
      return 0
    #==================================================================#
    #1D posteriors
    #==================================================================#

    #Loop over each parameter and determine the contigious and greedy
    #confidence levels and some statistics.

    #Add section for 1D marginal PDFs and sample plots
    tabid='onedmargtable'
    html_ompdf=html.add_collapse_section('1D marginal posterior PDFs',legend=legend,innertable_id=tabid)
    #Table matter
    if not noacf:
        html_ompdf_write= '<table id="%s"><tr><th>Histogram and Kernel Density Estimate</th><th>Samples used</th><th>Autocorrelation</th></tr>'%tabid
    else:
        html_ompdf_write= '<table id="%s"><tr><th>Histogram and Kernel Density Estimate</th><th>Samples used</th></tr>'%tabid
    #Add section for 1D confidence intervals
    tabid='onedconftable'
    html_ogci=html.add_collapse_section('1D confidence intervals (greedy binning)',legend=legend,innertable_id=tabid)
    html_ogci_write='<table id="%s" border="1"><tr><th/>'%tabid
    confidence_levels.sort()
    for cl in confidence_levels:
        html_ogci_write+='<th>%f</th>'%cl
    if injection:
        html_ogci_write+='<th>Injection Confidence Level</th>'
        html_ogci_write+='<th>Injection Confidence Interval</th>'
    html_ogci_write+='</tr>'

    onepdfdir=os.path.join(outdir,'1Dpdf')
    if not os.path.isdir(onepdfdir):
        os.makedirs(onepdfdir)

    sampsdir=os.path.join(outdir,'1Dsamps')
    if not os.path.isdir(sampsdir):
        os.makedirs(sampsdir)
    Nskip=0
    if 'chain' in pos.names:
        data,header=pos.samples()
        par_index=pos.names.index('cycle')
        chain_index=pos.names.index("chain")
        chains=unique(pos["chain"].samples)
        chainCycles = [sort(data[ data[:,chain_index] == chain, par_index ]) for chain in chains]
        chainNcycles = [cycles[-1]-cycles[0] for cycles in chainCycles]
        chainNskips = [cycles[1] - cycles[0] for cycles in chainCycles]
    elif 'cycle' in pos.names:
        cycles = sort(pos['cycle'].samples)
        Ncycles = cycles[-1]-cycles[0]
        Nskip = cycles[1]-cycles[0]

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

        #print "Binning %s to determine confidence levels ..."%par_name
        binParams={par_name:par_bin}

        toppoints,injectionconfidence,reses,injection_area,cl_intervals=bppu.greedy_bin_one_param(pos,binParams,confidence_levels)

        #oneDContCL,oneDContInj = bppu.contigious_interval_one_param(pos,binParams,confidence_levels)

        #Generate new BCI html table row
        BCItableline='<tr><td>%s</td>'%(par_name)
        cls=reses.keys()
        cls.sort()

        for cl in cls:
            BCItableline+='<td>%f</td>'%reses[cl]

        if injection is not None:
            if injectionconfidence is not None and injection_area is not None:

                BCItableline+='<td>%f</td>'%injectionconfidence
                BCItableline+='<td>%f</td>'%injection_area

            else:
                BCItableline+='<td/>'
                BCItableline+='<td/>'

        BCItableline+='</tr>'

        #Append new table line to section html
        html_ogci_write+=BCItableline

        #Generate 1D histogram/kde plots
        print "Generating 1D plot for %s."%par_name

        #Get analytic description if given
        pdf=cdf=None
        if analyticLikelihood:
            pdf = analyticLikelihood.pdf(par_name)
            cdf = analyticLikelihood.cdf(par_name)

        oneDPDFParams={par_name:50}
        rbins,plotFig=bppu.plot_one_param_pdf(pos,oneDPDFParams,pdf,cdf,plotkde=False)

        figname=par_name+'.png'
        oneDplotPath=os.path.join(onepdfdir,figname)
        plotFig.savefig(oneDplotPath)
        if(savepdfs): plotFig.savefig(os.path.join(onepdfdir,par_name+'.pdf'))
        plt.close(plotFig)

        if rbins:
            print "r of injected value of %s (bins) = %f"%(par_name, rbins)

        ##Produce plot of raw samples
        myfig=plt.figure(figsize=(4,3.5),dpi=200)
        pos_samps=pos[par_name].samples
        if not ("chain" in pos.names):
            # If there is not a parameter named "chain" in the
            # posterior, then just produce a plot of the samples.
            plt.plot(pos_samps,'k.',linewidth=0.0, markeredgewidth=0,figure=myfig)
            maxLen=len(pos_samps)
        else:
            # If there is a parameter named "chain", then produce a
            # plot of the various chains in different colors, with
            # smaller dots.
            data,header=pos.samples()
            par_index=pos.names.index(par_name)
            chain_index=pos.names.index("chain")
            chains=unique(pos["chain"].samples)
            chainData=[data[ data[:,chain_index] == chain, par_index ] for chain in chains]
            chainDataRanges=[range(len(cd)) for cd in chainData]
            maxLen=max([len(cd) for cd in chainData])
            for rng, data in zip(chainDataRanges, chainData):
                plt.plot(rng, data, marker=',',linewidth=0.0, markeredgewidth=0,figure=myfig)
            plt.title("Gelman-Rubin R = %g"%(pos.gelman_rubin(par_name)))

            #dataPairs=[ [rng, data] for (rng,data) in zip(chainDataRanges, chainData)]
            #flattenedData=[ item for pair in dataPairs for item in pair ]
            #maxLen=max([len(data) for data in flattenedData])
            #plt.plot(array(flattenedData),marker=',',linewidth=0.0,figure=myfig)


        injpar=pos[par_name].injval

        if injpar:
            if min(pos_samps)<injpar and max(pos_samps)>injpar:
                plt.axhline(injpar, color='r', linestyle='-.')
        myfig.savefig(os.path.join(sampsdir,figname.replace('.png','_samps.png')))
        if(savepdfs): myfig.savefig(os.path.join(sampsdir,figname.replace('.png','_samps.pdf')))
        plt.close(myfig)
        acfail=0
        if not (noacf):
            acffig=plt.figure(figsize=(4,3.5),dpi=200)
            if not ("chain" in pos.names):
                data=pos_samps[:,0]
                try:
                    (Neff, acl, acf) = bppu.effectiveSampleSize(data, Nskip)
                    lines=plt.plot(acf, 'k,', marker=',',linewidth=0.0, markeredgewidth=0, figure=acffig)
                    # Give ACL info if not already downsampled according to it
                    if nDownsample is None:
                        plt.title('Autocorrelation Function')
                    elif 'cycle' in pos.names:
                        last_color = lines[-1].get_color()
                        plt.axvline(acl/Nskip, linestyle='-.', color=last_color)
                        plt.title('ACL = %i   N = %i'%(acl,Neff))
                except FloatingPointError:
                    # Ignore
                    acfail=1
                    pass
            else:
                try:
                    acls = []
                    Nsamps = 0.0;
                    for rng, data, Nskip, Ncycles in zip(chainDataRanges, chainData, chainNskips, chainNcycles):
                        (Neff, acl, acf) = bppu.effectiveSampleSize(data, Nskip)
                        acls.append(acl)
                        Nsamps += Neff
                        lines=plt.plot(acf,'k,', marker=',',linewidth=0.0, markeredgewidth=0, figure=acffig)
                        # Give ACL info if not already downsampled according to it
                        if nDownsample is not None:
                            last_color = lines[-1].get_color()
                            plt.axvline(acl/Nskip, linestyle='-.', color=last_color)
                    if nDownsample is None:
                        plt.title('Autocorrelation Function')
                    else:
                        plt.title('ACL = %i  N = %i'%(max(acls),Nsamps))
                except FloatingPointError:
                    # Ignore
                    acfail=1
                    pass

            acffig.savefig(os.path.join(sampsdir,figname.replace('.png','_acf.png')))
            if(savepdfs): acffig.savefig(os.path.join(sampsdir,figname.replace('.png','_acf.pdf')))
            plt.close(acffig)

        if not noacf:
	  if not acfail:
	    acfhtml='<td width="30%"><img width="100%" src="1Dsamps/'+figname.replace('.png', '_acf.png')+'"/></td>'
	  else:
	    acfhtml='<td>ACF generation failed!</td>'
          html_ompdf_write+='<tr><td width="30%"><img width="100%" src="1Dpdf/'+figname+'"/></td><td width="30%"><img width="100%" src="1Dsamps/'+figname.replace('.png','_samps.png')+'"/></td>'+acfhtml+'</tr>'
        else:
            html_ompdf_write+='<tr><td width="30%"><img width="100%" src="1Dpdf/'+figname+'"/></td><td width="30%"><img width="100%" src="1Dsamps/'+figname.replace('.png','_samps.png')+'"/></td></tr>'


    html_ompdf_write+='</table>'

    html_ompdf.write(html_ompdf_write)

    html_ogci_write+='</table>'
    html_ogci.write(html_ogci_write)

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

        #print "Binning %s-%s to determine confidence levels ..."%(par1_name,par2_name)
        #Form greedy binning input structure
        greedy2Params={par1_name:par1_bin,par2_name:par2_bin}

        #Greedy bin the posterior samples
        try:
          toppoints,injection_cl,reses,injection_area=\
          bppu.greedy_bin_two_param(pos,greedy2Params,confidence_levels)
        except:
          # Failures may happen since some simburst set injval to nan
          continue

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
        try:
          greedy2ContourPlot=bppu.plot_two_param_kde_greedy_levels({'Result':pos},greedy2Params,[0.67,0.9,0.95],{'Result':'k'})
          greedy2contourpath=os.path.join(greedytwobinsdir,'%s-%s_greedy2contour.png'%(par1_name,par2_name))
          greedy2ContourPlot.savefig(greedy2contourpath)
          if(savepdfs): greedy2ContourPlot.savefig(greedy2contourpath.replace('.png','.pdf'))
          plt.close(greedy2ContourPlot)
        
          greedy2HistFig=bppu.plot_two_param_greedy_bins_hist(pos,greedy2Params,confidence_levels)
          greedy2histpath=os.path.join(greedytwobinsdir,'%s-%s_greedy2.png'%(par1_name,par2_name))
          greedy2HistFig.savefig(greedy2histpath)
          if(savepdfs): greedy2HistFig.savefig(greedy2histpath.replace('.png','.pdf'))
          plt.close(greedy2HistFig)
        except:
          pass
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

                myfig.savefig(twoDKdePath)
                if(savepdfs): myfig.savefig(twoDKdePath.replace('.png','.pdf'))
                plt.close(myfig)

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

    # Save posterior samples too...
    posfilename=os.path.join(outdir,'posterior_samples.dat')
    pos.write_to_file(posfilename)

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
    parser.add_option("-i","--inj",dest="injfile",help="SimBurst injection file",metavar="INJ.XML",default=None)
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

    parser.add_option("--no2D",action="store_true",default=False,help="Skip 2-D plotting.")
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
    parser.add_option("--downsample",action="store",default=None,help="(LALInferenceMCMC) approximate number of samples to record in the posterior",type="int")
    parser.add_option("--deltaLogL",action="store",default=None,help="(LALInferenceMCMC) Difference in logL to use for convergence test.",type="float")
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
    parser.add_option("--stats_only",action="store_true",default=False,dest="stats_only")
    parser.add_option("--archive",action="store",default=None,type="string",metavar="results.tar.gz",help="Create the given tarball with all results")
    parser.add_option("--psdfiles",action="store",default=None,type="string",metavar="H1,L1,V1",help="comma separater list of ASCII files with PSDs, one per IFO")
    (opts,args)=parser.parse_args()

    datafiles=[]
    if args:
      datafiles=datafiles+args
    if opts.data:
      datafiles=datafiles + opts.data
    
    if opts.fixedBurnin:
      fixedBurnins = [int(fixedBurnin) for fixedBurnin in opts.fixedBurnin]
    else:
      fixedBurnins = None
    if opts.archive=='None':
      opts.archive=None
    #List of parameters to plot/bin . Need to match (converted) column names.
    
    polParams=['psi','polarisation','polarization']
    skyParams=['ra','rightascension','declination','dec']
    timeParams=['time']
    ellParams=['alpha','polar_eccentricity','polar_angle']
    burstParams=['frequency','loghrss','quality','hrss','duration']
    phaseParams=['phase','phi_orb']
    #endTimeParams=['l1_end_time','h1_end_time','v1_end_time']
    endTimeParams=[]
    #statsParams=['logprior','logl','deltalogl','deltaloglh1','deltalogll1','deltaloglv1','deltaloglh2','deltaloglg1']
    statsParams=['logl']
    calibParams=['calpha_l1','calpha_h1','calamp_l1','calamp_h1']
    oneDMenu=polParams + skyParams + timeParams + statsParams+burstParams+ellParams+phaseParams+calibParams
   
    ifos_menu=['h1','l1','v1']
    from itertools import combinations
    for ifo1,ifo2 in combinations(ifos_menu,2):
      oneDMenu.append(ifo1+ifo2+'_delay')
    
    #oneDMenu=[]
    twoDGreedyMenu=[]
    if not opts.no2D:
        for b1,b2 in combinations(burstParams,2):
            twoDGreedyMenu.append([b1,b2])
        #for bu in burstParams:
        #   for sp in skyParams:
        #        twoDGreedyMenu.append([bu,sp])
        #for bu in burstParams:
        #    for ti in timeParams:
        #        twoDGreedyMenu.append([bu,ti])

        twoDGreedyMenu.append(['phi_orb','psi'])
        twoDGreedyMenu.append(['alpha','psi'])
        twoDGreedyMenu.append(['phi_orb','alpha'])
        twoDGreedyMenu.append(['loghrss','psi'])
        twoDGreedyMenu.append(['alpha','loghrss'])
        for i in calibParams[2:]:
          twoDGreedyMenu.append([i,'loghrss'])
        twoDGreedyMenu.append(['calpha_h1','calpha_l1'])
        twoDGreedyMenu.append(['calamp_h1','calamp_l1'])

    #twoDGreedyMenu=[['mc','eta'],['mchirp','eta'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['dist','m1'],['ra','dec']]
    #Bin size/resolution for binning. Need to match (converted) column names.
    greedyBinSizes={'time':1e-4,'ra':0.05,'dec':0.05,'polarisation':0.04,'rightascension':0.05,'declination':0.05, 'loghrss':0.01,'frequency':0.5,'quality':0.05,'phase':0.1,'phi_orb':0.1,'psi':0.04,'polarization':0.04,'alpha':0.01,'duration':0.0001,'calamp_l1':0.01,'calamp_h1':0.01,'calpha_h1':0.01,'calpha_l1':0.01,'polar_eccentricity':0.01}
    #for derived_time in ['h1_end_time','l1_end_time','v1_end_time','h1l1_delay','l1v1_delay','h1v1_delay']:
    #    greedyBinSizes[derived_time]=greedyBinSizes['time']
    #if not opts.no2D:
    #    for dt1,dt2 in combinations(['h1_end_time','l1_end_time','v1_end_time'],2):
    #      twoDGreedyMenu.append([dt1,dt2])
    #    for dt1,dt2 in combinations( ['h1l1_delay','l1v1_delay','h1v1_delay'],2):
    #      twoDGreedyMenu.append([dt1,dt2])
  
    #Confidence levels
    for loglname in ['logl','deltalogl','deltaloglh1','deltaloglv1','deltalogll1','logll1','loglh1','loglv1']:
        greedyBinSizes[loglname]=0.1
    confidenceLevels=[0.67,0.9,0.95,0.99]
    #2D plots list
    #twoDplots=[['mc','eta'],['mchirp','eta'],['mc', 'time'],['mchirp', 'time'],['m1','m2'],['mtotal','eta'],['distance','iota'],['dist','iota'],['RA','dec'],['ra', 'dec'],['m1','dist'],['m2','dist'],['mc', 'dist'],['psi','iota'],['psi','distance'],['psi','dist'],['psi','phi0'], ['a1', 'a2'], ['a1', 'iota'], ['a2', 'iota'],['eta','time'],['ra','iota'],['dec','iota'],['chi','iota'],['chi','mchirp'],['chi','eta'],['chi','distance'],['chi','ra'],['chi','dec'],['chi','psi']]
    twoDplots=twoDGreedyMenu
    cbcBayesBurstPostProc(
                        opts.outpath,datafiles,oneDMenu,twoDGreedyMenu,
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
                        li_flag=opts.lalinfmcmc,deltaLogL=opts.deltaLogL,fixedBurnins=fixedBurnins,nDownsample=opts.downsample,oldMassConvention=opts.oldMassConvention,
                        #followupMCMC options
                        fm_flag=opts.fm,
                        # Turn of ACF?
                        noacf=opts.noacf,
                        #Turn on 2D kdes
                        #twodkdeplots=opts.twodkdeplots,
                        twodkdeplots=False,
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
                        statsonly=opts.stats_only
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
        except:
            print 'Unable to send notification email'
#
