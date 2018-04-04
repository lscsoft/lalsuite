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

def makeInjObjs(injfile,event,posfiles):
    """
    Make a list of results objects from the posfiles and injections
    """
    getByTime=False
    resObjs=[]
    if injfile:
        import itertools
        injections = SimInspiralUtils.ReadSimInspiralFromFiles([injfile])
    i=0
    for posfile in posfiles:
        peparser=bppu.PEOutputParser('common')
        resObj=peparser.parse(open(posfile,'r'))
        pos=bppu.Posterior(resObj)
        if event is not None:
            injection=injections[event]
        else:
            time=pos['time'].mean
            injection=getInjByTime(time,injections)
        if injection is None:
            continue
        pos.set_injection(injection)
        resObjs.append(pos)
        i=i+1
    return resObjs

def makeSummaryFile(obj, params, outpath, confidencelevels,skyres=0.5):
    """
    Make a summary page with table of results and plots for each param in params
    """
    #Bin size/resolution for binning. Need to match (converted) column names.
    GreedyRes={'mc':0.0001,'m1':0.1,'m2':0.1,'mass1':0.1,'mass2':0.1,'mtotal':0.1,'eta':0.001,'iota':0.05,'time':5e-4,'distance':3.0,'dist':3.0,'mchirp':0.0001,'a1':0.02,'a2':0.02,'phi1':0.05,'phi2':0.05,'theta1':0.05,'theta2':0.05,'ra':0.01,'dec':0.01,'psi':0.01,'polarization':0.01}

    if 'distance' in obj.names:
        dist_name = 'distance'
    elif 'dist' in obj.names:
        dist_name = 'dist'

    for par in params:
        par=par.lower()
        print 'Binning %s to determine confidence intervals'%(par)
        
        try:
            obj[par]
        except KeyError:
            print 'No input chain for %s, skipping'%(par)
            continue
        try:
            GreedyRes[par]
        except KeyError:
            print 'No bin size for %s, skipping'%(par)
            continue
        binParams={par: GreedyRes[par]}

        toppoints,injectionconfidence,reses,injection_area,cl_intervals=bppu.greedy_bin_one_param(obj,binParams, confidencelevels)
        
        statfile=open(os.path.join(outpath,par+'_int.txt'),'w')
        for level in confidencelevels:
            print >>statfile,'%lf %lf'%(level, reses[level])
        if injectionconfidence is not None and injection_area is not None:
            print >>statfile,'%lf %lf'%(injectionconfidence,injection_area)
        else:
            print >>statfile,'0 0\n'
        print >>statfile,'12345 %lf'%(obj[par].stdev)
        print >>statfile,'67890 %lf'%(obj[par].stacc) 
        statfile.close()
    
    # Sky position
    top_ranked_pixels,sky_inj_cl,skyreses,injection_area=bppu.greedy_bin_sky(obj,skyres,confidencelevels)
    print "BCI for sky area:"
    print skyreses
    statfile=open(os.path.join(outpath,'sky_int.txt'),'w')
    fracs=sorted(skyreses.keys())
    skysizes=[skyreses[frac] for frac in fracs]
    for frac in fracs:
        print >>statfile,'%lf %lf'%(frac,skyreses[frac])
    if sky_inj_cl is not None and injection_area is not None:
        print >>statfile,'%lf %lf'%(sky_inj_cl,injection_area)
    else:
        print >>statfile,'0 0'
    statfile.close()

    # distance-iota
    greedy2params={dist_name:GreedyRes[dist_name], 'iota':GreedyRes['iota']}
    statfile=open(os.path.join(outpath,'dist_iota_int.txt'),'w')
    toppoints,injection_cl,reses,injection_area=bppu.greedy_bin_two_param(obj,greedy2params,confidencelevels)
    for frac in sorted(reses.keys()):
        print >>statfile,'%lf %lf'%(frac,reses[frac])
    if injection_cl is not None and injection_area is not None:
        print >>statfile,'%lf %lf'%(injection_cl,injection_area)
    else:
        print >>statfile,'0 0'


if __name__=='__main__':
    from optparse import OptionParser
    parser=OptionParser()
    parser.add_option("-i","--inj",help="Injection XML",metavar="INJ.XML",default=None,type="string")
    parser.add_option("-p","--parameter",help="Name of parameter to analyse, can be specified multiple times",metavar="mchirp",action="append",default=[],type="string")
    parser.add_option("-o","--outdir",help="output directory",metavar="PATH",type="string",default="./")
    parser.add_option("-e","--event",help="Injection event number",metavar="NUM",type="int",default=None,action="store")

    (opts,args)=parser.parse_args()
    if not os.path.isdir(opts.outdir):
        os.makedirs(opts.outdir)

    
    confidencelevels=[0.65,0.9,0.95]
    results_objects=makeInjObjs(opts.inj,opts.event,args)
    makeSummaryFile(results_objects[0],opts.parameter,opts.outdir,confidencelevels)

