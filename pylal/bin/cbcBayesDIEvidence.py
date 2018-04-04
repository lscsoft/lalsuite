#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcBayesDIEvidence.py: compute the direct-integration evidence
#       from a file of posterior samples.  See arXiv:0911.1777 for the
#       algorithm.
#
#       Copyright 2011
#       Will M. Farr <will.farr@ligo.org>
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
#
#

import numpy as np
from pylal import bayespputils as bp
from pylal import git_version
from optparse import OptionParser

__author__="Will M. Farr <will.farr@ligo.org>"
__version__="git id %s"%git_version.id
__date__=git_version.date

if __name__=='__main__':
    parser=OptionParser()
    parser.add_option("--data",dest="data",action="store",help="file of posterior samples",metavar="FILE")
    parser.add_option("--Nboxing",dest="Nboxing",action="store",default=64,type="int",metavar="N")
    parser.add_option("--bootstrap", dest='bootstrap',action='store',default=1,type='int',metavar='N')
    parser.add_option('--output',dest='output',action='store',default=None,metavar='FILE')

    (opts,args)=parser.parse_args()

    pos_parser=bp.PEOutputParser('common')

    f=open(opts.data, "r")
    try:
        pos=bp.Posterior(pos_parser.parse(f))
    finally:
        f.close()

    outfile=None
    if opts.output:
        outfile=open(opts.output,'w')
    try:
        for i in range(opts.bootstrap):
            if i == 0:
                log_ev=pos.di_evidence(boxing=opts.Nboxing)
            else:
                log_ev=pos.bootstrap().di_evidence(boxing=opts.Nboxing)
            print 'log(evidence) with Nboxing = %d is %.1f (evidence is %g)'%(opts.Nboxing,log_ev,np.exp(ev))
            if outfile:
                outfile.write('%g\n'%log_ev)
    finally:
        if outfile:
            outfile.close()
