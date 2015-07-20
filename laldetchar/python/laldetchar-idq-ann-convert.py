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


from optparse import *
import glob
import sys
import os
import pdb
import numpy
from laldetchar.idq import auxmvc_utils
from laldetchar import git_version

__author__ = 'Young-Min Kim <young-min.kim@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

def Rescale_attribute(MVSCTriggers, attr, new_min, new_max):
    variables=list(MVSCTriggers.dtype.names)

    old_max=0
    old_min=0

    for i in range(len(MVSCTriggers)):
        for var in variables:
            if ("_"+attr in var.strip()):
                if MVSCTriggers[var][i] < old_min:
                    old_min = MVSCTriggers[var][i]
                elif MVSCTriggers[var][i] > old_max:
                    old_max = MVSCTriggers[var][i]
    print "old max of %s : %f" % (attr, float(old_max))
    print "old min of %s : %f" % (attr, float(old_min))

    scale_factor = (new_max - new_min) / float(old_max - old_min)

    for i in range(len(MVSCTriggers)):
        for var in variables:
            if ("_"+attr in var.strip()):
                MVSCTriggers[var][i]=(MVSCTriggers[var][i] - old_min)*scale_factor + new_min
    print "%s normalization to (%f,%f) is done" % (attr, new_min, new_max)
        
    return MVSCTriggers


def Rescale_attribute2(MVSCTriggers, attr, new_min, new_max):
    variables=list(MVSCTriggers.dtype.names)

    old_max=0.0
    old_min=0.0

    for i in range(len(MVSCTriggers)):
        for var in variables:
            if ("_"+attr in var.strip()):
                if MVSCTriggers[var][i] < old_min:
                    old_min = MVSCTriggers[var][i]
                elif MVSCTriggers[var][i] > old_max:
                    old_max = MVSCTriggers[var][i]

    print "old max of %s : %f" % (attr, float(old_max))
    print "old min of %s : %f" % (attr, float(old_min))

    old_min=0.0
    if attr == "signif":
        old_max=150000.0
    elif attr == "dt":
        old_max=5.0
        old_min=-5.0
    elif attr == "freq":
        old_max=4096.0
    elif attr == "dur":
        old_max=4.0
    elif attr == "npts":
        old_max=7000.0

    print "set old max of %s : %f" % (attr, float(old_max))
    print "set old min of %s : %f" % (attr, float(old_min))

    scale_factor = (new_max - new_min) / float(old_max - old_min)

    for i in range(len(MVSCTriggers)):
        for var in variables:
            if ("_"+attr in var.strip()):
                MVSCTriggers[var][i]=(MVSCTriggers[var][i] - old_min)*scale_factor + new_min
    print "%s normalization to (%f,%f) is done" % (attr, new_min, new_max)
        
    return MVSCTriggers


def ApplyLog_dt(MVSCTriggers):
    variables=list(MVSCTriggers.dtype.names)
    channels=[]
    for var in variables:
        if '_dt' in var.strip():
            channels.append(var.split('_dt')[0])
    for ch in channels:
        for i in range(len(MVSCTriggers)):
            if ch+'_signif' not in variables:
                if MVSCTriggers[ch+'_dt'][i]:
                    #numpy.put(MVSCTriggers[var],[i],[numpy.sign(MVSCTriggers[ch+'_dt'][i])*(-1.0)*numpy.log(abs(MVSCTriggers[ch+'_dt'][i]))])
                    MVSCTriggers[ch+'_dt'][i]=numpy.sign(MVSCTriggers[ch+'_dt'][i])*(-1.0)*numpy.log(abs(MVSCTriggers[ch+'_dt'][i]))
            elif MVSCTriggers[ch+'_signif'][i] and not MVSCTriggers[ch+'_dt'][i]:
                #numpy.put(MVSCTriggers[var],[i],[numpy.sign(MVSCTriggers[ch+'_dt'][i])*1000.0])
                MVSCTriggers[ch+'_dt'][i]=numpy.sign(MVSCTriggers[ch+'_dt'][i])*1000.0
            else:
                #numpy.put(MVSCTriggers[var],[i],[numpy.sign(MVSCTriggers[ch+'_dt'][i])*(-1.0)*numpy.log(abs(MVSCTriggers[ch+'_dt'][i]))])
                MVSCTriggers[ch+'_dt'][i]=numpy.sign(MVSCTriggers[ch+'_dt'][i])*(-1.0)*numpy.log(abs(MVSCTriggers[ch+'_dt'][i]))
        
    return MVSCTriggers



description = """This program runs converting *.pat file to *.ann file for FANN training and evaluation."""

parser = OptionParser(version='Name: %%prog\n%s'%git_version.verbose_msg, 
                                usage='%prog [options]', 
                                description=description)
parser.add_option("-v","--verbose",action="store_true",help="verbose mode for printing run processing.")
#parser.add_option("-t","--pat-files", default=False, type="string", help="Provide training file names")
parser.add_option("","--transform-dt-function", default="log", type="string", help="Set up a function applying dt. e.g. inverse,log. Default is log.")
parser.add_option("","--normalization-attributes", default="signif,dt,freq,dur,npts", type="string", help="normalization attributes. Default attributes are signif,dt,freq,dur,npts.")
parser.add_option("","--min", default=0.0, type="float", help="new minimum value. Default is 0.0.")
parser.add_option("","--max", default=1.0, type="float", help="new maximum value. Default is 1.0.")

(opts,files)=parser.parse_args()

################ MAIN ##############################


#pat_files=glob.glob(opts.pat_files)
#pat_files = opts.pat_files.split(",")
channel_attributes = opts.normalization_attributes.split(",")
if opts.verbose:
    print "Converting files to ann type:"
    print files
    print "Normalization attributes:"
    print channel_attributes

for pat in files:
    Triggers=auxmvc_utils.ReadMVSCTriggers([pat],Classified=False)
    for attr in channel_attributes:
        if attr == 'dt':
            if opts.transform_dt_function == 'log':
                Triggers = Rescale_attribute2(Triggers,attr,-1.0,1.0)
                Triggers = ApplyLog_dt(Triggers)
                Triggers = Rescale_attribute(Triggers,attr,opts.min,opts.max)
                print "dt log transformation is done."
            else:
                Triggers = Rescale_attribute2(Triggers,attr,opts.min,opts.max)
        else:
            Triggers = Rescale_attribute2(Triggers,attr,opts.min,opts.max)

    file=open(pat.split("/")[-1].replace(".pat",".ann"),"w")
    print "Writing "+pat.split("/")[-1].replace(".pat",".ann")
    variables = list(Triggers.dtype.names)
    for va in ['GPS_s','GPS_ms','signif','SNR','unclean','glitch-rank','i']:
        if va in variables:
            variables.remove(va)
            if opts.verbose:
                print "removed: %s" % va
    if opts.verbose:
        print "the number of variables:%i" % len(variables)
    file.write(str(len(Triggers))+" "+str(len(variables))+" 1"+"\n")
    for i in range(len(Triggers)):
        file.write(" ".join([str(var) for var in list(Triggers[i])[5:-1]])+"\n")
        file.write(str(Triggers[i][-1])+"\n")
    file.close()
    if opts.verbose:
        print "Attributes normalization in "+pat.split("/")[-1]+" is done and re-wrote in "+pat.split("/")[-1].replace(".pat",".ann")
