#!/usr/bin/env python
#
# Copyright (C) 2009 Cristina Valeria Torres
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
This script queries the available Veto information given a gps time with
an asymetric windows specified at command line. This routine will
return a text string in MoinMoin as a table for inclusion in the
candidate checklist.
"""
__author__  = "Cristina Valeria Torres <cristina.torres@ligo.org>"
__prog__    = 'followupRatioTest.py'


import optparse
import sys
import os
from pylal import git_version
from numpy import any
from pylal.fu_utils import ratioTest

sys.path.append('@PYTHONLIBDIR@')

usage = """usage: %prog [options]"""

parser = optparse.OptionParser(usage,version=git_version.verbose_msg)
#Add all options to setup the query
parser.add_option("-R","--snr-ratio-test",action="store",type="string",\
                  metavar="PATH2PICKLE",default=None,
                  help="Set the location of the data (pickle) file\
used to perform the ratio check on the candidate file.")
parser.add_option("-i","--ifo1",action="store",type="string",\
                  metavar="IFO1",default=None,
                  help="Specify the name of the first IFO as\
L1,H1,H2,V1,G1, or T1")
parser.add_option("-j","--ifo2",action="store",type="string",\
                  metavar="IFO2",default=None,
                  help="Specify the name of the first IFO as\
L1,H1,H2,V1,G1, or T1")
parser.add_option("-k","--ifo3",action="store",type="string",\
                  metavar="IFO2",default=None,
                  help="Specify the name of the first IFO as\
L1,H1,H2,V1,G1, or T1")
                  
parser.add_option("-I","--snr1",action="store",type="string",\
                  metavar="SNRatIFO2",default=None,
                  help="Specify the SNR value for the IFO specified by\
 --IFO1.")
parser.add_option("-J","--snr2",action="store",type="string",\
                  metavar="SNRatIFO2",default=None,
                  help="Specify the SNR value for the IFO specified by\
 --IFO1.")
parser.add_option("-K","--snr3",action="store",type="string",\
                  metavar="SNRatIFO3",default=None,
                  help="Specify the SNR value for the IFO specified by\
 --IFO1.")
parser.add_option("-A","--time1",action="store",type="string",\
                  metavar="SNRatIFO3",default=None,
                  help="Specify the arrival time value for the IFO specified by\
 --IFO1.")
parser.add_option("-B","--time2",action="store",type="string",\
                  metavar="SNRatIFO3",default=None,
                  help="Specify the arrival time value for the IFO specified by\
 --IFO2.")
parser.add_option("-C","--time3",action="store",type="string",\
                  metavar="SNRatIFO3",default=None,
                  help="Specify the arrival value for the IFO specified by\
 --IFO3.")
parser.add_option("-p","--output-format",action="store",type="string",\
                      metavar="OUTPUT_TYPE",default="MOINMOIN", \
                      help="The valid output options here are \
LIST(python var), MOINMOIN, and HTML.")
parser.add_option("-o","--output-file",action="store",type="string",\
                  metavar="OUTPUTFILE.wiki",default=None,\
                  help="This sets the name of the file to save the DQ\
 result into.")
                  
######################################################################

(opts,args) = parser.parse_args()


pickleFile=opts.snr_ratio_test
ifo1=opts.ifo1
ifo2=opts.ifo2
ifo3=opts.ifo3
snr1=opts.snr1
snr2=opts.snr2
snr3=opts.snr3
time1=opts.time1
time2=opts.time2
time3=opts.time3
outputFile=opts.output_file
outputType=opts.output_format
#
test=ratioTest()
test.setPickleLocation(pickleFile)
#Create ifo,snr,time listing if the opts are not(None)
ifoPairs=list()
if not any([((x==None) or (x.lower()=='none')) for x in [ifo1,snr1,time1]]):
    ifoPairs.append([ifo1,snr1,time1])
if not any([((x==None) or (x.lower()=='none')) for x in [ifo2,snr2,time2]]):
    ifoPairs.append([ifo2,snr2,time2])
if not any([((x==None) or (x.lower()=='none')) for x in [ifo3,snr3,time3]]):    
    ifoPairs.append([ifo3,snr3,time3])
#
if len(ifoPairs)<2:
    raise Exception, "Insufficient IFO input arguments!"

rawResults=test.checkPairs(ifoPairs)
#
if outputType.upper().strip() == "LIST":
    results=rawResults
if outputType.upper().strip() == "MOINMOIN":
    results=test.generateMOINMOINTable(rawResults)
if outputType.upper().strip() == "HTML":
    results=test.generateHTMLTable(rawResults)

if outputFile == None:
    sys.stdout.write("%s"%(results))
else:
    fp=file(outputFile,"w")
    fp.writelines("%s"%(results))
    fp.close()


