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
This script queries the available DQ information given a gps time with
an asymetric windows specified at command line. This routine will
return a text string in MoinMoin as a table for inclusion in the
candidate checklist.
"""
__author__  = "Cristina Valeria Torres <cristina.torres@ligo.org>"
__prog__    = 'followupQueryDQ.py'


import optparse
import sys
import os
from pylal import git_version
from pylal.fu_utils import followupDQV

sys.path.append('@PYTHONLIBDIR@')

usage = """usage: %prog [options]"""

parser = optparse.OptionParser(usage,version=git_version.verbose_msg)
#Add all options to setup the query
parser.add_option("-I","--ifo-list",action="store",type="string",\
                  metavar="IFOLIST",default="H1,H2,V1,L1",\
                  help="Give a comma separated list of ifos such as \
H1,H2,L1,V1.")
parser.add_option("-X","--segment-url",action="store",type="string",\
                      metavar="URL",default=None,\
                      help="Using this argument specify a URL the LDBD \
server that you want to query DQ Veto segment information from for\
example ldbd://metaserver.phy.syr.edu:30015")
parser.add_option("-w","--window",action="store",type="string",\
                      metavar="frontWin,backWin",default="30,15",\
                      help="Using this argument you can specify a \
asymetric window around the trigger time for performing DQ queries. \
The two times should be postive values seperated by a comma. If only \
one time is specified then the window is assumed to be symetric. \
Example: --window='100,10'")
parser.add_option("-t","--trigger-time",action="store",type="string", \
                      metavar="GPStime", default=None,\
                      help="Using this argument you can specify the \
GPS time of the trigger to check the data quality flags on.")
parser.add_option("-p","--output-format",action="store",type="string",\
                      metavar="OUTPUT_TYPE",default="MOINMOIN", \
                      help="The valid output options here are \
LIST(python var), MOINMOIN, and HTML.")
parser.add_option("-o","--output-file",action="store",type="string",\
                  metavar="OUTPUTFILE.wiki",default=None,\
                  help="This sets the name of the file to save the DQ\
 result into.")
parser.add_option("-b","--estimate-background",action="store_true",\
                  default=False,help="Use this flag to \
invoke the generation of a DQ background to include in the output \
of the DQ information for the time in question.")
parser.add_option("-l","--background-location",action="store",\
                   type="string",metavar="myDiskURL", \
                   default="automatic",help="Specify the disk path to the \
stored location of a static DQ background file.  This is epoch \
specific.")
parser.add_option("-B","--blind",action="store_true",\
                  default=False,help="Set this flag so that your query \
to the segment database is inspiral injection blinded.  i.e. Your \
query will ignore know injection flags related to inspiral search.")

######################################################################

(opts,args) = parser.parse_args()


server=opts.segment_url
triggerTime=opts.trigger_time
outputType=opts.output_format
outputFile=opts.output_file
ifos=opts.ifo_list.upper().split(",")
backgroundLocation=str(opts.background_location).strip()
estimateBackground=opts.estimate_background

#If ifo args seem wrong
if sum([len(x) != 2 for x in ifos]):
    sys.stderr.write("The args passed to --ifo-list are incorrectly formatted! %s\n"%opts.ifo_list)
    sys.exit(1)
    
if len(opts.window.split(',')) == 1:
    frontWindow=backWindow=opts.window
if len(opts.window.split(',')) == 2:
    frontWindow,backWindow=opts.window.split(',')
    if len(backWindow) == 0:
        backWindow=frontWindow
x=followupDQV(server,blinded=opts.blind)
x.fetchInformationDualWindow(triggerTime,frontWindow,backWindow,ifoList=ifos)
if estimateBackground:
    cancelBackgroundEstimation=False
    if backgroundLocation == "automatic":
        backgroundLocation=x.figure_out_pickle("automatic")
        x.resetPicklePointer(backgroundLocation)
        #Is automatically determined pickle present?
        if not os.path.isfile(x.getPicklePointer()):
            cancelBackgroundEstimation=True
    elif backgroundLocation != "":
        x.resetPicklePointer(backgroundLocation)
        #Check background file exists!
        if not os.path.isfile(x.getPicklePointer()):
            cancelBackgroundEstimation=True
    if not cancelBackgroundEstimation:
        x.estimateDQbackground()
    else:
        sys.stderr.write("%s does not exist!\n"%(backgroundLocation))
        sys.stderr.write("Generate one with followupGenerateDQBackground.py.\n")
        sys.stderr.write("Skipping background use...\n")        
result=""
if outputType.upper().strip() == "LIST":
    result=x.generateResultList()
if outputType.upper().strip() == "MOINMOIN":
    result=x.generateMOINMOINTable("DQ")
if outputType.upper().strip() == "HTML":
    result=x.generateHTMLTable("DQ")

if outputFile == None:
    sys.stdout.write("%s"%(result))
else:
    fp=file(outputFile,"w")
    fp.writelines("%s"%(result))
    fp.close()
