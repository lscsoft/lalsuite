#
# Copyright (C) 2004, 2005 Cristina V. Torres
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
from optparse import OptionParser
import getopt
import os
import string
import sys
import time

# Simple python code to parse replace into a new DOT file via the dag
# file

#Simple interface where we specify the DAG this DOT was created from"
parser = OptionParser()

parser.add_option("-t","--dot",dest="dotFile",
                  default="",
                  help="This option is used to specify the path to the DOT file we wish to reparse",
                  metavar="FILE.dot")
parser.add_option("-g","--dag",dest="dagFile",
                   default="",
                   help="This options is used to specify the path to the DAG file which the DOT file was created for.",
                   metavar="FILE.dag")
parser.add_option("-o","--out",dest="outFile",
                   default="./out.dot",
                   help="This option specified the complete pathname for the resulting reparsed DOT information.\n Remember to make PS of the DOT file use dot -Tps OUT.dot -o OUT.ps",
                   metavar="OUT.dot")
(options,args)=parser.parse_args()

mydotfile=str(options.dotFile)
mydagfile=str(options.dagFile)
myoutfile=str(options.outFile)
if mydotfile == '' or mydagfile == '':
    print "Error with input file args!."
    os.abort()

dot_fp=open(mydotfile,'r')
dotData=dot_fp.readlines()
dot_fp.close()

dag_fp=open(mydagfile,'r')
dagData=dag_fp.readlines()
dag_fp.close()

#Scan the DAG file to make serial num and submit file key pairs
dagJobs=[]
color="white"
for entry in dagData:
    if entry.__contains__('JOB '):
        tmpVar=entry.split(' ')
        dagJobs.append([tmpVar[1],os.path.basename(tmpVar[2]).replace('\n',''),str(color)])

#Now we scan the DOT file we must quote all strings which are not
#labels
newDotData=[]
for entry in dotData:
    for key in dagJobs:
        editLine=entry
        if editLine.__contains__(key[0]):
            entry=str(editLine)
            if str(entry).lower().__contains__('(done)'):
                color="green"
            elif str(entry).lower().__contains__('(r)'):
                color="yellow"
            elif str(entry).lower().__contains__('(i)'):
                color="blue"
            else:
                color="red"
            editLine2=editLine.replace('="'+key[0],'="'+key[1])
            editLine3=editLine2.replace(key[0],'"'+key[0]+'"')
            editLine4=editLine3.replace("]",' style=filled fillcolor='+color+']')
            entry=editLine4
    newDotData.append(entry.replace('\n',''))
            
newDot_fp=open(myoutfile,'w')
for entry in newDotData:
    newDot_fp.write(entry+'\n')
newDot_fp.close()    
        
 
