#!/usr/bin/python

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
for entry in dagData:
    if entry.__contains__('JOB '):
        tmpVar=entry.split(' ')
        dagJobs.append([tmpVar[1],os.path.basename(tmpVar[2]).replace('\n','')])

#Now we scan the DOT file we must quote all strings which are not
#labels
newDotData=[]
for entry in dotData:
    for key in dagJobs:
        editLine=entry
        if editLine.__contains__(key[0]):
            editLine2=editLine.replace('="'+key[0],'="'+key[1])
            editLine3=editLine2.replace(key[0],'"'+key[0]+'"')
            entry=editLine3
    newDotData.append(entry.replace('\n',''))
            
newDot_fp=open(myoutfile,'w')
for entry in newDotData:
    newDot_fp.write(entry+'\n')
newDot_fp.close()    
        
