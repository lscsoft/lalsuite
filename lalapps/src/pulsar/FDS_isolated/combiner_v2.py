#!/usr/bin/env python


## combiner_v2.py based on combiner.py
## -adapted to new result file format, incl. f1dot values
## -produces input file for zellepolka from *ZIPPED* EaH result files
## -optional inclusion of each result filename to each row of data
##  for later frequency shifting to fixed fiducial time
 
## Import standard modules to the python path
import string, os,sys, shutil, time, os.path, zipfile
from optparse import OptionParser


## command line options
parser = OptionParser()
parser.add_option("-w", "--workdir", dest="workdir",
                  help="Directory where the result file will be.",
		  default="./")
parser.add_option("-t", "--targetdir", dest="targetdir",
                  help="Directory under which the data files are.",
                  default="./targetdir/")
parser.add_option("-i", "--tmpdir", dest="tmpdir",
                  help="Directory where the result files will be "
		  "under workdir.",
		  default="separate_tmpdir/")
parser.add_option("-r", "--resultfile", dest="resultfile",
                  help="Name of the result file.",
		  default="Combined")
parser.add_option("-n", "--includefnames", dest="includefnames",
                  help="Ask your intention whether to include each result filename.",
		  default=1)
parser.add_option("-d", "--datastretchID", dest="datastretchID",
                  help="Give data-stretch ID of files to process.",
                  default="5515")
parser.add_option("-c", "--candThrPerData", dest="candThrPerData",
                  help="Give candidate threshold per data-stretch.",
                  default="24923")
parser.add_option("-v", "--verbose", dest="vrbflg",
                  help="Ask your intention for safety.",
		  default=1)
(options, args) = parser.parse_args()

datastretchID=options.datastretchID
candThrPerData=options.candThrPerData
print "Number of candidates kept per data-stretch: "+candThrPerData
inclfnames=options.includefnames
vrbflg=options.vrbflg

if inclfnames == 1:
    print "'Einstein at Home' - result filenames will be included into combined file."

## set-up and normalize directory names.
workdir=os.path.abspath(options.workdir)+"/"
targetdir=os.path.abspath(options.targetdir)+"/"
tmpdir=os.path.abspath(workdir+options.tmpdir)+"/"
resultfile=os.path.abspath(workdir+options.resultfile)
## check if there are such directories.
if os.path.isdir(workdir) is False:
    print "Directory ",workdir ," does not exist."
if os.path.isdir(targetdir) is False:
    print "Directory ",targetdir ," does not exist."

## print out the directories.
if vrbflg == 1:
    print "working directory: "+workdir
    print "target directory: "+targetdir
    print "result file name: "+resultfile

try: os.mkdir(tmpdir)
except OSError, err:
    import errno
    print "Warning:", err
    uinput=raw_input("Are you sue to delete "+tmpdir+".: [y/n]")
    if uinput is "y":
       os.system("rm -rf "+tmpdir)
       os.mkdir(tmpdir)
    else:
       sys.exit()		   


##---------------------------------------------------------------
## Functions
##---------------------------------------------------------------


def addfooterinfo(targetfile="",info=""):
   """
	add info to the targetfile footer. 
   """	
   __funcname="addfooterinfo"  
   __errmsg="Error"+" in "+__funcname
   if (targetfile is ""):
      print __errmsg, "Enough numbers of arguments are not given." 
      sys.exit()    
   if os.path.isfile(targetfile) is False:
      print __errmsg, targetfile, "does not exit!!" 
      sys.exit() 
##      
   __tmpfilename=targetfile+"_tmp_addfilename"
   shutil.move(targetfile,__tmpfilename)
   __command=["sed '$a", info, "' < ", __tmpfilename," > ", targetfile]
   __comexe=string.join(__command,"")
   os.system(__comexe)
   os.remove(__tmpfilename)
   
   

def addinfoToEachRaw(targetfile="",info=""):
    """
       add info to each raw. 
    """
    __funcname="addinfoToEachRaw"  
    __errmsg="Error"+" in "+__funcname
    if (targetfile is ""):
        print __errmsg, "Enough numbers of arguments are not given." 
        sys.exit()    
    if os.path.isfile(targetfile) is False:
        print __errmsg, targetfile, "does not exit!!" 
        sys.exit() 
##      
    __tmpfilename=targetfile+"_tmp_addfilename"
    shutil.move(targetfile,__tmpfilename)
    __command=["sed -e 's/^/",
               info, " /' < ", __tmpfilename,
               " > ", targetfile]
    __comexe=string.join(__command,"")
    os.system(__comexe)
    os.remove(__tmpfilename)
    

def deleteinfo(targetfile="",label=""):
    """
       delete info whose first colomn has a label label.
    """
    __funcname="deleteinfo"  
    __errmsg="Error"+" in "+__funcname
    if (targetfile is ""):
        print __errmsg, "Enough numbers of arguments are not given." 
        sys.exit()    
    if os.path.isfile(targetfile) is False:
        print __errmsg, targetfile, "does not exit!!" 
        sys.exit() 
##      
    __tmpfilename=targetfile+"_tmp_addfilename"
    shutil.move(targetfile,__tmpfilename)
    __command=["sed -e '/^",
               label, "/d' < ", __tmpfilename,
               " > ", targetfile]
    __comexe=string.join(__command,"")
    os.system(__comexe)
    os.remove(__tmpfilename)




def appendfile(fromfile="",tofile=sys.stdout):
    """
       Append fromfile to tofile.
    """
    if fromfile is "":
       print "Error"
       sys.exit()	
    os.system("cat "+fromfile+" >> "+tofile)
    

##---------------------------------------------------------------
## Ask intention
if os.path.isfile(resultfile):		
   print "Current size of the result file."
   os.system("du -hs "+resultfile)
   if vrbflg == 1:
       uinput=raw_input("Do you want to delete [d] the existing result file or append [a] results to it?:[d/a]")
       if uinput is "d":
           os.remove(resultfile)
       elif uinput is "a":
           pass
       else:
           sys.exit()
   
print "size of the target directory before unzip (if it is a zipped file)"
os.system("du -hs "+targetdir+" | tail -n 1")
if vrbflg == 1:
    uinput=raw_input("Are you sure to proceed?:[y/n]")
    if uinput is "y":
        pass
    else:
        sys.exit()		   


##---------------------------------------------------------------
## Main 

##datastretchIDs=[5515, 5538, 5613, 5653, 5783, 5813, 5828, 5946, 5955, 6102, 6120, 6126, 6130, 6341, 6497, 6514, 6537]

## This loop separates files.
for parentdir, childdirs, files in os.walk(targetdir):
    for file in files:              # For each file under the targetdir,
        if file is not []:          # If there is a file,
        ## work only data-stretch after data-stretch ##
            myfp = os.popen("/export1/people/hpletsch/opt/lal/src/lalapps/src/pulsar/FDS_isolated/dumpWUparams -i %s" % file)
            thisdatastretch=myfp.read(4)
            myfp.close()
            if thisdatastretch == datastretchID:
                myfp2 = os.popen("/export1/people/hpletsch/opt/lal/src/lalapps/src/pulsar/FDS_isolated/dumpWUparams2 -i %s -c %s" % (file,str(candThrPerData)) )
                candThrPerWU=myfp2.readline()
                myfp2.close()
                copiedfile=workdir+file  
                targetfile=os.path.join(parentdir,file)   
                shutil.copy2(targetfile,copiedfile) # Copy file to working dir.
                f=open(copiedfile,'r')
                buff=f.read(2)
                f.close()
                print str(candThrPerWU)+" "+copiedfile
                if buff == 'PK':
                    zipfilename=copiedfile+"_tmp.zip"   
                    shutil.move(copiedfile,zipfilename)
                    os.system("unzip -p "+zipfilename+" | sort -n -r +4 | more | head -"+str(candThrPerWU)+" > "+copiedfile)   # Unzip the copied file
                    os.remove(zipfilename)
                    shutil.move(copiedfile, tmpdir)
                else:
                    copiedfile2=copiedfile+"2"
                    os.system("more "+copiedfile+" | sort -n -r +4 | more | head -"+str(candThrPerWU)+" > "+copiedfile2)   # Unzip the copied file
                    shutil.move(copiedfile2, tmpdir+file)
                    pass  # if the file already unzipped, do nothing. 

##This loop add various information to files.

fileid=0
for parentdir, childdirs, files in os.walk(tmpdir):
    for file in files:              # For each file under the targetdir,
        if file is not []:          # If there is a file, 
            targetfile=os.path.join(parentdir,file)
            deleteinfo(targetfile=targetfile,label="%")
            addinfoToEachRaw(targetfile=targetfile,info=str(fileid))
            if inclfnames == 1:
                addinfoToEachRaw(targetfile=targetfile,info=str(file))
            fileid+=1           
            appendfile(fromfile=targetfile,tofile=resultfile)
addfooterinfo(targetfile=resultfile,info="%DONE") 
os.system("rm -rf "+tmpdir)
