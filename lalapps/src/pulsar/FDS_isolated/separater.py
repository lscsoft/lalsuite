#!/usr/bin/env python


## From each file, take contents between textfrom and textto.
textfrom="%1"
textto="%coincidences"


## Import standard modules to the python path
import string, os,sys, shutil, time, os.path
from optparse import OptionParser


## command line options
parser = OptionParser()
parser.add_option("-w", "--workdir", dest="workdir",
                  help="Directory where the result directory will be.",
		  default="/home/yousuke/workspace/S3EatH/worskapce/")
parser.add_option("-t", "--targetdir", dest="targetdir",
                  help="Directory under which the data files will be.",
		  default="/home/yousuke/workspace/S3EatH/S3/")
parser.add_option("-r", "--resultdir", dest="resultdir",
                  help="Directory where the result files will be "
		  "under workdir.",
		  default="datadir/")
(options, args) = parser.parse_args()


workdir=options.workdir
targetdir=options.targetdir
resultdir=workdir+options.resultdir
os.mkdir(resultdir)


##---------------------------------------------------------------
## Functions
##---------------------------------------------------------------

def addheaderinfo(targetfile="",info="",headerlabel="%"):
   """
	add info to the targetfile header. The file header should be 
	indicated by headerlabel. 
   """	
   __errmsg="Error"
   if (targetfile is ""):
      print __errmsg, "Enough numbers of arguments are not given." 
      sys.exit()    
   if os.path.isfile(targetfile) is False:
      print __errmsg, targetfile, "does not exit!!" 
      sys.exit() 
##      
   __tmpfilename=targetfile+"_tmp_addfilename"
   shutil.move(targetfile,__tmpfilename)
   __command=["sed -e s/",headerlabel,"/",headerlabel,
              info, "___/ < ", __tmpfilename,
              " > ", targetfile]
   __comexe=string.join(__command,"")
   os.system(__comexe)
   os.remove(__tmpfilename)



def fileseparater(targetfile="",textto="",resultdir="",niter=2,sectionlabel="%"):
   """
      split the sections in the targetfile into files and put them on 
      resultdir. The sections are defined by sectionlabel.
   """
   __errmsg="Error"
   if (targetfile is "" or textto is "" or resultdir is ""):
      print __errmsg, "Enough numbers of arguments are not given." 
      sys.exit()    
   if os.path.isfile(targetfile) is False:
      print __errmsg, targetfile, "does not exit!!" 
      sys.exit() 
##
   __tmpfilename=targetfile+"_tmp_fileseparater"
   shutil.move(targetfile,__tmpfilename)
   for i in range(niter-1):
       __textfrom=str(i+1)
       __textto="%"+str(i+2)
       __resultfile=targetfile+"___"+__textfrom
       __command=["sed -e /",
		  sectionlabel,     
		  __textfrom,
		  "/,/",
		  __textto,
		  "/p -e d ", 
		  __tmpfilename,
		  " | grep -v ", 
		  __textto, 
		  " > ",
		  __resultfile]    
       __comexe=string.join(__command,"")
       os.system(__comexe)
       shutil.move(__resultfile,resultdir)
   __textfrom=str(niter)
   __resultfile=targetfile+"___"+__textfrom
   __command=["sed -e /",
              sectionlabel,
              __textfrom,
	      "/,/",
              textto,
	      "/p -e d ", 
              __tmpfilename,
              " | grep -v ", 
              textto, 
              " > ",
              __resultfile]    
   __comexe=string.join(__command,"")
   os.system(__comexe)
   shutil.move(__resultfile,resultdir)
   os.remove(__tmpfilename)
 



##---------------------------------------------------------------
## Ask intention
if os.path.isdir(resultdir):		
   print "Current size of the result directory."
   os.system("du -h "+resultdir)
print "size of the target directory before unzip"
os.system("du -h "+targetdir+" | tail -n 1")
uinput=raw_input("Are you sure to proceed?:[y/n]")
if uinput is "y":
   pass
else:
   sys.exit()		   


##---------------------------------------------------------------
## Main 

## This loop separates files.
for parentdir, childdirs, files in os.walk(targetdir):
    for file in files:              # For each file under the targetdir,
	if file is not []:          # If there is a file, 
	   copiedfile=workdir+file  
	   targetfile=os.path.join(parentdir,file)   
	   shutil.copy2(targetfile,copiedfile) # Copy file to working dir.
	   zipfilename=copiedfile+"_tmp.zip"   
	   shutil.move(copiedfile,zipfilename)
	   os.system("unzip -q "+zipfilename)   # Unzip the copied file
	   os.remove(zipfilename)	
	   fileseparater(targetfile=copiedfile,textto=textto,resultdir=resultdir)
## This loop add filenames to files header.
for parentdir, childdirs, files in os.walk(resultdir):
    for file in files:              # For each file under the targetdir,
	if file is not []:          # If there is a file, 
	   targetfile=os.path.join(parentdir,file)   
	   addheaderinfo(targetfile=targetfile,info=file)
