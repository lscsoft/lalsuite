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
                  help="Directory where the result file will be.",
		  default="/home/yousuke/workspace/S3EatH/worskapce/")
parser.add_option("-t", "--targetdir", dest="targetdir",
                  help="Directory under which the data files will be.",
		  default="/home/yousuke/workspace/S3EatH/S3/")
parser.add_option("-r", "--resultfile", dest="resultfile",
                  help="Name of the result file.",
		  default="Combined")
(options, args) = parser.parse_args()


workdir=options.workdir
targetdir=options.targetdir
resultfile=workdir+options.resultfile



##---------------------------------------------------------------
## Functions
##---------------------------------------------------------------

def addfilename(targetfile="",niter=2):
    """
      Replace the line starting as %n (n=1,2,...,niter) with %filename_n
    """
    __tmpfilename=targetfile+"_tmp_addfilename"
    for i in range(niter):
        shutil.move(targetfile,__tmpfilename)
        label=str(i+1)
	targetbase=os.path.basename(targetfile)
        __command=["sed -e s/%",label,"/%",
		   targetbase, "_",label,"/ < ", __tmpfilename,
                   " > ", targetfile]
        __comexe=string.join(__command,"")
        os.system(__comexe)
	os.remove(__tmpfilename)


def datamodifier(targetfile="",textfrom="",textto=""):
   """
      From the targetfile, take contents between textfrom and textto.
   """
   __errmsg="Error"
   if (targetfile is "" or textfrom is "" or textto is ""):
      print __errmsg, "Enough numbers of arguments are not given." 
      sys.exit()    
   if os.path.isfile(targetfile) is False:
      print __errmsg, targetfile, "does not exit!!" 
      sys.exit() 
   __tmpfilename=targetfile+"_tmp_datamodifier"
   shutil.move(targetfile,__tmpfilename)
   __command=["sed -e /",
              textfrom,
	      "/,/",
              textto,
	      "/p -e d ", 
              __tmpfilename,
              " | grep -v ", 
              textto, 
              " > ",
              targetfile]
   __comexe=string.join(__command,"")
   os.system(__comexe)
   os.remove(__tmpfilename)
 


def filecombiner(fromfile="",tofile=sys.stdout):
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
   os.system("du -h "+resultfile)
print "size of the target file before unzip"
os.system("du -h "+targetdir+" | tail -n 1")
uinput=raw_input("Are you sure to proceed?:[y/n]")
if uinput is "y":
   pass
else:
   sys.exit()		   


##---------------------------------------------------------------
## Main 

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
	   datamodifier(targetfile=copiedfile,textfrom=textfrom,textto=textto)
	   addfilename(targetfile=copiedfile)
	   filecombiner(fromfile=copiedfile,tofile=resultfile)
	   os.remove(copiedfile)