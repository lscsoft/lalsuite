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
		  default="./")
parser.add_option("-t", "--targetdir", dest="targetdir",
                  help="Directory under which the data files will be.",
		  default="./targetdir/")
parser.add_option("-r", "--resultfile", dest="resultfile",
                  help="Name of the result file.",
		  default="Combined")
(options, args) = parser.parse_args()



## set-up directories and files
workdir=os.path.abspath(options.workdir)+"/"
targetdir=os.path.abspath(options.targetdir)+"/"
resultfile=os.path.abspath(workdir+options.resultfile)
## check if there are such directories.
if os.path.isdir(workdir) is False:
    print "Directory ",workdir ," does not exist."
if os.path.isdir(targetdir) is False:
    print "Directory ",targetdir ," does not exist."
## print out the directories.
print "working directory: "+workdir
print "target directory: "+targetdir
print "result file name: "+resultfile




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


def textcutter(targetfile="",textfrom="",textto=""):
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
   __tmpfilename=targetfile+"_tmp_textcutter"
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
   __command=["sed -e '$a",
              info, "' < ", __tmpfilename,
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
    __command=["sed -e 's/$/   ",
               info, "/' < ", __tmpfilename,
               " > ", targetfile]
    __comexe=string.join(__command,"")
    os.system(__comexe)
    os.remove(__tmpfilename)

    
##---------------------------------------------------------------
## Ask intention
if os.path.isfile(resultfile):		
   print "Current size of the result file."
   os.system("du -h "+resultfile)
   uinput=raw_input("Do you want to delete [d], or append [a]?:[d/a]")
   if uinput is "d":
       os.system("rm -f "+resultfile)
   elif uinput is "a":
       pass
   else:
       sys.exit()
       

print "size of the target file before unzip"
os.system("du -h "+targetdir+" | tail -n 1")
uinput=raw_input("Are you sure to proceed?:[y/n]")
if uinput is "y":
   pass
else:
   sys.exit()		   


##---------------------------------------------------------------
## Main 
fileid=1
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
           ## take sections from textfrom to textto.
	   textcutter(targetfile=copiedfile,textfrom=textfrom,textto=textto)
           addinfoToEachRaw(targetfile=copiedfile,info=str(fileid))
           fileid+=1
           ## write the file name at each section header.
	   addfilename(targetfile=copiedfile)
           ## append the file into one result file.
	   appendfile(fromfile=copiedfile,tofile=resultfile)
	   os.remove(copiedfile)

