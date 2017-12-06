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
		  default="./")
parser.add_option("-t", "--targetdir", dest="targetdir",
                  help="Directory under which the data files will be.",
                  default="./separate_targetdir/")
parser.add_option("-r", "--resultdir", dest="resultdir",
                  help="Directory where the result files will be "
		  "under workdir.",
		  default="separate_resultdir/")
(options, args) = parser.parse_args()


## set-up and normalize directory names.
workdir=os.path.abspath(options.workdir)+"/"
targetdir=os.path.abspath(options.targetdir)+"/"
resultdir=os.path.abspath(workdir+options.resultdir)+"/"
## check if there are such directories.
if os.path.isdir(workdir) is False:
    print "Directory ",workdir ," does not exist."
if os.path.isdir(targetdir) is False:
    print "Directory ",targetdir ," does not exist."
## print out the directories.
print "working directory: "+workdir
print "target directory: "+targetdir
print "result directory: "+resultdir


try: os.mkdir(resultdir)
except OSError, err:
    import errno
    print "Warning:", err
    uinput=raw_input("Are you sue to delete "+resultdir+".: [y/n]")
    if uinput is "y":
       os.system("rm -rf "+resultdir)
       os.mkdir(resultdir)
    else:
       sys.exit()		   



##---------------------------------------------------------------
## Functions
##---------------------------------------------------------------

def addheaderinfo(targetfile="",info=""):
   """
	add info to the targetfile header. 
   """	
   __funcname="addheaderinfo"  
   __errmsg="Error"+" in "+__funcname
   if (targetfile is ""):
      print __errmsg, "Enough numbers of arguments are not given." 
      sys.exit()    
   if os.path.isfile(targetfile) is False:
      print __errmsg, targetfile, "does not exist!!" 
      sys.exit() 
##      
   __tmpfilename=targetfile+"_tmp_addfilename"
   shutil.move(targetfile,__tmpfilename)
   __command=["sed -e '1i",
              info, "' < ", __tmpfilename,
              " > ", targetfile]
   __comexe=string.join(__command,"")
   os.system(__comexe)
   os.remove(__tmpfilename)



def addheaderinfo2(targetfile="",info="",headerlabel="%"):
   """
	add info to the targetfile header. The file header should be 
	indicated by headerlabel. 
   """	
   __funcname="addheaderinfo"  
   __errmsg="Error"+" in "+__funcname
   if (targetfile is ""):
      print __errmsg, "Enough numbers of arguments are not given." 
      sys.exit()    
   if os.path.isfile(targetfile) is False:
      print __errmsg, targetfile, "does not exist!!" 
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
               info, "    /' < ", __tmpfilename,
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


    



def fileseparater(targetfile="",textto="",resultdir="",niter=2,sectionlabel="%"):
   """
      split the sections in the targetfile into files and put them on 
      resultdir. The sections are defined by sectionlabel.
   """
   __funcname="addheaderinfo"  
   __errmsg="Error"+" in "+__funcname
   if (targetfile is "" or textto is "" or resultdir is ""):
      print __errmsg, "Enough numbers of arguments are not given." 
      sys.exit()    
   if os.path.isfile(targetfile) is False:
      print __errmsg, targetfile, "does not exist!!" 
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
	   os.system("unzip -qa "+zipfilename)   # Unzip the copied file
	   os.remove(zipfilename)	
	   fileseparater(targetfile=copiedfile,textto=textto,resultdir=resultdir)
## This loop add various information to files.
fileid=1
for parentdir, childdirs, files in os.walk(resultdir):
    for file in files:              # For each file under the targetdir,
	if file is not []:          # If there is a file, 
	   targetfile=os.path.join(parentdir,file)
	   deleteinfo(targetfile=targetfile,label="%")
           addinfoToEachRaw(targetfile=targetfile,info=str(fileid))
           fileid+=1           
	   addheaderinfo(targetfile=targetfile,info="%"+file)
	   addfooterinfo(targetfile=targetfile,info="%DONE")

