#!/usr/bin/env python


## From each file, take contents between textfrom and textto.
textfrom="%1"
textto="%coincidences"


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
parser.add_option("-v", "--verbose", dest="vrbflg",
                  help="Ask your intention for safety.",
		  default=1)
(options, args) = parser.parse_args()

vrbflg=options.vrbflg

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


    



def fileseparater(targetfile="",textto="",tmpdir="",niter=2,sectionlabel="%"):
   """
      split the sections in the targetfile into files and put them on 
      tmpdir. The sections are defined by sectionlabel.
   """
   __funcname="addheaderinfo"  
   __errmsg="Error"+" in "+__funcname
   if (targetfile is "" or textto is "" or tmpdir is ""):
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
       shutil.move(__resultfile,tmpdir)
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
   shutil.move(__resultfile,tmpdir)
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
   os.system("du -h "+resultfile)
   if vrbflg == 1:
       uinput=raw_input("Do you want to delete [d] the existing result file or append [a] results to it?:[d/a]")
       if uinput is "d":
           os.remove(resultfile)
       elif uinput is "a":
           pass
       else:
           sys.exit()
   
print "size of the target directory before unzip (if it is a zipped file)"
os.system("du -h "+targetdir+" | tail -n 1")
if vrbflg == 1:
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
           f=open(copiedfile,'r')
           buff=f.read(2);
           f.close()
           if buff == '%1':
               pass  # if the file starts with '%1', then we assume the file ascii, do nothing. 
           else:     # if not, we assume it a zip file and unzip it.
               zipfilename=copiedfile+"_tmp.zip"   
               shutil.move(copiedfile,zipfilename)
               os.system("unzip -qa "+zipfilename)   # Unzip the copied file
               os.remove(zipfilename)
	   fileseparater(targetfile=copiedfile,textto=textto,tmpdir=tmpdir)
## This loop add various information to files.
fileid=1
for parentdir, childdirs, files in os.walk(tmpdir):
    for file in files:              # For each file under the targetdir,
	if file is not []:          # If there is a file, 
	   targetfile=os.path.join(parentdir,file)
	   deleteinfo(targetfile=targetfile,label="%")
           addinfoToEachRaw(targetfile=targetfile,info=str(fileid))
           fileid+=1           
	   appendfile(fromfile=targetfile,tofile=resultfile)


addfooterinfo(targetfile=resultfile,info="%DONE")
os.system("rm -rf "+tmpdir)

