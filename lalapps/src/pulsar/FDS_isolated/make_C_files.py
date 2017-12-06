#!/usr/bin/env python2
"""
M.Alessandra Papa, July 05
"""

# import standard modules to the python path
import sys, os, shutil, math,random
import getopt, re, string,popen2
import ConfigParser
sys.path.append('/usr/lib/python2')

pi=3.14159265358979323844

# Function usage
def usage():
  msg = """\
Usage: allsky_pulsar_pipe.py [options]

  -h, --help               display this message
  -S  --input-dir          input dir where tot Confidence files are
  -f, --start-freq         first frequency 
  -N  --how-many           how many files
  -W  --ouput-dir     Output dir
"""
  print >> sys.stderr, msg


# ------------------- parse the command line options -------------------- #
# initialise command line argument variables
#------------------------------------------------------------------------ #

input_dir = None
output_dir = None
how_many = -1
start_freq = -1.0

shortop = "hN:S:W:f:"
longop = [
   "help",
   "how-many=",
   "input-dir=",
   "output-dir=",
   "start-freq=",
   ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

for o, a in opts:
  if o in ("-f", "--start-freq"):
    start_freq = float(a)
  elif o in ("-S", "--input-dir"):
    input_dir = a      
  elif o in ("-W", "--output-dir"):
    output_dir = a      
  elif o in ("-N", "--how-many"):
    how_many = int(a)      
  elif o in ("-h", "--help"):
    usage()
    sys.exit(0)
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

print " ==========="
print "Give it an input-dir where the code linearfit.run is"
print "and where the subdir Ctot lives."
print "In output_dir you will find a file called h095"
print "with 3 columns: freq, dC/dh0, h095"
print "and a directory, NewCtot, that contains new Confidence.data files"
print "with the tolerance computed correctly from the binomial sigma."
print " ==========="

if start_freq == -1.0:
  print >> sys.stderr, "No job starting freq specified."
  print >> sys.stderr, "Use --start-freq to specify it."
  sys.exit(1)
if not input_dir:
  print >> sys.stderr, "No input directory specified."
  print >> sys.stderr, "Use --input-dir to specify it (w/out slash)."
  sys.exit(1)
if not output_dir:
  print >> sys.stderr, "No output directory specified."
  print >> sys.stderr, "Use --output-dir to specify it (w/out slash)."
  sys.exit(1)
if how_many == -1:
  print >> sys.stderr, "Did not specify how many files"
  print >> sys.stderr, "Use --how-many to specify it."
  sys.exit(1)


# -------------------------------------------------------------------------------- #

# ------------------------- Set up working directory ----------------------------- # 


ifile = 0
coincidence_band=.06
Npolka=20

# make output directory
if not os.path.exists(output_dir):
  print output_dir
  try: os.mkdir(output_dir)
  except OSError, err:
    import errno
    print "Warning:", err
  os.mkdir(output_dir)
  sys.stdout.flush()


subdir=''.join([output_dir,'/NewCtot'])
print subdir
if not os.path.exists(subdir):
  print subdir,' does not exist. Creating...'
  # make output sub-directory
  try: os.mkdir(subdir)
  except OSError, err:
    import errno
    print "Warning:", err  
  os.mkdir(subdir)
  sys.stdout.flush()

h095_out_filename=''.join([output_dir,'/h095'])
h095_file=open(h095_out_filename,mode='w')

print how_many
while ifile < how_many:

  #define this ere because it enters in confidence.data file name
  freq=float(start_freq) + ifile * float(coincidence_band) * Npolka
  ifile=ifile+1

  res_in=''.join([input_dir,'/Ctot/Confidence.data-',str(freq)])
  res_out=''.join([subdir,'/Confidence.data-',str(freq)])
  #print '====='
  #print res_in
  #print res_out

  if os.path.exists(res_in):

    print res_in, ' file exists. Reading ...'

    Cdata_file=open(res_in,mode='r')
    NewCdata_file=open(res_out,mode='w')
    
    line_list=Cdata_file.readlines()
    length=len(line_list)
    #print 'length=',length
    #print line_list

    iline=0
    while iline < length:

      line=line_list[iline]
      #print line
      [sNinj,stol,sh0,sdh0,sconfidence]=line.split(None,5)
      confidence=float(sconfidence)
      Ninj=int(sNinj)
      tol=math.sqrt(confidence*(1.0-confidence)/Ninj)

      if Ninj > 3000 :
        if confidence > 0.925:
          print >>NewCdata_file,Ninj,tol,sh0,sdh0,confidence
      iline = iline+1
      #endif Ninj > 1000
    #end iline
      
    Cdata_file.close()
    NewCdata_file.close()

    cmd=''.join(['./linearfit.run ',res_out,' > dmp'])
    os.system(cmd)
    #print cmd

    dmp_file=open('dmp',mode='r')
    line_list=dmp_file.readlines()
    length=len(line_list)
    line=line_list[length-1]
    #print line
    [sh095,sSlope,sCerr,schi2,sDhh,sDhl,sDhp,smu,sstd,sprob]=line.split(None,10)
    dmp_file.close()

#sh095 95% confidence
#sSlope slope of the best fit liear fit
#sCerr 1 sigma error on the confidence
#schi2 chi square value
#sDhh,sDhl lower and upper C=95%  intercepts for the +/- 1 sigma C fits
#sDhp percentage total error on h095 (Dhh+Dhl)/h095

    print >>h095_file,freq,sh095,sSlope,sCerr,schi2,sDhh,sDhl,sDhp,smu,sstd,sprob
  else:
    print 'File ',res_in,'non esiste!'
  #endif this Confidence_tot_freq file exists
#end the loop over different freq files

h095_file.close()


# ------------------------------------------------------------------------------------ #
