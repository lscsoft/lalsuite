#!/usr/bin/env python2.2
"""
pulsar_pipe.in - standalone pulsar pipeline driver script
X. Siemens
"""

# import standard modules to the python path
import sys, math
import getopt
sys.path.append('/usr/lib/python2.2')

#./ConfidenceErrors.py -s 160.0 -e 760.0 -b 1.2

# ------------------- parse the command line options -------------------- #

file_basename = None
starting_freq = 0
end_freq = 0
band = None

shortop = "s:e:b:"
longop = [
   "start-freq=",
   "end-freq=",
   "band=",
   ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  print 'Exiting'
  sys.exit(1)
  
for o, a in opts:
  if o in ("-s", "--start-freq"):
    starting_freq = float(a)
  elif o in ("-e", "--end-freq"):
    end_freq = float(a)      
  elif o in ("-b", "--band"):
    band = float(a)      
  else:
    print >> sys.stderr, "Unknown option:", o
    sys.exit(1)

# -------------------------------------------------------------------------------- #

res_file=open('ConfidenceErrors-'+str(starting_freq)+'-'+str(end_freq),'w')

print 'Writing to ConfidenceErrors-'+str(starting_freq)+'-'+str(end_freq)
# First loop over files:

f=starting_freq
while f < end_freq:
  
  file='Confidence.data-'+str(f)
  conf_file=open(file,mode='r')

  print 'Opening Confidence.data-'+str(f)

  # loop ove lines in each file:find out how many lines there are in the file
  no_lines=0
  for line in conf_file:
    no_lines = no_lines +1
  conf_file.close()  

  total = no_lines

  print '    total lines:',str(total)
  
  conf_file=open(file,mode='r')
  no_lines=0  
  for line in conf_file:
    #get second to last line
    if no_lines == (total-2):
      [Ninj1,tol1,h01,dh01,conf1]=line.split(None,5)
    #get last line
    if no_lines == (total-1):
      [Ninj2,tol2,h02,dh02,conf2]=line.split(None,5)
    no_lines = no_lines +1
  conf_file.close()  


  #OK now have last and second to last lines of file
  Ninj=float(Ninj2)
  if float(Ninj1) < float(Ninj2):
    Ninj=float(Ninj1)
    
  DeltaC = 1.0/math.sqrt(Ninj)
  dC = float(conf2)+1.0/math.sqrt(float(Ninj2))-(float(conf1)-1.0/math.sqrt(float(Ninj1)))
  dh = float(h02)-float(h01)
  
  if dC == 0 or dh == 0 or total == 0 or total == 1:
    print 'Cannot determine error from last two lines of file: ', file
    dC = 1
    dh = 0
  
  hprime=dh/dC

  deltah=hprime*DeltaC
  
#  print >>res_file, f, Ninj1, Ninj2, conf1, conf2, h01, h02, deltah, deltah/float(h02)
  print >>res_file, f, h02, deltah,deltah/float(h02) 

  f=f+band
     
res_file.close()

# -------------------------------------------------------------------------------- #

