#!/usr/local/bin/python2.3
"""
This script has been written to produce a file containing a table
describing the frequency bands and template banks used.
"""

# import required modules
import sys
import os
import getopt
import re
import string
import tempfile
import math
import ConfigParser

# append the lalapps python path
sys.path.append('@PYTHONLIBDIR@')

# program usage
def usage():
  msg = """\
Usage: GenerateFreqMeshFile [options]

  -h, --help               display this message
  -S, --sourcefile         the name of the source file
  -s, --source             the name of the source
  -B, --splitband          the size of the split up band
  -b, --meshband           the size of the mesh band
  -d, --pdet               the primary detector code
  -e, --sdet               the secondary detector code
  -D, --meshdir            the directory where the meshes are stored
  -n, --nodes              the number of nodes to split up band (optional)
  -r, --minres             the minimum frequency splitting resolution
  -o, --outfile            the name of the output file
  """
  print >> sys.stderr, msg

# parse the command line options to figure out what we should do
shortop = "hS:s:B:b:n:d:e:D:o:r:"
longop = [
"help",
"sourcefile=",
"source=",
"splitband=",
"meshband=",
"pdet=",
"sdet=",
"meshdir=",
"nodes=",
"minres=",
"outfile=",
]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

# default options
sourcefile = None
source = None
splitband = None
meshband = None
p_det = None
s_det = None
meshdir = None
nodes = 0
minres = 0.01
outfile = None

# process options
for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-S", "--sourcefile"):
    sourcefile = a
  elif o in ("-s", "--source"):
    source = a
  elif o in ("-B", "--splitband"):
    splitband = a
  elif o in ("-b", "--meshband"):
    meshband = a
  elif o in ("-d", "--pdet"):
    p_det = a
  elif o in ("-e", "--sdet"):
    s_det = a
  elif o in ("-D", "--meshdir"):
    meshdir = a
  elif o in ("-n", "--nodes"):
    nodes = a
  elif o in ("-r", "--minres"):
    minres = a
  elif o in ("-o", "--outfile"):
    outfile = a
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

if not sourcefile:
  print >> sys.stderr, "No sourcefile file specified."
  print >> sys.stderr, "Use --sourcefile FILE to specify location."
  sys.exit(1)

if not source:
  print >> sys.stderr, "No source specified."
  print >> sys.stderr, "Use --source PATH to specify a source."
  sys.exit(1)

if ((splitband==None)&(nodes==None)):
  print >> sys.stderr, "No search/injection bandwidth or nodes specified."
  print >> sys.stderr, "Use --splitband REAL8 or --nodes INT4 to specify a value."
  sys.exit(1)

if not meshband:
  print >> sys.stderr, "No mesh bandwidth specified."
  print >> sys.stderr, "Use --meshband REAL8 to specify a value."
  sys.exit(1)

if not outfile:
  print >> sys.stderr, "No output file specified."
  print >> sys.stderr, "Use --outfile FILE to specify a location."
  sys.exit(1)

if not p_det:
  print >> sys.stderr, "No primary detector specified."
  print >> sys.stderr, "Use --pdet STRING to specify a value."
  sys.exit(1)

if not s_det:
  print >> sys.stderr, "No secondary detector specified."
  print >> sys.stderr, "Use --sdet STRING to specify a value."
  sys.exit(1)

if not meshdir:
  print >> sys.stderr, "No mesh directory specified."
  print >> sys.stderr, "Use --meshdir PATH to specify a location."
  sys.exit(1)

######################################################################
# read source from sourcefile

try: sf = open(sourcefile, 'r')
except:
  print "ERROR : cannot open sourcefile"
  sys.exit(1)

for line in sf:
  temp = line.rstrip().split()
  if (temp[0]==source):
    f_A = (float)(temp[10])
    f_err_A = (float)(temp[11])
    f_B = (float)(temp[12])
    f_err_B = (float)(temp[13])
  
sf.close()

# calculate frequency band boundaries
f_min_A = f_A-f_err_A
f_max_A = f_A+f_err_A
f_min_B = f_B-f_err_B
f_max_B = f_B+f_err_B

######################################################################
# calculate the search/injections bands

f_min = []
f_max = []
band = []
i=0
j=0

if nodes!=0:
  nodes_A = (int)(math.floor((float)(nodes)*(f_err_A/(f_err_A+f_err_B))+0.5))
  nodes_B = (int)(math.floor((float)(nodes)*(f_err_B/(f_err_B+f_err_A))+0.5))

# for a double band
if f_B!=0.0:

  # calculate step sizes
  if nodes==0:
    step = float(splitband)
    n_step = math.floor((((float)(f_max_B)-(float)(f_min_B))/(float)(step))+0.5)
  else:
    step = float(minres)*math.floor(((float((float(f_max_B)-float(f_min_B))/float(nodes_B)))/float(minres))+0.5)
    nodes_B = (float(f_max_B)-float(f_min_B))/(float)(step)
    n_step = math.floor((float(nodes_B))+0.5)
   
  while i<n_step:
    f_max.append(f_max_B-(i*step))
    if f_max[i]>(f_min_B+step):
      band.append(step)
    else:
      band.append(f_max[i]-f_min_B)
    f_min.append(f_max[i]-band[i])
    i=i+1

# calculate step sizes
if nodes==0:
  step = float(splitband)
  n_step = (int)(math.floor(((float)(f_max_A)-(float)(f_min_A))/(float)(step)))
  #print "n_step is %d" %(n_step)
else:
  step = float(minres)*math.floor(((float((float(f_max_A)-float(f_min_A))/float(nodes_A)))/float(minres))+0.5)
  #print "step is %f" %(step)
  nodes_A = (float(f_max_A)-float(f_min_A))/(float)(step)
  #print "nodes_A is %f" %(nodes_A)
  n_step = math.floor((float(nodes_A))+0.5)
  #print "n_step is %d" %(n_step)
   
while j<n_step:
  f_max.append(f_max_A-(j*step))
  if f_max[i]>(f_min_A+step):
    band.append(step)
  else:
    band.append(f_max[i]-f_min_A)
  f_min.append(f_max[i]-band[i])
  

  i=i+1
  j=j+1

Nband = (int)(i)
#print "Nband is %d" %(Nband)

######################################################################
# calculate the mesh bands

fmesh_min = []
fmesh_max = []
bandmesh = []
i=0
j=0

# for a double band
if f_B!=0.0:

  # calculate step sizes
  step = float(meshband)
  n_step = math.floor((((float)(f_max_B)-(float)(f_min_B))/(float)(step))+0.5)
     
  while i<n_step:
    fmesh_max.append(f_max_B-(i*step))
    if fmesh_max[i]>(f_min_B+step):
      bandmesh.append(step)
    else:
      bandmesh.append(fmesh_max[i]-f_min_B)
    fmesh_min.append(fmesh_max[i]-bandmesh[i])

    i=i+1

# calculate step sizes
step = float(meshband)
n_step = math.ceil(((float)(f_max_A)-(float)(f_min_A))/(float)(step))
#print "n_step for bands is %d" %(n_step)
  
while j<n_step:
  fmesh_max.append(f_max_A-(j*step))
  if fmesh_max[i]>(f_min_A+step):
    bandmesh.append(step)
  else:
    bandmesh.append(fmesh_max[i]-f_min_A)
  fmesh_min.append(fmesh_max[i]-bandmesh[i])
  i=i+1
  j=j+1

Nmesh = i
#print "Nmesh is %d" %(Nmesh)

######################################################################
# Now we will organise the final table of results

pmesh = []
smesh = []
i=0
of = open(outfile, 'w')

# loop over the number of search/inj frequency bands
while i<Nband:

  # loop over the number of mesh bands
  j=0
  #print "Nmesh is %d" %(Nmesh)
  while j<Nmesh:

   
    # if the mesh band is correct
    if (fmesh_max[j]>=f_max[i])&(fmesh_min[j]<f_max[i]):   
      pmesh.append("%s/%s/mesh_%s_%s_%.6f.mesh" %(meshdir,p_det,p_det,source,fmesh_max[j]))
      smesh.append("%s/%s/mesh_%s_%s_%.6f.mesh" %(meshdir,s_det,s_det,source,fmesh_max[j]))
      #print "%f %f %f" %(fmesh_max[j],fmesh_min[j],f_max[i])
      #print "j is %d" %(j)
  
    j=j+1
    #print "j is %d" %(j)
    
  # output the results to file
  #of.write("%.3f %.3f %.3f\n" %(f_min[i], f_max[i], band[i]))
  of.write("%.3f %.3f %.3f %s %s\n" %(f_min[i], f_max[i], band[i], pmesh[i], smesh[i]))

  i=i+1
  #print "i is %d Nband is %d" %(i,Nband)

of.close()  

sys.exit(0)

