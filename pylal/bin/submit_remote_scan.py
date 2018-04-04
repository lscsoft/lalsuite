#!/usr/bin/env python

__prog__ = "submit_remote_scan"
__title__ = "script to submit remote qscans via Cm messages"
__author__ = "Romain Gouaty"

import os
import tarfile
import shutil
import sys
import commands
import subprocess
import time
from optparse import *

from pylal import git_version

sys.path.append('@PYTHONLIBDIR@')

usage = """ %prog [options]
"""

parser = OptionParser(usage, version=git_version.verbose_msg)

parser.add_option("-g","--gps-time",action="store",type="string",\
    metavar=" GPS",help="GPS time to be used for the scan")

parser.add_option("-f","--config-file",action="store",type="string",\
    metavar=" PATH2FILE",help="Path to qscan configuration file")

parser.add_option("-t","--qscan-type",action="store",type="string",\
    metavar=" TYPE",help="qscan type")

parser.add_option("","--remote-output",action="store",type="string",\
    metavar=" PATH",help="output path for the remote scan")

#parser.add_option("-s","--sender",action="store",type="string",\
#    metavar=" STRING",help="name of the omega sender")

parser.add_option("","--remote-receiver",action="store",type="string",\
    metavar=" STRING",default="Omega",help="name of the remote omega receiver")

#parser.add_option("-r","--local-receiver",action="store",type="string",\
#    metavar=" STRING",help="name of the local omega receiver")

parser.add_option("-o","--output-path",action="store",type="string",\
    default="", metavar=" PATH",\
    help="path where the results will be stored")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#################################
# Main program

# Get a universal unique identifyer
uuid_temp = subprocess.Popen("uuidgen",stdout=subprocess.PIPE).communicate()[0]
uuid = uuid_temp.strip().replace("-","")

# Build a "UNAME"
system_name = subprocess.Popen("uname",stdout=subprocess.PIPE).communicate()[0]
hardware_platform = subprocess.Popen("uname -i",shell=True,stdout=subprocess.PIPE).communicate()[0]
uname = system_name.strip() + "-" + hardware_platform.strip()

# Get the "OMEGADROOT"
omegadroot = subprocess.Popen("echo $OMEGADROOT",shell=True,stdout=subprocess.PIPE).communicate()[0]


# Prepare the arguments to submit the scan
omegaSenderArg = "OmegaSender_" + uuid + " " + opts.remote_receiver + " CBC " + opts.qscan_type + " " + opts.gps_time + " " + opts.config_file + " " + opts.remote_output

# Prepare the executables 
omegadsend = omegadroot.strip() + "/" + uname + "/OmegaDSend.exe"
omegadreceive = omegadroot.strip() + "/" + uname + "/OmegaDReceive.exe"

# source to the correct environment
#source_env = subprocess.call("source /archive/home/romain/virgoApp/dot_bash.sh", shell=True)
#if source_env != 0:
#  print >> sys.stderr, "cm environment could not be sourced"
#  sys.exit(1)

#ip_info = subprocess.Popen("/sbin/ifconfig", shell=True, stdout=subprocess.PIPE).communicate()[0]
#print >> sys.stderr, ip_info

#cm_names = subprocess.Popen("/archive/home/romain/virgoApp/Cm/v8r4/Linux-x86_64/cm.exe names", shell=True, stdout=subprocess.PIPE).communicate()[0]
#print >> sys.stdout, "\"cm names\" returns: " + cm_names


print >> sys.stdout, "Running command: " + omegadreceive + " OmegaReceiver_" + uuid + " ./"

# Be ready to receive a message back from Cascina
print >> sys.stdout, "starting local receiver..."
receiver = subprocess.Popen(omegadreceive + " OmegaReceiver_" + uuid + " ./", shell=True)

# Wait for 10 seconds before launching the Omega sender
#time.sleep(10)

# Send a Cm message to Cascina
sender = subprocess.call(omegadsend + " " + omegaSenderArg, shell=True)
if sender != 0:
  print >> sys.stderr, "OmegaDSend returned an error when run with these arguments: " + omegaSenderArg
  sys.exit(1)
else:
  print >> sys.stdout, "The following Cm message has been sent to Cascina: " + "OmegaDSend " + omegaSenderArg

# Wait for the receiver to finish and check for the result
receiver_result = os.waitpid(receiver.pid, 0)
#print receiver_result
if receiver_result[1] != 0:
  print >> sys.stderr, "OmegaDReceive command failed: OmegaReceiver_" + uuid
  sys.exit(1)
else:
  print >> sys.stdout, "OmegaDReceive finished succesfully: OmegaReceiver_" + uuid

# Untar the result 
result_output_file = "omegaOut_" + uuid + ".tar.gz"
if os.path.exists(result_output_file):
  if tarfile.is_tarfile(result_output_file):
    tar = tarfile.open(result_output_file,"r:gz")
    file_list = tar.getnames()
    for file in file_list:
      tar.extract(file)
    tar.close()
    os.remove(result_output_file)
  else:
    print >> sys.stderr, "File " + result_output_file + " is not a valid tar file!!"
    sys.exit(1)
else:
  print >> sys.stderr, "File " + result_output_file + " could not be found!!"
  sys.exit(1)

# Move the result to the outputpath
result_directory = "scanResults_" + uuid
if not os.path.exists(result_directory):
  print >> sys.stderr, "Directory " + result_directory + " could not be found!!"
  sys.exit(1)
else:
  if os.path.exists(opts.output_path):
    print >> sys.stdout, "Cleaning directory " + opts.output_path + "..."
    shutil.rmtree(opts.output_path)
  print >> sys.stdout, "moving results to " + opts.output_path + "..."
  shutil.move(result_directory,opts.output_path)  

