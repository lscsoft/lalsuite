import sys
import os
from optparse import *
import exceptions
import subprocess
import ConfigParser as cp

usage = """
        lvS5stat.dag generation script.
        """

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage, version="%prog CVS $Id$ " )
  parser.add_option("-d","--generate-dag",action="store_true",default=False,\
      help="generate dag" )
  parser.add_option("-f", "--config-file",action="store",type="string",\
    metavar=" FILE",help="use configuration file FILE")
  
  (options,args) = parser.parse_args()

  return options, sys.argv[1:]

opts, args = parse_command_line()  

# make a log directory if not made 
config=cp.ConfigParser()
config.read(opts.config_file)
# Execute the septime dag
septime=config.get('septime','lalapps_write_septime')
subprocess.call([str(septime)],shell=True)
subprocess.call(["condor_submit_dag septime_zero_lag.dag"],shell=True)
subprocess.call(["condor_submit_dag septime_slide.dag"],shell=True)
subprocess.call(["condor_submit_dag septime_injection.dag"],shell=True)

dagstatus = 1

while dagstatus==1:

  subprocess.call(["sleep 600"],shell=True)

  subprocess.call(["tail -1 septime_zero_lag.dag.dagman.out | awk '{print $10 \" \" $11}' > lastline_septime_zerolag.txt"],shell=True)
  subprocess.call(["tail -1 septime_slide.dag.dagman.out | awk '{print $10 \" \" $11}' > lastline_septime_slide.txt"],shell=True)
  subprocess.call(["tail -1 septime_injection.dag.dagman.out | awk '{print $10 \" \" $11}' > lastline_septime_injection.txt"],shell=True)

  zerolagfile = open("lastline_septime_zerolag.txt", "r")
  slidefile = open("lastline_septime_slide.txt", "r")
  injfile = open("lastline_septime_injection.txt", "r")
  for zline,sline,iline in zip(zerolagfile, slidefile, injfile):
    zexit=zline.rstrip()
    sexit=sline.rstrip()
    iexit=iline.rstrip()
    if zexit == "STATUS 0" and sexit == "STATUS 0" and iexit == "STATUS 0":
      print "The septime dags are complete!"
      dagstatus=0
    else:
      dagstatus=1
  zerolagfile.close()
  slidefile.close()
  injfile.close()


firstcoire=config.get('first-coire','lalapps_write_first_coire')
subprocess.call([str(firstcoire)],shell=True)
subprocess.call(["condor_submit_dag first_coire.dag"],shell=True)
subprocess.call(["condor_submit_dag first_coire_slide.dag"],shell=True)
subprocess.call(["condor_submit_dag first_coire_injection.dag"],shell=True)

firstcoiredagstatus = 1

while firstcoiredagstatus==1:

  subprocess.call(["sleep 600"],shell=True)

  subprocess.call(["tail -1 first_coire.dag.dagman.out | awk '{print $10 \" \" $11}' > lastline_first_coire.txt"],shell=True)
  subprocess.call(["tail -1 first_coire_slide.dag.dagman.out | awk '{print $10 \" \" $11}' > lastline_first_coire_slide.txt"],shell=True)
  subprocess.call(["tail -1 first_coire_injection.dag.dagman.out | awk '{print $10 \" \" $11}' > lastline_first_coire_injection.txt"],shell=True)

  firstcoirefile = open("lastline_first_coire.txt", "r")
  firstcoireslidefile = open("lastline_first_coire_slide.txt", "r")
  firstcoireinjfile = open("lastline_first_coire_injection.txt", "r")
  for zval,sval,ival in zip(firstcoirefile, firstcoireslidefile, firstcoireinjfile):
    zcode=zval.rstrip()
    scode=sval.rstrip()
    icode=ival.rstrip()
    if zcode == "STATUS 0" and scode == "STATUS 0" and icode == "STATUS 0":
      print "The first coire dags are complete!"
      firstcoiredagstatus=0
    else:
      firstcoiredagstatus=1
  firstcoirefile.close()
  firstcoireslidefile.close()
  firstcoireinjfile.close()

secondcoire=config.get('second-coire','lalapps_write_second_coire')
subprocess.call([str(secondcoire)],shell=True)
subprocess.call(["condor_submit_dag -no_submit second_coire.dag"],shell=True)
# Since the second*coire*slide*dag is run locally instead of on the nodes, set maxjobs to 1.
subprocess.call(["condor_submit_dag -no_submit -maxjobs 1 second_coire_slide.dag"],shell=True)
subprocess.call(["condor_submit_dag -no_submit second_coire_injection.dag"],shell=True)

corse=config.get('corse','lalapps_corse')
subprocess.call([str(corse)],shell=True)
# Since the corse*dag is run locally instead of on the nodes, set maxjobs to 1.
subprocess.call(["condor_submit_dag -no_submit -maxjobs 1 corse_all_data.dag"],shell=True)

use_eff=config.get('efficiency','use_efficiency')
subprocess.call([str(use_eff)],shell=True)
# Run locally, so set maxjobs to 1.
subprocess.call(["condor_submit_dag -no_submit -maxjobs 1 use_efficiency.dag"],shell=True)

dag_file=open("lvS5stat.dag","w")
dag_file.write("JOB 1 second_coire.dag.condor.sub" + "\n")
dag_file.write("JOB 2 second_coire_slide.dag.condor.sub" + "\n")
dag_file.write("JOB 3 second_coire_injection.dag.condor.sub" + "\n")
dag_file.write("JOB 4 corse_all_data.dag.condor.sub" + "\n")
dag_file.write("JOB 5 use_efficiency.dag.condor.sub" + "\n")
dag_file.write("PARENT 1 2 3 CHILD 4" + "\n")
dag_file.write("PARENT 4 CHILD 5" + "\n")
dag_file.close()

print "\nCreated a DAG file which can be submitted by executing"
print "\n   condor_submit_dag lvS5stat.dag"
print "\nfrom a condor submit machine (e.g. hydra.phys.uwm.edu)"

subprocess.call(["condor_submit_dag lvS5stat.dag"],shell=True)

dagstatus = 1

while dagstatus==1:

  subprocess.call(["sleep 600"],shell=True)

  subprocess.call(["tail -1 lvS5stat.dag.dagman.out | awk '{print $10 \" \" $11}' > lastline_lvS5stat.txt"],shell=True)

  dagfile = open("lastline_lvS5stat.txt", "r")
  for line in dagfile:
    dexit=line.rstrip()
    if dexit == "STATUS 0":
      print "The lvS5stat dag is complete!"
      dagstatus=0
    else:
      dagstatus=1
  dagfile.close()

subprocess.call(["lalapps_generate_upper_limits -f upper_limit.ini"],shell=True)

subprocess.call(["condor_submit_dag upper_limit.dag"],shell=True)

