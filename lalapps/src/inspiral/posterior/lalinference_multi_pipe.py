#!/usr/bin/env @PYTHONPROG@

# DAG generation code for running LALInference pipeline
# (C) 2012 John Veitch
# 2013 Salvatore Vitale: extended to work with several ini files

from lalinference import lalinference_pipe_utils as pipe_utils
import ConfigParser
from optparse import OptionParser,OptionValueError
import sys

usage=""" %prog [options] config1.ini config2.ini ... configN.ini
Setup a Condor DAG file to run the LALInference pipeline based on
N config.ini files.

The user must specify either an injection file to analyse, with the --inj option,
a list of SnglInspiralTable or CoincInspiralTable triggers with the --<x>-triggers options,
or an ASCII list of GPS times with the --gps-time-file option.

The user must also specify and ini file which will contain the main analysis config.

"""

import os
def vararg_callback(option, opt_str, value, parser):
    assert value is None
    value = []
    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False
    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(arg)
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)

parser=OptionParser(usage)
parser.add_option("-r","--run-path",default=None,action="store",type="string",help="Directory to run pipeline in (default: $PWD)",metavar="RUNDIR")
parser.add_option("-p","--daglog-path",default=None,action="store",type="string",help="Path to directory to contain DAG log file. SHOULD BE LOCAL TO SUBMIT NODE",metavar="LOGDIR")
parser.add_option("-g","--gps-time-file",action="store",type="string",default=None,help="Text file containing list of GPS times to analyse",metavar="TIMES.txt")
parser.add_option("-t","--single-triggers",action="store",type="string",default=None,help="SnglInspiralTable trigger list",metavar="SNGL_FILE.xml")
parser.add_option("-C","--coinc-triggers",action="store",type="string",default=None,help="CoinInspiralTable trigger list",metavar="COINC_FILE.xml")
parser.add_option("-L","--lvalert",action="store",type="string",default=None,help="LVAlert coinc file",metavar="coinc_G0000.xml")
parser.add_option("--gid",action="store",type="string",default=None,help="Optional GraceDB ID for submitting results")
parser.add_option("-I","--injections",action="store",type="string",default=None,help="List of injections to perform and analyse",metavar="INJFILE.xml")
parser.add_option("-P","--pipedown-db",action="store",type="string",default=None,help="Pipedown database to read and analyse",metavar="pipedown.sqlite")
parser.add_option("-F","--folder-names",dest="fnames",action="callback", callback=vararg_callback,help="Space separated list of folders that will be created, corresponding to the TIGER parameters that are being tested or GR. The order has to be the same used with the ini files!",default=None,metavar="GR phi1")
parser.add_option("--condor-submit",action="store_true",default=False,help="Automatically submit the condor dag")
parser.add_option("-x", "--dax",action="store_true",default=False, help="Delete the ligo_data_find jobs and populate frame LFNs in the DAX -- WARNING: not tested in multi_pipe! Don't assume it will work.")
parser.add_option("-G", "--grid-site",action="store",type="string",metavar="SITE", help="Specify remote site in conjunction with --dax option. e.g. --grid-site=creamce for Bologna cluster.\
Supported options are: creamce and local",default=None)

(opts,args)=parser.parse_args()

if len(args)>1:
  print 'Using %s ini files\n'%len(args)
elif len(args)==1:
  inifile=args[0]

inits=args
ninits=len(inits)
fnames=opts.fnames
nfnames=len(fnames)

if not ninits==nfnames:
  print "You seem to be using %d parser files and %d foldernames. These two numbers must be the same. Exiting...\n"%(ninits,nfnames)
  sys.exit(1)

fnames_dic={}
for fname in fnames:
  fnames_dic[fnames.index(fname)]=str(fname)

glob_hyp=fnames
hyp_str=" "

for hy in glob_hyp:
  hyp_str+=hy+" "

cp=ConfigParser.ConfigParser()
cp.optionxform = str

first_dag=True
common_path=opts.run_path

for inifile in inits:
  cp.readfp(open(inifile))
  if opts.run_path is not None:
    cp.set('paths','basedir',opts.run_path)
  if opts.daglog_path is not None:
    cp.set('paths','daglogdir',opts.daglog_path)
  else:
    cp.set('paths','daglogdir',opts.run_path)
  if opts.gps_time_file is not None:
    cp.set('input','gps-time-file',opts.gps_time_file)
  if opts.single_triggers is not None:
    cp.set('input','sngl-inspiral-file',opts.single_triggers)
  if opts.injections is not None:
    cp.set('input','injection-file',opts.injections)
  if opts.coinc_triggers is not None:
    cp.set('input','coinc-inspiral-file',opts.coinc_triggers)
  if opts.lvalert is not None:
    cp.set('input','lvalert-file',opts.lvalert)
  if opts.gid is not None:
    cp.set('input','gid',opts.gid)
  if opts.pipedown_db is not None:
    cp.set('input','pipedown-db',opts.pipedown_db)

  if opts.run_path is not None:
    cp.set('paths','basedir',os.path.join(str(common_path),str(fnames_dic[inits.index(inifile)])))
    # Create the DAG from the configparser object
  if first_dag:
    dag=pipe_utils.LALInferencePipelineDAG(cp,site=opts.grid_site)
    first_dag=False
    dag.write_sub_files()
    dag2=dag
  else:
    dag2=pipe_utils.LALInferencePipelineDAG(cp,first_dag=False,previous_dag=dag,site=opts.grid_site)
    dag2.write_sub_files()    
    dag=dag2

# Create the DAG from the configparser object
dag2.set_dag_file(os.path.join(common_path,'common_dag'))
dag2.write_dag()
dag2.write_script()
# End of program
print 'Successfully created DAG file.'
print 'Now run condor_submit_dag %s\n'%(dag2.get_dag_file())

if opts.condor_submit:
    import subprocess
    from subprocess import Popen
           
    x = subprocess.Popen(['condor_submit_dag',dag.get_dag_file()])
    x.wait()
    if x.returncode==0:
      print 'Submitted DAG file'
    else:
      print 'Unable to submit DAG file'
