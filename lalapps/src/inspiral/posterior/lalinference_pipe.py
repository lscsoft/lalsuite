# DAG generation code for running LALInference pipeline
# (C) 2012 John Veitch, Vivien Raymond

from lalinference import lalinference_pipe_utils as pipe_utils
from lalapps import inspiralutils
import ConfigParser
from optparse import OptionParser,OptionValueError
import sys
import ast
import os
import uuid
from glue import pipeline

usage=""" %prog [options] config.ini
Setup a Condor DAG file to run the LALInference pipeline based on
the config.ini file.

The user can specify either an injection file to analyse, with the --inj option,
a list of SnglInspiralTable or CoincInspiralTable triggers with the --<x>-triggers options,
a GraceDB ID with the --gid option,
or an ASCII list of GPS times with the --gps-time-file option.

If none of the above options are given, the pipeline will analyse the entire
stretch of time between gps-start-time and gps-end-time, specified in the config.ini file.

The user must also specify and ini file which will contain the main analysis config.

"""
parser=OptionParser(usage)
parser.add_option("-r","--run-path",default=None,action="store",type="string",help="Directory to run pipeline in (default: $PWD)",metavar="RUNDIR")
parser.add_option("-p","--daglog-path",default=None,action="store",type="string",help="Path to directory to contain DAG log file. SHOULD BE LOCAL TO SUBMIT NODE",metavar="LOGDIR")
parser.add_option("-g","--gps-time-file",action="store",type="string",default=None,help="Text file containing list of GPS times to analyse",metavar="TIMES.txt")
parser.add_option("-t","--single-triggers",action="store",type="string",default=None,help="SnglInspiralTable trigger list",metavar="SNGL_FILE.xml")
parser.add_option("-C","--coinc-triggers",action="store",type="string",default=None,help="CoinInspiralTable trigger list",metavar="COINC_FILE.xml")
parser.add_option("--gid",action="store",type="string",default=None,help="GraceDB ID")
parser.add_option("-I","--injections",action="store",type="string",default=None,help="List of injections to perform and analyse",metavar="INJFILE.xml")
parser.add_option("-B","--burst_injections",action="store",type="string",default=None,help="SimBurst table for LIB injections",metavar="INJFILE.xml")
parser.add_option("-P","--pipedown-db",action="store",type="string",default=None,help="Pipedown database to read and analyse",metavar="pipedown.sqlite")
parser.add_option("--condor-submit",action="store_true",default=False,help="Automatically submit the condor dag")
parser.add_option("--pegasus-submit",action="store_true",default=False,help="Automatically submit the pegasus dax")
parser.add_option("-x", "--dax",action="store_true",default=False, help="Delete the ligo_data_find jobs and populate frame LFNs in the DAX")
parser.add_option("-G", "--grid-site",action="store",type="string",metavar="SITE", help="Specify remote site in conjunction with --dax option. e.g. --grid-site=creamce for Bologna cluster.\
Supported options are: creamce and local",default=None)

(opts,args)=parser.parse_args()

if len(args)!=1:
  parser.print_help()
  print 'Error: must specify one ini file'
  sys.exit(1)

inifile=args[0]

cp=ConfigParser.ConfigParser()
fp=open(inifile)
cp.optionxform = str
cp.readfp(fp)

approx='approx'
if cp.has_option('engine','approx'):
  pass
elif cp.has_option('engine','approximant'):
  approx='approximant'
else:
  print "Error: was expecting an 'approx' filed in the [engine] section\n"
  sys.exit(1)

apps=cp.get('engine',approx)
apps=apps.split(',')

samps=cp.get('analysis','engine')
samps=samps.split(',')

rundir_root=os.path.abspath(opts.run_path)
if opts.daglog_path is not None:
  daglogdir=os.path.abspath(opts.daglog_path)
elif opts.run_path is not None:
  daglogdir=os.path.abspath(opts.run_path)
else:
  daglogdir=os.path.abspath(cp.get('paths','basedir'))

use_roq=False
if cp.has_option('paths','roq_b_matrix_directory'):
  from numpy import genfromtxt
  path=cp.get('paths','roq_b_matrix_directory')
  if not os.path.isdir(path):
    print "The ROQ directory %s does not seem to exist\n"%path
    sys.exit(1)
  use_roq=True
  roq_paths=os.listdir(path)
  roq_params={}
  mc_ranges={}
  def key(item): # to order the ROQ bases
    return float(item[1]['seglen'])

  print "WARNING: Overwriting user choice of flow, srate, seglen,mc_min, mc_max and q-min"
  for roq in roq_paths:
    params=os.path.join(path,roq,'params.dat')
    roq_params[roq]=genfromtxt(params,names=True)
    mc_ranges[roq]=[float(roq_params[roq]['chirpmassmin']),float(roq_params[roq]['chirpmassmax'])]
  ordered_roq_paths=[item[0] for item in sorted(roq_params.items(), key=key)][::-1]
  i=0
  for roq in ordered_roq_paths:
    if i>0:
      # change min, just set to the max of the previous one since we have already aligned it in the previous iteration of this loop
      #mc_ranges[roq][0]+= (mc_ranges[roq_lengths[i-1]][1]-mc_ranges[roq][0])/2.
      mc_ranges[roq][0]=mc_ranges[ordered_roq_paths[i-1]][1]
    if i<len(roq_paths)-1:
      mc_ranges[roq][1]-= (mc_ranges[roq][1]- mc_ranges[ordered_roq_paths[i+1]][0])/2.
    i+=1
else:
  roq_paths=[None]

fp.close()


outerdaglog=os.path.join(daglogdir,'lalinference_multi_'+str(uuid.uuid1())+'.log')
outerdag=pipeline.CondorDAG(outerdaglog,dax=opts.dax)
outerdag.set_dag_file(os.path.join(rundir_root,'multidag'))


for sampler in samps:

  for app in apps:

    for roq in roq_paths:

      if not os.path.isdir(os.path.join(rundir_root,sampler,app)):
        os.makedirs(os.path.join(rundir_root,sampler,app))
      opts.run_path=os.path.abspath(os.path.join(rundir_root,sampler,app))

      inifile=args[0]

      cp=ConfigParser.ConfigParser()
      fp=open(inifile)
      cp.optionxform = str
      cp.readfp(fp)
      fp.close()
      cp.set('engine',approx,app)
      cp.set('analysis','engine',sampler)
      if roq is None:
        for p in dict(cp.items('paths')).keys():
          #current value
          if 'webdir' in p or 'url' in p or 'basedir' in p or 'daglogdir' in p:
            out=cp.get('paths',p)
            # append approximant prefix
            cp.set('paths',p,os.path.join(out,sampler,app))
      else:
        # do the appropriate hacks for ROQ
        if not os.path.isdir(os.path.join(rundir_root,sampler,app,roq)):
          os.makedirs(os.path.join(rundir_root,sampler,app,roq))
        opts.run_path=os.path.abspath(os.path.join(rundir_root,sampler,app,roq))
        for p in dict(cp.items('paths')).keys():
          #current value
          if 'webdir' in p or 'url' in p or 'basedir' in p or 'daglogdir' in p:
            out=cp.get('paths',p)
            # append approximant prefix
            cp.set('paths',p,os.path.join(out,sampler,app,roq))

        path=cp.get('paths','roq_b_matrix_directory')
        thispath=os.path.join(path,roq)
        cp.set('paths','roq_b_matrix_directory',thispath)
        mc_min=mc_ranges[roq][0]
        mc_max=mc_ranges[roq][1]
        flow=int(roq_params[roq]['flow'])
        srate=int(2.*roq_params[roq]['fhigh'])
        seglen=int(roq_params[roq]['seglen'])
        # params.dat uses the convention q>1 so our q_min is the inverse of their qmax
        q_min=1./float(roq_params[roq]['qmax'])
        cp.set('engine','srate',str(srate))
        cp.set('engine','seglen',str(seglen))
        tmp=cp.get('lalinference','flow')
        tmp=eval(tmp)
        for i in tmp.keys():
          tmp[i]=flow
        cp.set('lalinference','flow',str(tmp))
        cp.set('engine','chirpmass-min',str(mc_min))
        cp.set('engine','chirpmass-max',str(mc_max))
        cp.set('engine','q-min',str(q_min))
      if opts.condor_submit and opts.pegasus_submit:
          print 'Error: Please only specify one of --condor-submit or --pegasus-submit'
          sys.exit(1)

      if opts.run_path is not None:
        cp.set('paths','basedir',os.path.abspath(opts.run_path))

      if not cp.has_option('paths','basedir'):
        print 'Warning: No --run-path specified, using %s'%(os.getcwd())
        cp.set('paths','basedir',os.path.abspath(os.getcwd()))

      if opts.daglog_path is not None:
        cp.set('paths','daglogdir',os.path.abspath(opts.daglog_path))
      elif opts.run_path is not None:
        cp.set('paths','daglogdir',os.path.abspath(opts.run_path))
      else:
        cp.set('paths','daglogdir',os.path.abspath(cp.get('paths','basedir')))

      local_work_dir=cp.get('paths','daglogdir')

      if opts.gps_time_file is not None:
        cp.set('input','gps-time-file',os.path.abspath(opts.gps_time_file))

      if opts.single_triggers is not None:
        cp.set('input','sngl-inspiral-file',os.path.abspath(opts.single_triggers))

      if opts.injections is not None:
        cp.set('input','injection-file',os.path.abspath(opts.injections))

      if opts.burst_injections is not None:
        if opts.injections is not None:
          print "ERROR: cannot pass both inspiral and burst tables for injection\n"
          sys.exit(1)
        cp.set('input','burst-injection-file',os.path.abspath(opts.burst_injections))

      if opts.coinc_triggers is not None:
        cp.set('input','coinc-inspiral-file',os.path.abspath(opts.coinc_triggers))

      #if opts.lvalert is not None:
      #  cp.set('input','lvalert-file',os.path.abspath(opts.lvalert))

      if opts.gid is not None:
        cp.set('input','gid',opts.gid)

      if opts.pipedown_db is not None:
        cp.set('input','pipedown-db',os.path.abspath(opts.pipedown_db))


      # Create the DAG from the configparser object
      dag=pipe_utils.LALInferencePipelineDAG(cp,dax=opts.dax,site=opts.grid_site)
      if((opts.dax) and not cp.has_option('lalinference','fake-cache')):
        # Create a text file with the frames listed
        pfnfile = dag.create_frame_pfn_file()
        peg_frame_cache = inspiralutils.create_pegasus_cache_file(pfnfile)
      else:
        peg_frame_cache = '/dev/null'

      # A directory to store the DAX temporary files
      execdir=os.path.join(local_work_dir,'lalinference_pegasus_'+str(uuid.uuid1()))
      olddir=os.getcwd()
      os.chdir(cp.get('paths','basedir'))
      if opts.grid_site is not None:
        site='local,'+opts.grid_site
      else:
        site=None
      # Create the DAX scripts
      if opts.dax:
        dag.prepare_dax(tmp_exec_dir=execdir,grid_site=site,peg_frame_cache=peg_frame_cache)
        # Ugly hack to replace pegasus.transfer.links=true in the pegasus.properties files created by pipeline.py
        # Turns off the creation of links for files on the local file system. We use pegasus.transfer.links=false
        # to make sure we have a copy of the data in the runing directory (useful when the data comes from temporary
        # low latency buffer).
        if cp.has_option('analysis','pegasus.transfer.links'):
          if cp.get('analysis','pegasus.transfer.links')=='false':
            lines=[]
            with open('pegasus.properties') as fin:
              for line in fin:
                line = line.replace('pegasus.transfer.links=true', 'pegasus.transfer.links=false')
                lines.append(line)
            with open('pegasus.properties','w') as fout:
              for line in lines:
                fout.write(line)
        if cp.has_option('analysis','accounting_group'):
          lines=[]
          with open('sites.xml') as fin:
            for line in fin:
              if '<profile namespace="condor" key="getenv">True</profile>' in line:
                line=line+'    <profile namespace="condor" key="accounting_group">'+cp.get('analysis','accounting_group')+'</profile>\n'
              lines.append(line)
          with open('sites.xml','w') as fout:
            for line in lines:
              fout.write(line)
        if cp.has_option('analysis','accounting_group_user'):
          lines=[]
          with open('sites.xml') as fin:
            for line in fin:
              if '<profile namespace="condor" key="getenv">True</profile>' in line:
                line=line+'    <profile namespace="condor" key="accounting_group_user">'+cp.get('analysis','accounting_group_user')+'</profile>\n'
              lines.append(line)
          with open('sites.xml','w') as fout:
            for line in lines:
              fout.write(line)

      full_dag_path=os.path.join(cp.get('paths','basedir'),dag.get_dag_file())
      dagjob=pipeline.CondorDAGManJob(full_dag_path,dir=rundir_root)
      dagnode=pipeline.CondorDAGManNode(dagjob)
      outerdag.add_node(dagnode)

      dag.write_sub_files()
      dag.write_dag()
      dag.write_script()
      os.chdir(olddir)

if(opts.dax):
  # Create a text file with the frames listed
  pfnfile = outerdag.create_frame_pfn_file()
  peg_frame_cache = inspiralutils.create_pegasus_cache_file(pfnfile)

outerdag.write_sub_files()
outerdag.write_dag()
outerdag.write_script()

# End of program
print 'Successfully created DAG file.'
if not opts.dax:
  print 'Now run condor_submit_dag %s\n'%(outerdag.get_dag_file())

if opts.condor_submit:
  import subprocess
  from subprocess import Popen

  x = subprocess.Popen(['condor_submit_dag',outerdag.get_dag_file()])
  x.wait()
  if x.returncode==0:
    print 'Submitted DAG file'
  else:
    print 'Unable to submit DAG file'

if opts.pegasus_submit:
  import subprocess
  from subprocess import Popen

  os.chdir(os.path.abspath(cp.get('paths','basedir')))
  x = subprocess.Popen('./pegasus_submit_dax')
  x.wait()
  if x.returncode==0:
    print 'Submitted DAX file'
  else:
    print 'Unable to submit DAX file'
