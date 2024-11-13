##python
# End-to-end LALInference test pipeline
# (C) 2014 John Veitch

from __future__ import print_function

import os
import sys
import uuid
from optparse import OptionParser

from configparser import ConfigParser

from glue import pipeline

from lalinference import lalinference_pipe_utils as pipe_utils
from lalinference.lalinference_pipe_utils import mkdirs

usage=""" %prog [options] config.ini
Setup a DAG to run an end-to-end lalinference test:
 1) Generate samples from prior
 2) Analyse a set of injections drawn from the prior
 3) Run P vs P test on results
"""

parser=OptionParser(usage)
parser.add_option("-r","--run-path",default='./',action="store",type="string",help="Directory to run pipeline in (default: $PWD)",metavar="RUNDIR")
parser.add_option("-p","--daglog-path",default=None,action="store",type="string",help="Path to directory to contain DAG log file. SHOULD BE LOCAL TO SUBMIT NODE",metavar="LOGDIR")
parser.add_option("-I","--injections",default=None,action="store",type="string",help="Path to injection file, will bypass the prior-sampling stage",metavar="injections.xml")
parser.add_option('-N','--trials',action='store',type='int',metavar='NUM',help='Number of prior samples to analyse')


(opts,args)=parser.parse_args()

if len(args)==0:
        parser.print_help()
        sys.exit(1)

inifile=args[0]

# Set up the configuration for the sub-dags

prior_cp=ConfigParser()
prior_cp.optionxform = str
try:
  prior_cp.read_file(open(inifile))
except AttributeError:
  prior_cp.readfp(open(inifile))

main_cp=ConfigParser()
main_cp.optionxform = str
try:
  main_cp.read_file(open(inifile))
except AttributeError:
  main_cp.readfp(open(inifile))


rundir=os.path.abspath(opts.run_path)

if opts.daglog_path is not None:
  prior_cp.set('paths','daglogdir',os.path.join(os.path.abspath(opts.daglog_path),'prior'))
  main_cp.set('paths','daglogdir',os.path.join(os.path.abspath(opts.daglog_path),'main'))
  daglogdir=os.path.abspath(opts.daglog_path)
else:
  prior_cp.set('paths','daglogdir',os.path.join(os.path.abspath(opts.run_path),'prior'))
  main_cp.set('paths','daglogdir',os.path.join(os.path.abspath(opts.run_path),'main'))
  daglogdir=os.path.abspath(opts.run_path)

webdir=main_cp.get('ppanalysis','webdir')
priordir=os.path.join(rundir,'prior')
maindir=os.path.join(rundir,'main')
priorwebdir=os.path.join(webdir,'prior')
mainwebdir=os.path.join(webdir,'injections')
prior_cp.set('paths','basedir',priordir)
main_cp.set('paths','basedir',maindir)
prior_cp.set('paths','webdir',priorwebdir)
main_cp.set('paths','webdir',mainwebdir)

outerlogdir=os.path.join(daglogdir,'log')
mkdirs(outerlogdir)
mkdirs(priordir)
mkdirs(maindir)
mkdirs(priorwebdir)
mkdirs(mainwebdir)

# Add the prior options to the sub dag
if prior_cp.get('analysis','engine')=='lalinferencenest':
  prior_cp.set('engine','sampleprior',str(20*opts.trials)) # more samples drawn since not all will end up in posterior
  prior_cp.set('engine','zeroLogLike','')
  prior_cp.set('engine','nlive',str(20*opts.trials))
if prior_cp.get('analysis','engine')=='lalinferenceburst':
  prior_cp.set('engine','sampleprior',str(20*opts.trials)) # more samples drawn since not all will end up in posterior
  prior_cp.set('engine','zeroLogLike','')
  prior_cp.set('engine','nlive',str(20*opts.trials))
elif prior_cp.get('analysis','engine')=='lalinferencemcmc':
  prior_cp.set('engine','Neff',str(opts.trials))
  prior_cp.set('engine','zeroLogLike','')
elif prior_cp.get('analysis','engine')=='lalinferencebambi':
  prior_cp.set('engine','zeroLogLike','')
  prior_cp.set('engine','nlive',str(opts.trials))
elif prior_cp.get('analysis','engine')=='lalinferencebambimpi':
  prior_cp.set('engine','zeroLogLike','')
  prior_cp.set('engine','nlive',str(opts.trials))

# Create a DAG to contain the other scripts
outerdaglog=os.path.join(daglogdir,'lalinference_injection_test_'+str(uuid.uuid1())+'.log')
outerdag=pipeline.CondorDAG(outerdaglog)
outerdag.set_dag_file(os.path.join(rundir,'priortest'))

# Run code with prior sampling
trig_time=1085855789
fake_event=pipe_utils.Event(trig_time=trig_time)
tfpath=os.path.join(rundir,'time.txt')
tfile=open(tfpath,'w')
print('%i\n'%(trig_time), file=tfile)
tfile.close()
prior_cp.set('input','gps-time-file',tfpath)

priordag=pipe_utils.LALInferencePipelineDAG(prior_cp)
priordag.set_dag_file(os.path.join(priordir,'lalinference_priorsample'))
priordagjob=pipeline.CondorDAGManJob(priordag.get_dag_file(),dir=priordir)
priordagnode=pipeline.CondorDAGManNode(priordagjob)
# Find the output file
pagenode=filter(lambda n:isinstance(n,pipe_utils.ResultsPageNode), priordag.get_nodes())[0]
priorfile=pagenode.get_pos_file()

# Convert prior samples to injections
convertsub=os.path.join(rundir,'samples2injections.sub')
converterr=os.path.join(outerlogdir,'samples2injection-$(cluster)-$(process)-$(node).err')
convertout=os.path.join(outerlogdir,'samples2injection-$(cluster)-$(process)-$(node).out')

if opts.injections:
        injfile=os.path.abspath(opts.injections)
else:
        injfile=os.path.join(rundir,'priorsamples.xml')
approx=prior_cp.get('engine','approx')
prior2injexe=prior_cp.get('condor','pos_to_sim_burst')
prior2injjob=pipeline.CondorDAGJob('vanilla',prior2injexe)
if main_cp.has_option('analysis','accounting_group'):
  prior2injjob.add_condor_cmd('accounting_group',main_cp.get('analysis','accounting_group'))
prior2injjob.set_sub_file(convertsub)
prior2injjob.set_stderr_file(converterr)
prior2injjob.set_stdout_file(convertout)
prior2injjob.add_condor_cmd('getenv','True')
prior2injnode=pipeline.CondorDAGNode(prior2injjob)
prior2injnode.add_var_opt('output',injfile)
prior2injnode.add_var_opt('num-of-injs',str(opts.trials))
prior2injnode.add_var_opt('approx',approx)
prior2injnode.add_var_arg(priorfile)
prior2injnode.add_parent(priordagnode)

# Create the pipeline based on the injections
#main_cp.set('input','injection-file',injfile)
main_cp.set('input','gps-start-time',str(trig_time-1000))
main_cp.set('input','gps-end-time',str(trig_time+1000))
maindag=pipe_utils.LALInferencePipelineDAG(main_cp)
maindag.set_dag_file(os.path.join(maindir,'lalinference_pipeline'))
maindagjob=pipeline.CondorDAGManJob(maindag.get_dag_file(),dir=maindir)
maindagnode=pipeline.CondorDAGManNode(maindagjob)
maindag.config.set('input','burst-injection-file',injfile)
for i in range(int(opts.trials)):
  ev=pipe_utils.Event(trig_time=trig_time,event_id=i)
  e=maindag.add_full_analysis(ev)

outerdag.add_node(maindagnode)

if not opts.injections:
    outerdag.add_node(priordagnode)
    outerdag.add_node(prior2injnode)
    maindagnode.add_parent(prior2injnode)

# Get a list of posterior samples files
resultspagenodes=filter(lambda n: isinstance(n, pipe_utils.ResultsPageNode), maindag.get_nodes())
posteriorfiles=[n.get_pos_file() for n in resultspagenodes]

# Analyse results of injection runs to generate PP plot
ppsub=os.path.join(rundir,'ppanalysis.sub')
pperr=os.path.join(outerlogdir,'ppanalysis-$(cluster)-$(process)-$(node).err')
ppout=os.path.join(outerlogdir,'ppanalysis-$(cluster)-$(process)-$(node).out')
ppexe=prior_cp.get('condor','ppanalysis')
ppjob=pipeline.CondorDAGJob('vanilla',ppexe)
ppjob.set_sub_file(ppsub)
ppjob.set_stderr_file(pperr)
ppjob.set_stdout_file(ppout)
ppjob.add_condor_cmd('getenv','True')
if main_cp.has_option('analysis','accounting_group'):
  ppjob.add_condor_cmd('accounting_group',main_cp.get('analysis','accounting_group'))

ppnode=pipeline.CondorDAGNode(ppjob)
ppnode.add_var_opt('injXML',injfile)
if main_cp.has_option('ppanalysis','webdir'):
  outdir=main_cp.get('ppanalysis','webdir')
else:
  outdir=os.path.join(rundir,'ppanalysis')

mkdirs(outdir)
ppnode.add_var_opt('outdir',outdir)
for f in posteriorfiles:
  ppnode.add_var_arg(f)

ppnode.add_parent(maindagnode)
outerdag.add_node(ppnode)


outerdag.write_sub_files()
outerdag.write_dag()
outerdag.write_script()

if not opts.injections:
    priordag.write_sub_files()
    priordag.write_dag()
    priordag.write_script()

maindag.write_sub_files()
maindag.write_dag()
maindag.write_script()

# End of program
print('Successfully created DAG file.')
print('Now run condor_submit_dag %s\n'%(outerdag.get_dag_file()))


