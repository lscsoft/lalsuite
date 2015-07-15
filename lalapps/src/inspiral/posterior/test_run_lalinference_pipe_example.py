#!/usr/bin/env python

import argparse
import shutil
import distutils.spawn
import os
import lalinference

parser = argparse.ArgumentParser(description="Runs a test run of lalinference using lalinference_pipe_example.ini.")

parser.add_argument('-i','--ini_file', type=str, nargs='?',
                    default='lalinference_pipe_example.ini',
                    help='lalinference_pipe ini file to process.')

parser.add_argument('-x','--injection_file', type=str, nargs='?',
                    default='fiducialBNS.xml',
                    help='injection file.')

parser.add_argument('-e','--engine', type=str, nargs='?',
                    default='lalinferencemcmc',
                    help='lalinference engine to run with.')

parser.add_argument('-o','--output', type=str, nargs='?',
                    default=None,
                    help='output directory.')

args = parser.parse_args()


if 'UNCLEAN' in lalinference.InferenceVCSId:
    default_outputdir=os.getenv('HOME')+'/public_html/lalinference_testrun/'+lalinference.InferenceVCSId+'_UNCLEAN/'+args.engine+'/'
else:
    default_outputdir=os.getenv('HOME')+'/public_html/lalinference_testrun/'+lalinference.InferenceVCSId+'/'+args.engine+'/'


if args.output == None:
    args.output=default_outputdir

backup_file=args.output+'/'+os.path.basename(args.ini_file)+'.bak'
ini_file=args.output+'/'+os.path.basename(args.ini_file)

os.makedirs(args.output)

shutil.copy(args.injection_file,args.output)
shutil.copy(args.ini_file,backup_file)

try:
    shutil.copy(args.ini_file,ini_file)
except:
    print ini_file+' will be modified.'

os.chdir(args.output)


path_keys = {'datafind': 'ligo_data_find',
            'mergescript': 'lalapps_nest2pos',
            'resultspage': 'cbcBayesPostProc.py',
            'segfind': 'ligolw_segment_query',
            'ligolw_print': 'ligolw_print',
            'coherencetest': 'lalapps_coherence_test',
            'lalinferencenest': 'lalinference_nest',
            'lalinferencemcmc': 'lalinference_mcmc',
            'lalinferencebambi': 'lalinference_bambi',
            'skyarea': 'run_sky_area.py',
            'mpirun': 'mpirun',
            'mpiwrapper': 'lalinference_mpi_wrapper',
            'gracedb': 'gracedb',
            'ppanalysis': 'cbcBayesPPAnalysis.py',
            'pos_to_sim_inspiral': 'cbcBayesPosToSimInspiral.py'}


def replace(line):
    for key in path_keys.keys():
        if key+'=/' in line:
            albert_path=line.split('=')[-1]
            exec_path=distutils.spawn.find_executable(path_keys[key])
            if exec_path==None:
                exec_path=distutils.spawn.find_executable('true')
                print 'Could not find executable for '+path_keys[key]+'. Will use '+exec_path+' instead.'
            new_path=exec_path+'\n'
            return line.replace(albert_path,new_path)
    if 'engine=' in line:
        return line.replace(line.split('=')[-1],args.engine)
    if 'nparallel=' in line:
        return '# '+line
    if 'webdir=' in line:
        return line.replace(line.split('=')[-1],os.getcwd()+'/webdir/')
    if 'baseurl=' in line:
        return line.replace(line.split('=')[-1],'file://'+os.getcwd()+'./webdir/')
    if 'fake-cache=' in line:
        return line.replace(line,"fake-cache={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}")
    if 'ignore-science-segments=' in line:
        return 'ignore-science-segments=True'
    if 'dataseed=' in line:
        return line.replace('#','').strip()+'\n'
    if 'margphi=' in line:
        return line.replace('#','').strip()+'\n'
    if 'parname-max' in line:
        return line+'distance-max=500\n'
    if 'deltaLogL=' in line:
        return line.replace('#','').strip()+'\n'
    if 'accounting_group=' in line:
        return line.replace('#','').strip()+'\n'
    return line


with open(ini_file,'w') as fout:
    with open(backup_file,'r') as fin:
        for line in fin:
            fout.write(replace(line))

lalinferenceargs = [ 'lalinference_pipe'
		     , '-I'
		     , args.injection_file
		     , '-r'
		     , './run'
		     , '-p'
		     , './daglog'
		     , ini_file
		     , '--dax'
		     , '--grid-site'
		     , 'local'
		     , '--pegasus-submit'
                     #, '--condor-submit'
		     ]

os.execlp('lalinference_pipe', *lalinferenceargs)
