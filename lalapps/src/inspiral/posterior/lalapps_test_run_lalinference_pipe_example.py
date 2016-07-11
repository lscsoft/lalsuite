
import argparse
import shutil
import distutils.spawn
import os
import subprocess
import glob
import lalinference

parser = argparse.ArgumentParser(description="Runs test runs of lalinference using lalinference_pipe_example.ini.")

parser.add_argument('-i','--ini_file', type=str, nargs='?',
                    default='lalinference_pipe_example.ini',
                    help='lalinference_pipe ini file to process.')

parser.add_argument('-x','--injection_file', type=str, nargs='?',
                    default='fiducialBNS.xml',
                    help='injection file.')

parser.add_argument('--gracedb', action='store_true',
                    default=False,
                    help='Runs the analysis for the GraceDB test event T169545.')

parser.add_argument('--bbh-injection', type=str, nargs='?',
                    default='',
                    help='injection file for optional BBH analysis.')

parser.add_argument('-e','--engine', type=str, nargs='?',
                    default='lalinferencemcmc',
                    help='lalinference engine to run with.')

parser.add_argument('-o','--output', type=str, nargs='?',
                    default=None,
                    help='output directory.')

parser.add_argument('--condor-submit', action='store_true',
                    default=False,
                    help='submits the test suite.')

args = parser.parse_args()


if 'UNCLEAN' in lalinference.InferenceVCSId:
    default_outputdir=os.getenv('HOME')+'/lalinference_testrun/'+lalinference.InferenceVCSId+'_UNCLEAN/'+args.engine+'/'
else:
    default_outputdir=os.getenv('HOME')+'/lalinference_testrun/'+lalinference.InferenceVCSId+'/'+args.engine+'/'


if args.output == None:
    args.output=default_outputdir
else:
    args.output=os.path.abspath(args.output)

args.injection_file=os.path.abspath(args.injection_file)
if args.bbh_injection:
    args.bbh_injection=os.path.abspath(args.bbh_injection)

backup_file=args.output+'/'+os.path.basename(args.ini_file)+'.bak'
ini_file=args.output+'/'+os.path.basename(args.ini_file)

os.makedirs(args.output)

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
            'lalinferencedatadump': 'lalinference_datadump',
            'bayesline': 'BayesLine',
            'skyarea': 'run_sky_area',
            'mpirun': 'mpirun',
            'mpiwrapper': 'lalinference_mpi_wrapper',
            'gracedb': 'gracedb',
            'ppanalysis': 'cbcBayesPPAnalysis.py',
            'pos_to_sim_inspiral': 'cbcBayesPosToSimInspiral.py'}

def replace(line):
    if line[0]=='#':
        return line
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
    if 'accounting_group=' in line:
        return line.replace('#','').strip()+'\n'
    return line

with open(ini_file,'w') as fout:
    with open(backup_file,'r') as fin:
        for line in fin:
            fout.write(replace(line))


############################################################

os.makedirs(args.output+'/fiducialBNS/')
os.chdir(args.output+'/fiducialBNS/')

shutil.copy(args.injection_file,args.output+'/fiducialBNS/')
shutil.copy(ini_file,args.output+'/fiducialBNS/'+os.path.basename(ini_file)+'.bak')
shutil.copy(ini_file,args.output+'/fiducialBNS/')

def replace_fiducial(line):
    if 'webdir=' in line:
        return line.replace(line.split('=')[-1],os.getcwd()+'/webdir/')
    if 'baseurl=' in line:
        return line.replace(line.split('=')[-1],'file://'+os.getcwd()+'/webdir/')
    if 'fake-cache=' in line:
        return line.replace(line,"fake-cache={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}")
    if 'ignore-science-segments=' in line:
        return 'ignore-science-segments=True\n'
    if 'dataseed=' in line:
        return line.replace('#','').strip()+'\n'
    if 'margphi=' in line:
        return line.replace('#','').strip()+'\n'
    if 'amporder=' in line:
        return 'amporder=0\n'
    if 'parname-max' in line:
        return line+'distance-max=500\n'
    if 'deltaLogL=' in line:
        return line.replace('#','').strip()+'\n'
    return line

with open(args.output+'/fiducialBNS/'+os.path.basename(ini_file),'w') as fout:
    with open(args.output+'/fiducialBNS/'+os.path.basename(ini_file)+'.bak','r') as fin:
        for line in fin:
            fout.write(replace_fiducial(line))

lalinferenceargs = [ 'lalinference_pipe'
		     , '-I'
		     , args.injection_file
		     , '-r'
		     , './run'
		     , '-p'
		     , './daglog'
		     , args.output+'/fiducialBNS/'+os.path.basename(ini_file)
                     ]

if args.condor_submit:
    lalinferenceargs.append('--condor-submit')

subprocess.call(lalinferenceargs)

###################################################

def replace_GraceDB(line):
    if 'webdir=' in line:
        return line.replace(line.split('=')[-1],os.getcwd()+'/webdir/')
    if 'baseurl=' in line:
        return line.replace(line.split('=')[-1],'file://'+os.getcwd()+'/webdir/')
    if 'ifos=' in line:
        return "ifos=['H1','L1']\n"
    if 'upload-to-gracedb=' in line:
        return 'upload-to-gracedb=True\npegasus.transfer.links=false\n'
    if 'ignore-science-segments=' in line:
        return 'ignore-science-segments=True\n'
    if 'bayesline=' in line:
        return line.replace('#','').strip()+'\n'
    if 'skyarea=' in line:
        return line.replace('#','').strip()+'\n'
    if 'types=' in line:
        return "types={'H1':'H1_HOFT_C00','L1':'L1_HOFT_C00','V1':'V1Online'}\n"
    if 'channels=' in line:
        return "channels={'H1':'H1:GDS-CALIB_STRAIN','L1':'L1:GDS-CALIB_STRAIN','V1':'V1:FAKE_h_16384Hz_4R'}\n"
    if 'margphi=' in line:
        return line.replace('#','').strip()+'\n'
    if 'amporder=' in line:
        return 'amporder=0\npsdFit=\ndifferential-buffer-limit=1000000\n'
    if 'parname-max' in line:
        return line+'distance-max=300\n'
    return line

if args.gracedb:

    os.makedirs(args.output+'/GraceDB/')
    os.chdir(args.output+'/GraceDB/')

    shutil.copy(ini_file,args.output+'/GraceDB/'+os.path.basename(ini_file)+'.bak')
    shutil.copy(ini_file,args.output+'/GraceDB/')


    with open(args.output+'/GraceDB/'+os.path.basename(ini_file),'w') as fout:
        with open(args.output+'/GraceDB/'+os.path.basename(ini_file)+'.bak','r') as fin:
            for line in fin:
                fout.write(replace_GraceDB(line))

    lalinferenceargs = [ 'lalinference_pipe'
                         , '--gid'
                         , 'T169545'
                         , '-r'
                         , './run'
                         , '-p'
                         , './daglog'
                         , args.output+'/GraceDB/'+os.path.basename(ini_file)
                         ]

    if args.condor_submit:
        lalinferenceargs.append('--condor-submit')

    subprocess.call(lalinferenceargs)

############################################################

def replace_fiducial_bbh(line):
    if 'webdir=' in line:
        return line.replace(line.split('=')[-1],os.getcwd()+'/webdir/')
    if 'baseurl=' in line:
        return line.replace(line.split('=')[-1],'file://'+os.getcwd()+'/webdir/')
    if 'fake-cache=' in line:
        return line.replace(line,"fake-cache={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}")
    if 'ignore-science-segments=' in line:
        return 'ignore-science-segments=True\n'
    if 'dataseed=' in line:
        return line.replace('#','').strip()+'\n'
    if 'disable-spin=' in line:
        return '#disable-spin=\n'
    if 'margphi=' in line:
        return '#margphi=\n'
    if 'margtime=' in line:
        return line.replace('#','').strip()+'\n'
    if 'amporder=' in line:
        return 'amporder=-1\nfref=0\n'
    if 'parname-max' in line:
        return line+'distance-max=2000\n'
    if 'deltaLogL=' in line:
        return 'deltaLogL=7\n'
    return line

if args.bbh_injection != '':

    os.makedirs(args.output+'/fiducialBBH/')
    os.chdir(args.output+'/fiducialBBH/')

    shutil.copy(args.bbh_injection,args.output+'/fiducialBBH/')
    shutil.copy(ini_file,args.output+'/fiducialBBH/'+os.path.basename(ini_file)+'.bak')
    shutil.copy(ini_file,args.output+'/fiducialBBH/')

    with open(args.output+'/fiducialBBH/'+os.path.basename(ini_file),'w') as fout:
        with open(args.output+'/fiducialBBH/'+os.path.basename(ini_file)+'.bak','r') as fin:
            for line in fin:
                fout.write(replace_fiducial_bbh(line))

    lalinferenceargs = [ 'lalinference_pipe'
                         , '-I'
                         , args.bbh_injection
                         , '-r'
                         , './run'
                         , '-p'
                         , './daglog'
                         , args.output+'/fiducialBBH/'+os.path.basename(ini_file)
                         ]

    if args.condor_submit:
        lalinferenceargs.append('--condor-submit')

    subprocess.call(lalinferenceargs)
