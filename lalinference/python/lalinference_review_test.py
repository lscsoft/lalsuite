import argparse
import shutil
import distutils.spawn
import os
import subprocess
import socket
import glob
import lalinference
import ConfigParser

prefix=''
try:
    prefix=os.environ['LALINFERENCE_DATADIR']
except KeyError:
    prefix=''

parser = argparse.ArgumentParser(description="Runs review tests of lalinference using lalinference_pipe_example.ini.")

parser.add_argument('-i','--ini_file', type=str, nargs='?',
                    default=os.path.join(prefix,'lalinference_pipe_example.ini'),
                    help='lalinference_pipe ini file to process.')

parser.add_argument('--bns-injection', type=str, nargs='?',
                    default=False,
                    const=os.path.join(prefix,'fiducialBNS.xml'),
                    help='injection file for fiducial BNS analysis.')

parser.add_argument('--gracedb', action='store_true',
                    default=False,
                    help='Runs the analysis for the GraceDB test event T169545.')

parser.add_argument('--analytic-tests', action='store_true',
                    default=False,
                    help='Run on the unmodal and bimodial Gaussian and Rosenbrock test functions.')

parser.add_argument('--analytic-csv-dir', type=str,
                    default=prefix,
                    help='Directory containing the CSVs describing the analytic tests.')

parser.add_argument('--pptest', action='store_true',
                    default=False,
                    help='Runs a P-P analysis.')

parser.add_argument('--bbh-injection', type=str, nargs='?',
                    default=False,
                    const=os.path.join(prefix,'fiducialBBH.xml'),
                    help='injection file for fiducial BBH analysis.')

parser.add_argument('-e','--engine', type=str, nargs='?',
                    default='lalinferencemcmc,lalinferencenest',
                    help='lalinference engine to run with.')

parser.add_argument('-o','--output', type=str, nargs='?',
                    default=None,
                    help='output directory.')

parser.add_argument('--condor-submit', action='store_true',
                    default=False,
                    help='submits the test suite.')

args = parser.parse_args()


if 'UNCLEAN' in lalinference.InferenceVCSInfo.vcsId:
    default_outputdir=os.getenv('HOME')+'/lalinference_testrun/'+lalinference.InferenceVCSInfo.vcsId+'_UNCLEAN/'+args.engine+'/'
else:
    default_outputdir=os.getenv('HOME')+'/lalinference_testrun/'+lalinference.InferenceVCSInfo.vcsId+'/'+args.engine.replace(',','_')+'/'

if args.output == None:
    web_outputdir=default_outputdir
    args.output=os.path.abspath(default_outputdir)
else:
    web_outputdir=os.path.abspath(args.output)
    args.output=os.path.abspath(args.output)

os.makedirs(args.output)

if args.bns_injection:
    args.bns_injection=os.path.abspath(args.bns_injection)
if args.bbh_injection:
    args.bbh_injection=os.path.abspath(args.bbh_injection)

analytic_csv_dir = os.path.abspath(args.analytic_csv_dir)
if args.analytic_tests:
    csvs = glob.glob(analytic_csv_dir+"/*csv")
    for csv in csvs:
        shutil.copy(csv,args.output+'/'+os.path.basename(csv))

os.chdir(args.output)

lalinf_prefix=''
try:
    lalinf_prefix=os.environ['LALINFERENCE_PREFIX']
except KeyError:
    print 'LALINFERENCE_PREFIX variable not defined, could not find LALInference installation.'
    sys.exit()

def init_ini_file(file=args.ini_file):
    cp=ConfigParser.SafeConfigParser()
    fp=open(file)
    cp.optionxform = str
    cp.readfp(fp)
    fp.close()

    cp.set('condor','lalsuite-install',lalinf_prefix)
    cp.set('analysis','engine',args.engine)
    cp.remove_option('analysis','nparallel')

    return cp

############################################################

def set_fiducial_bns(cp):

    cp.set('lalinference','fake-cache',"{'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}")
    cp.set('analysis','dataseed','1234')
    cp.set('engine','0noise','')

    cp.set('paths','webdir',web_outputdir+'/fiducialBNS/webdir/')
    cp.set('lalinference','flow',"{'H1':40,'L1':40,'V1':40}")
    cp.set('engine','approx','SEOBNRv4_ROMpseudoFourPN')
    cp.set('resultspage','deltaLogP','5')
    cp.set('engine','comp-max','3.5')
    cp.set('engine','comp-min','0.8')

    cp.set('engine','neff','500')
    cp.set('engine','nlive','512')

    return cp

if args.bns_injection:

    os.makedirs(args.output+'/fiducialBNS/')
    os.chdir(args.output+'/fiducialBNS/')

    shutil.copy(args.bns_injection,args.output+'/fiducialBNS/')
    BNS_ini_file=os.path.join(args.output,'fiducialBNS','BNS.ini')

    cpBNS=set_fiducial_bns(init_ini_file())
    with open(BNS_ini_file,'w') as cpfile:
        cpBNS.write(cpfile)

    lalinferenceargs = [ 'lalinference_pipe'
    		     , '-I'
    		     , args.bns_injection
    		     , '-r'
    		     , './run'
    		     , '-p'
    		     , './daglog'
    		     , BNS_ini_file ]

    if args.condor_submit:
        lalinferenceargs.append('--condor-submit')

    subprocess.call(lalinferenceargs)

###################################################

def set_GraceDB(cp):

    cp.set('analysis','upload-to-gracedb','True')
    cp.set('paths','webdir',web_outputdir+'/GraceDB/webdir/')

    cp.set('datafind','types',"{'H1':'H1_HOFT_C00','L1':'L1_HOFT_C00','V1':'V1Online'}")
    cp.set('data','channels',"{'H1':'H1:GDS-CALIB_STRAIN','L1':'L1:GDS-CALIB_STRAIN','V1':'V1:FAKE_h_16384Hz_4R'}")

    cp.set('engine','distance-max','300')

    return cp

if args.gracedb:

    os.makedirs(args.output+'/GraceDB/')
    os.chdir(args.output+'/GraceDB/')

    GDB_ini_file=os.path.join(args.output,'GraceDB','GDB.ini')

    cpGDB=set_GraceDB(init_ini_file())
    with open(GDB_ini_file,'w') as cpfile:
        cpGDB.write(cpfile)

    lalinferenceargs = [ 'lalinference_pipe'
                         , '--gid'
                         , 'T169545'
                         , '-r'
                         , './run'
                         , '-p'
                         , './daglog'
                         , GDB_ini_file ]

    if args.condor_submit:
        lalinferenceargs.append('--condor-submit')

    subprocess.call(lalinferenceargs)

############################################################

def set_fiducial_bbh(cp):

    cp.set('lalinference','fake-cache',"{'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}")
    cp.set('analysis','dataseed','1234')
    cp.set('engine','0noise','')

    cp.set('paths','webdir',web_outputdir+'/fiducialBBH/webdir/')
    cp.set('lalinference','flow',"{'H1':40,'L1':40,'V1':40}")
    cp.set('engine','approx','IMRPhenomPv2pseudoFourPN')
    cp.set('analysis','roq','True')
    cp.remove_option('engine','disable-spin')
    cp.set('resultspage','deltaLogP','6')
    cp.set('engine','distance-max','2000')

    cp.set('engine','neff','500')
    cp.set('engine','nlive','512')

    return cp

if args.bbh_injection:

    os.makedirs(args.output+'/fiducialBBH/')
    os.chdir(args.output+'/fiducialBBH/')

    shutil.copy(args.bbh_injection,args.output+'/fiducialBBH/')
    BBH_ini_file=os.path.join(args.output,'fiducialBBH','BBH.ini')

    cpBBH=set_fiducial_bbh(init_ini_file())
    with open(BBH_ini_file,'w') as cpfile:
        cpBBH.write(cpfile)

    lalinferenceargs = [ 'lalinference_pipe'
                         , '-I'
                         , args.bbh_injection
                         , '-r'
                         , './run'
                         , '-p'
                         , './daglog'
                         , BBH_ini_file ]

    if args.condor_submit:
        lalinferenceargs.append('--condor-submit')

    subprocess.call(lalinferenceargs)

############################################################
def set_analytic_test(cp, test_func):

    cp.set('paths','webdir',web_outputdir+'/'+test_func+'/webdir/')
    cp.set('lalinference','fake-cache',"{'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}")
    cp.set('analysis','dataseed','1234')

    cp.set('input','analyse-all-time','True')
    cp.set('input','gps-start-time','0')
    cp.set('input','gps-end-time','2')
    cp.set('input','segment-overlap','0')
    cp.set('input','psd-length','1')

    cp.set('engine','seglen','1')
    cp.set('engine','approx','SpinTaylorT4')
    cp.set('engine',test_func,'')
    cp.set('engine','neff','10000')
    cp.set('engine','nlive','512')

    cp.set('resultspage','deltaLogP','7')
    if test_func != "rosenbrockLikelihood":
        csv = args.output+'/'+test_func+'_means.csv'
        if os.path.isfile(csv):
            cp.set('resultspage','meanVectors',csv)
        csv = args.output+'/'+'test_correlation_matrix.csv'
        if os.path.isfile(csv):
            cp.set('resultspage','covarianceMatrix',csv)

    return cp

if args.analytic_tests:
    test_funcs = ['correlatedGaussianLikelihood', 'bimodalGaussianLikelihood', 'rosenbrockLikelihood']
    for test_func in test_funcs:
        os.makedirs(args.output+'/' + test_func + '/')
        os.chdir(args.output+'/' + test_func + '/')


        shutil.copy(args.bbh_injection,args.output+'/'+test_func+'/')
        analytic_ini_file=os.path.join(args.output,test_func,'analytic.ini')

        cpanalytic=set_analytic_test(init_ini_file(), test_func)
        with open(analytic_ini_file,'w') as cpfile:
            cpanalytic.write(cpfile)

        lalinferenceargs = [ 'lalinference_pipe'
                             , '-r'
                             , './run'
                             , '-p'
                             , './daglog'
                             , analytic_ini_file ]

        if args.condor_submit:
            lalinferenceargs.append('--condor-submit')

        subprocess.call(lalinferenceargs)

############################################################

def set_pptest(cp):

    cp.set('paths','webdir',web_outputdir+'/pptest/webdir/')
    cp.set('lalinference','fake-cache',"{'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO','V1':'LALSimAdVirgo'}")
    cp.set('analysis','dataseed','1234')

    cp.remove_option('engine','margphi')
    cp.set('engine','margtime','')
    cp.set('engine','amporder','-1')
    cp.set('engine','fref','0')
    cp.set('engine','distance-max','2000')

    cp.set('resultspage','deltaLogP','7')

    return cp

if args.pptest:

    os.makedirs(args.output+'/pptest/')
    os.chdir(args.output+'/pptest/')

    pptest_ini_file=os.path.join(args.output,'pptest','pptest.ini')

    cppptest=set_pptest(init_ini_file())
    with open(pptest_ini_file,'w') as cpfile:
        cppptest.write(cpfile)

    lalinferenceargs = [ 'lalinference_pp_pipe'
                         , '-r'
                         , './run'
                         , '-N'
                         , '100'
                         , pptest_ini_file ]

    subprocess.call(lalinferenceargs)

    if args.condor_submit:
        condor_submit_dag = ['condor_submit_dag',args.output+'/pptest/run/priortest.dag']
        subprocess.call(condor_submit_dag)
