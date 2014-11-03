#!/usr/bin/python

import argparse, sys, os, shutil, subprocess, math, random

parser = argparse.ArgumentParser(description='Script that generates configuration files for TwoSpect Monte Carlo injections and runs TwoSpect')
parser.add_argument('--twospectExe', required=True, type=str, help='Path to TwoSpect executable')
parser.add_argument('--dir', required=True, type=str, help='Base directory where the job subdirectories are located')
parser.add_argument('--jobnum', required=True, type=int, help='Job number, also the base directory/subdirectory path')
parser.add_argument('--IFO', required=True, type=str, help='Interferometer to use (specify multiple as CSV)')
parser.add_argument('--t0', required=True, type=int, help='Start time of search')
parser.add_argument('--fmin', required=True, type=float, help='Minimum frequency of injection (Hz)')
parser.add_argument('--inputSFTs', type=str, help='Input SFT files')
parser.add_argument('--gaussianNoiseWithSFTgaps', action='store_true', help='Input SFT files')
parser.add_argument('--injPol', type=int, choices=[0, 1, 2], default=2, help='Polarizations, 0 = linear, 1 = random, 2 = circular (%(default)s)')
parser.add_argument('--minEcc', type=float, default=0, help='Minimum of eccentricity of orbit (%(default)s)')
parser.add_argument('--maxEcc', type=float, default=0, help='Maximum of eccentricity of orbit (%(default)s)')
parser.add_argument('--eccDist', type=int, choices=[-1, 0, 1], default=0, help='Ecc. distribution, 0 = linear, 1 = log, -1 = inv. log (%(default)s)')
parser.add_argument('--minSpindown', type=float, default=0, help='Minimum spindown of source (%(default)s)')
parser.add_argument('--maxSpindown', type=float, default=0, help='Maximum spindown of source (%(default)s)')
parser.add_argument('--spindownDist', type=int, choices=[0], default=0, help='Spindown distribution, 0 = linear (%(default)s)')
parser.add_argument('--Tobs', type=int, default=40551300, help='Observation time of the search in seconds (%(default)s)')
parser.add_argument('--Tsft', type=int, default=1800, help='Coherence length of each SFT (%(default)s)')
parser.add_argument('--SFToverlap', type=int, default=900, help='Overlap of each SFT in seconds (%(default)s)')
parser.add_argument('--fspan', type=float, default=0.25, help='Frequency span of injection band (%(default)s Hz)')
parser.add_argument('--h0min', type=float, help='Minimum strain value to inject')
parser.add_argument('--h0max', type=float, help='Maximum strain value to inject')
parser.add_argument('--h0dist', type=int, choices=[-1, 0, 1], default=1, help='Distribution of h0, 0 = linear, 1 = log, -1 = inv. log (%(default)s)')
parser.add_argument('--h0val', type=float, help='Constant strain value to inject')
parser.add_argument('--skylocations', type=int, default=0, help='Number of sky locations to search, 0 exact location (%(default)s)')
parser.add_argument('--timestampsFile', type=str, help='File containing start times of SFTs to produce (specify multiple as CSV)')
parser.add_argument('--segmentFile', type=str, help='File containing <start end> of segments to produce SFTs (specify multiple as CSV)')
parser.add_argument('--injskyra', type=float, help='Right ascension of injection')
parser.add_argument('--injskydec', type=float, help='Declination of injection')
parser.add_argument('--injPmin', type=float, default=7200, help='Minimum period of injection')
parser.add_argument('--injPmax', type=float, default=8110260, help='Maximum period of injection')
parser.add_argument('--periodDist', type=int, choices=[-1, 0, 1], default=0, help='Orbit period dist. 0 = linear, 1 = log, -1 = inv. log (%(default)s)')
parser.add_argument('--injDfmin', type=float, default=0.0, help='Minimum modulation depth of injection in Hz (%(default)s)')
parser.add_argument('--injDfmax', type=float, default=0.1, help='Maximum modulation depth of injection (%(default)s)')
parser.add_argument('--injDfExpAllow', type=float, default=1.0, help='(Not yet implemented)')
parser.add_argument('--ihsfactor', type=int, default=5, help='IHS folding factor (%(default)s)')
parser.add_argument('--seed', type=int, default=0, help='Seed value + jobnum for producing injections and Gaussian noise (%(default)s + jobnum)')
parser.add_argument('--scox1', action='store_true', help='Inject Sco X-1 signals')
parser.add_argument('--weightedIHS', action='store_true', help='Use weighted IHS statistic')
parser.add_argument('--ulonly', action='store_true', help='Only produce ULs from IHS statistic, no follow up')
parser.add_argument('--templateTest', action='store_true', help='Brute force template test')
parser.add_argument('--ihsfar', type=float, default=1.0, help='IHS false alarm rate (%(default)s)')
parser.add_argument('--ihsfomfar', type=float, default=1.0, help='IHS figure of merit false alarm rate (%(default)s)')
parser.add_argument('--tmplfar', type=float, default=1.0, help='Template statistic false alarm rate (%(default)s)')
parser.add_argument('--tmplLength', type=int, default=500, help='Maximum length of a template (%(default)s)')
parser.add_argument('--markBadSFTs', action='store_true', help='Mark and remove bad SFTs')
parser.add_argument('--keepOnlyTopNumIHS', type=int, help='Keep only the top N number of IHS outliers')
parser.add_argument('--useCorrectNScosi', action='store_true', help='Use the correct NS cosi orientation for coherent analysis')
parser.add_argument('--useCorrectNSpsi', action='store_true', help='Use the correct NS psi orientation for coherent analysis')
args = parser.parse_args()

IFOList = args.IFO.split(',')
numIFOs = len(IFOList)
numberTimestampFiles = 0
numberSegmentFiles = 0
if args.timestampsFile:
    timestampFileList = args.timestampsFile.split(',')
    numberTimestampFiles = len(timestampFileList)
if args.segmentFile:
    segmentFileList = args.segmentFile.split(',')
    numberSegmentFiles = len(segmentFileList)

if not args.inputSFTs and args.gaussianNoiseWithSFTgaps: sys.exit('Must provide SFT file when gaussianNoiseWithSFTgaps is specifed')
if (args.inputSFTs and (numberTimestampFiles>0 or numberSegmentFiles>0)) or (numberTimestampFiles>0 and numberSegmentFiles>0): sys.exit('Only choose sftFile OR sftFile and gaussianNoiseWithSFTgaps OR timestampsfile OR segmentfile')
if numberTimestampFiles>0 and numIFOs!=numberTimestampFiles: sys.exit('Number of IFOs and timestamp files must be the same')
if numberSegmentFiles>0 and numIFOs!=numberSegmentFiles: sys.exit('Number of IFOs and segment files must be the same')
if args.minEcc<0: sys.exit('Minimum eccentricity must be 0 or larger')
if args.maxEcc>=1: sys.exit('Minimum eccentricity must be 0 or larger')
if (args.h0min and not args.h0max) or (args.h0max and not args.h0min): sys.exit('Must specify both h0min and h0max')
if (args.injskyra and not args.injskydec) or (args.injskydec and not args.injskyra): sys.exit('Must specify both injskyra and injskydec')
if args.h0val and (args.h0min or args.h0max): sys.exit('h0val cannot be specified with h0min or h0max')
if args.skylocations<0: sys.exit('skylocations must be 0 or greater')
if args.ihsfactor<1: sys.exit('ihsfactor must be 1 or greater')
if args.injDfExpAllow<1.0: sys.exit('injDfExpAllow must be 1 or greater')
if args.ihsfar>1.0: sys.exit('ihsfar must be less than or equal to 1')
if args.ihsfomfar>1.0: sys.exit('ihsfomfar must be less than or equal to 1')
if args.tmplfar>1.0: sys.exit('tmplfar must be less than or equal to 1')
if args.tmplLength<1: sys.exit('tmplLength must be 1 or greater')

scoX1P = 68023.70
scoX1asini = 1.44

scriptPID = os.getpid()
os.mkdir('/local/user/egoetz/{}'.format(scriptPID))

ifokey = []
skygridfile = ''
skygridcreated = 0
for ifoval in IFOList:
    if ifoval == 'LHO' or ifoval == 'H1': ifokey.append('H1')
    elif ifoval == 'LLO' or ifoval == 'L1': ifokey.append('L1')
    elif ifoval == 'Virgo' or ifoval == 'V1': ifokey.append('V1')
    elif ifoval == 'H2': ifokey.append('H2')
    elif ifoval == 'H2r': ifokey.append('H2r')
    else: sys.exit('Interferometer '+ifoval+' not recognized!')
    if args.skylocations>=1 and skygridcreated!=1:
        resetH2r = 0
        if ifokey[-1]=='H2r':
            ifokey[-1] = 'H2'
            resetH2r = 1
        skygridfile = '/local/user/egoetz/{}/skygrid-{}.dat'.format(scriptPID, ifokey[-1])
        subprocess.check_call(['/atlas/user/atlas3/egoetz/lalsuite-master/lalapps/src/pulsar/TwoSpect/skygridsetup','--fmin={} --fspan={} --IFO={} --Tcoh={} --SFToverlap={} --t0={} --Tobs={} --v2 --outfilename={}'.format(args.fmin, args.fspan, ifokey[-1], args.Tsft, args.SFToverlap, args.t0, args.dur, skygridfile)])
        if resetH2r==1: ifokey[-1] = 'H2r'
        skygridcreated = 1
        
Pmin = args.injPmin
Pmax = args.injPmax
dfmin = args.injDfmin
dfmax = args.injDfmax
if args.scox1:
    Pmin = args.Tobs/(round(args.Tobs/scoX1P)+1.0);
    Pmax = args.Tobs/(round(args.Tobs/scoX1P)-1.0);
    dfmin = 2*math.pi*args.fmin*(scoX1asini-3.0*0.18)/(scoX1P+3.0*.0432);
    dfmax = 2*math.pi*(args.fmin+args.fspan)*(scoX1asini+3.0*0.18)/(scoX1P-3.0*.0432);

random.seed(args.seed+args.jobnum)

for ii in range(0, 10):
    h0 = 0
    if args.h0val: h0 = args.h0val
    else:
        if args.h0dist==0: h0 = (args.h0max-args.h0min)*random.random() + args.h0min
        elif args.h0dist==1: h0 = 10**((math.log10(args.h0max)-math.log10(args.h0min))*random.random()) * args.h0min
        else: h0 = (args.h0max + args.h0min) - 10**((math.log10(args.h0max)-math.log10(args.h0min))*random.random()) * args.h0min
    
    psi = math.pi*random.random()
    phi0 = 2*math.pi*random.random()

    alpha = 0
    delta = 0
    if not args.injskyra and not args.injskydec:
        alpha = 2*math.pi*random.random()
        delta = math.acos(random.uniform(-1,1))-0.5*math.pi
    else:
        alpha = args.injskyra
        delta = args.injskydec

    f0 = args.fmin + args.fspan*random.random()

    P = 0
    df = 0
    asini = 0
    ecc = 0
    argp = 0
    f1dot = 0.0
    if args.scox1:
        alpha = alpha = 4.275699238500
        delta = -0.272973858335
        P = scoX1P + random.gauss(0, 0.0432)
        asini = scoX1asini + random.gauss(0, 0.18)
        df = 2*math.pi*f0*asini/P
    else:
        if args.periodDist==0: P = (args.injPmax - args.injPmin)*random.random() + args.injPmin
        elif args.periodDist==1: P = 10**((math.log10(args.injPmax)-math.log10(args.injPmin))*random.random()) * args.injPmin
        else: P = (args.injPmax + args.injPmin) - 10**((math.log10(args.injPmax)-math.log10(args.injPmin))*random.random()) * args.injPmin
        
        df = (args.injDfmax - args.injDfmin)*random.random() + args.injDfmin
        while df-0.5/args.Tsft<1e-6 or df > args.injDfExpAllow*P/(2*args.Tsft*args.Tsft): df = (args.injDfmax - args.injDfmin)*random.random() + args.injDfmin

        asini = df*P/(2*math.pi*f0)

        if args.maxEcc>0:
            if args.eccDist==0: ecc = (args.maxEcc - args.minEcc)*random.random() + args.minEcc
            elif args.eccDist==1: ecc = 10**((math.log10(args.maxEcc)-math.log10(args.minEcc))*random.random()) * args.minEcc
            else: ecc = (args.maxEcc + args.minEcc) - 10**((math.log10(args.maxEcc)-math.log10(args.minEcc))*random.random()) * args.minEcc
            argp = 2*math.pi*random.random()

        if args.minSpindown!=0 or args.maxSpindown!=0: f1dot = (args.maxSpindown - args.minSpindown)*random.random() + args.minSpindown

    cosi = 1
    if args.injPol==0: cosi = 0
    elif args.injPol==1: cosi = random.uniform(-1,1)
    
    mfdconfig = open('/local/user/egoetz/{}/mfdconfig'.format(scriptPID),'w')
    mfdconfig.write("""\
Alpha {}
Delta {}
h0 {}
cosi {}
psi {}
phi0 {}
Freq {}
orbitasini {}
orbitEcc {}
orbitTpSSB 900000000
orbitPeriod {}
orbitArgp {}
f1dot {}
refTime 900000000""".format(alpha,delta,h0,cosi,psi,phi0,f0,asini,ecc,P,argp,f1dot))
    mfdconfig.close()

    injections = open(args.dir+'/'+str(args.jobnum)+'/injections.dat','a')
    injections.write('{} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(alpha, delta, h0, cosi, psi, phi0, f0, asini, ecc, P, argp, f1dot, df))
    injections.close()

    twospectconfig = open('/local/user/egoetz/'+str(scriptPID)+'/twospectconfig','w')
    twospectconfig.write("""\
fmin {}
fspan {}
t0 {}
Tobs {}
Tcoh {}
SFToverlap {}
ihsfar {}
ihsfomfar {}
tmplfar {}
Pmin {}
Pmax {}
dfmin {}
dfmax {}
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 1
maxTemplateLength {}
outdirectory {}/{}
ephemEarth /home/egoetz/opt/lscsoft-master/share/lalpulsar/earth00-19-DE405.dat.gz
ephemSun /home/egoetz/opt/lscsoft-master/share/lalpulsar/sun00-19-DE405.dat.gz
FFTplanFlag 1
fastchisqinv TRUE
useSSE TRUE
outfilename logfile_{}.txt
ULfilename uls_{}.dat
configCopy input_copy_{}.conf
injectionSources @/local/user/egoetz/{}/mfdconfig
ihsfactor {}
""".format(args.fmin, args.fspan, args.t0, args.Tobs, args.Tsft, args.SFToverlap, args.ihsfar, args.ihsfomfar, args.tmplfar, Pmin, Pmax, dfmin, dfmax, args.tmplLength, args.dir, args.jobnum, ii, ii, ii, scriptPID, args.ihsfactor))

    if not args.inputSFTs: twospectconfig.write('injRandSeed '+str(random.randint(1, 1000000))+'\n')
    if args.markBadSFTs: twospectconfig.write('markBadSFTs TRUE\n')
    if args.ulonly: twospectconfig.write('IHSonly TRUE\n')
    if args.keepOnlyTopNumIHS: twospectconfig.write('keepOnlyTopNumIHS {}\n'.format(args.keepOnlyTopNumIHS))

    if args.skylocations==0: twospectconfig.write('skyRegion ({},{})\n'.format(alpha, delta))
    elif args.skylocations>0:
        RAvalues = []
        DECvalues = []
        distances = []
        skyfile = open(skygridfile,'r')
        for line in skyfile:
            values = map(float, line.split())
            RAvalues.append(values[0])
            DECvalues.append(values[1])
            distances.append(math.acos(math.sin(abs(values[1]-0.5*math.pi))*math.sin(abs(args.injskydec-0.5*math.pi))*math.cos(values[0]-args.injskyra)+math.cos(abs(values[1]-0.5*math.pi))*math.cos(abs(args.injskydec-0.5*math.pi))))
        skyfile.close()
        sortedindex = [jj[0] for jj in sorted(enumerate(distances), key=lambda x:x[1])]
        skyfile2 = open('/local/user/egoetz/{}/skygrid2.dat\n'.format(scriptPID),'w')
        for jj in range(0, args.skylocations-1): skyfile2.write('{} {}\n'.format(RAvalues[sortedindex[jj]], DECvalues[sortedindex[jj]]))
        skyfile2.close()
        twospectconfig.write('skyRegionFile /local/user/egoetz/{}/skygrid2.dat\n'.format(scriptPID))

    if args.templateTest: twospectconfig.writelines('bruteForceTemplateTest TRUE\n','templateTestF {}\n'.format(f0), 'templateTestP {}\n'.format(P), 'templateTestDf {}\n'.format(df))

    if args.inputSFTs: twospectconfig.write('inputSFTs {}\n'.format(args.inputSFTs))

    if args.gaussianNoiseWithSFTgaps: twospectconfig.write('gaussianNoiseWithSFTgaps TRUE\n')

    if args.useCorrectNScosi: twospectconfig.write('assumeNScosi {}\n'.format(cosi))
    if args.useCorrectNSpsi: twospectconfig.write('assumeNSpsi {}\n'.format(psi))

    twospectconfig.write('IFO ')
    for jj in range(0, numIFOs):
        twospectconfig.write(IFOList[jj])
        if jj<numIFOs-1: twospectconfig.write(',')
    twospectconfig.write('\n')

    if numberTimestampFiles>0:
        twospectconfig.write('timestampsFile ')
        for jj in range(0, numIFOs):
            twospectconfig.write(timestampsFileList[jj])
            if jj<numIFOs-1: twospectconfig.write(',')
    elif numberSegmentFiles>0:
        twospectconfig.write('segmentFile ')
        for jj in range(0, numIFOs):
            twospectconfig.write(segmentFileList[jj])
            if jj<numIFOs-1: twospectconfig.write(',')
    
    twospectconfig.close()

    subprocess.check_call([args.twospectExe, '@/local/user/egoetz/{}/twospectconfig'.format(scriptPID)])

shutil.rmtree('/local/user/egoetz/{}'.format(scriptPID))

   
