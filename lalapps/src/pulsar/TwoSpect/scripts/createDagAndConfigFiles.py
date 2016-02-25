#!/usr/bin/python

import argparse, sys, os, subprocess, math

def writeConfigFile(args, fmin, directorynumber):
    configFile = open('{}_{}/in/{}'.format(args.dir, args.IFO.replace(',',''), directorynumber),'w')
    configFile.write("""\
inputSFTs \"{}\"
IFO {}
fmin {}
fspan {}
t0 {}
Tobs {}
Tsft {}
SFToverlap {}
ihsfar {}
ihsfomfar {}
tmplfar {}
Pmin {}
Pmax {}
dfmin {}
dfmax {}
blksize {}
avesqrtSh {}
minTemplateLength {}
maxTemplateLength {}
outdirectory {}_{}/output/{}
skyRegionFile {}_{}/in/skygrid.{}
FFTplanFlag {}
fastchisqinv TRUE
vectorMath {}
ihsfactor {}
randSeed {}
""".format(args.inputSFTs, args.IFO, fmin, args.fspan, args.t0, args.Tobs, args.Tsft, args.SFToverlap, args.ihsfar, args.ihsfomfar, args.tmplfar, args.Pmin, args.Pmax, args.dfmin, args.dfmax, args.blksize, args.avesqrtSh, args.minTemplateL, args.maxTemplateL, args.dir, args.IFO.replace(',',''), directorynumber, args.dir, args.IFO.replace(',',''), directorynumber, args.FFTplan, args.vectorMath, args.ihsfactor, int(float(fmin)*args.Tsft + 0.5)+directorynumber))
    if args.linedetection: configFile.write('lineDetection {}\n'.format(args.linedetection))
    if args.keepTopIHS: configFile.write('keepOnlyTopNumIHS {}\n'.format(args.keepTopIHS))
    if args.markBadSFTs: configFile.write('markBadSFTs TRUE\n')
    if args.assumeNScosi: configFile.write('assumeNScosi {}\n'.format(args.assumeNScosi))
    if args.assumeNSpsi: configFile.write('assumeNSpsi {}\n'.format(args.assumeNSpsi))
    if args.cosiSignCoherent: configFile.write('cosiSignCoherent {}\n'.format(args.cosiSignCoherent))
    if args.ulonly: configFile.write('IHSonly TRUE\n')
    configFile.close()

def writeJobToDagFile(args, directorynumber):
    dagfile = open('{}_{}/dag'.format(args.dir, args.IFO.replace(',','')), 'a')
    dagfile.write('JOB A{} {}_{}/condor\n'.format(directorynumber, args.dir, args.IFO.replace(',','')))
    dagfile.write('VARS A{} PID=\"{}\"\n'.format(directorynumber, directorynumber))
    dagfile.close()
    

parser = argparse.ArgumentParser(description='Script that sets up the condor/dag file for TwoSpect analysis',fromfile_prefix_chars='@')
parser.add_argument('--dir', required=True, type=str, help='Output directory path. Detector name will be appended')
parser.add_argument('--logfile', required=True, type=str, help='Output path/filename of log file, detector name appended')
parser.add_argument('--program', required=True, type=str, help='Path and file name of executable')
parser.add_argument('--accountingGroup', required=True, type=str, help='Label of accounting group, ex: ligo.sim.s6.cw.allskybinary.twospect')
parser.add_argument('--memRequest', type=int, default=1550, help='Amount of memory to request of the nodes in MB')
parser.add_argument('--fstart', required=True, type=float, help='Frequency in Hz to start making config files')
parser.add_argument('--fspan', required=True, type=float, help='Frequency span in Hz of each band being searched')
parser.add_argument('--numberBands', type=int, default=1, help='Number of bands to produce >= 1')
parser.add_argument('--IFO', required=True, type=str, help='Interferometer to use (specify multiple as CSV)')
parser.add_argument('--inputSFTs', required=True, type=str, help='Path and file name of inputSFTs (can be glob-style)')
parser.add_argument('--t0', required=True, type=int, help='Start time of search')
parser.add_argument('--Tobs', type=int, default=40551300, help='Observation time of the search in seconds (%(default)s)')
parser.add_argument('--Tsft', type=int, default=1800, help='Coherence length of each SFT (%(default)s)')
parser.add_argument('--SFToverlap', type=int, default=900, help='Overlap of each SFT in seconds (%(default)s)')
parser.add_argument('--Pmin', type=float, default=7200, help='Minimum period of search (%(default)s)')
parser.add_argument('--Pmax', type=float, default=8110260, help='Maximum period of search (%(default)s)')
parser.add_argument('--dfmin', type=float, default=0.000277, help='Minimum modulation depth of search in Hz (%(default)s)')
parser.add_argument('--dfmax', type=float, default=0.1, help='Maximum modulation depth of search (%(default)s)')
parser.add_argument('--ihsfar', type=float, default=1.0e-14, help='IHS false alarm rate (%(default)s)')
parser.add_argument('--ihsfomfar', type=float, default=1.0, help='IHS figure of merit false alarm rate (%(default)s)')
parser.add_argument('--tmplfar', type=float, default=1.0e-18, help='Template statistic false alarm rate (%(default)s)')
parser.add_argument('--blksize', type=int, default=101, help='Block size of running median (%(default)s)')
parser.add_argument('--minTemplateL', type=int, default=1, help='Minimum template length in pixels (%(default)s)')
parser.add_argument('--maxTemplateL', type=int, default=500, help='Maximum template length in pixels (%(default)s)')
parser.add_argument('--avesqrtSh', type=float, default=1.0e-22, help='Expected noise background of SFTs (%(default)s)')
parser.add_argument('--FFTplan', type=int, default=3, help='Plan flag for the FFTs: 0, 1, 2, or 3 (%(default)s)')
parser.add_argument('--maxSkyLocPerJob', type=int, default=200, help='Maximum sky locations per job (%(default)s)')
parser.add_argument('--linedetection', type=float, help='Line detection threshold setting')
parser.add_argument('--keepTopIHS', type=int, help='Keep only the top N IHS candidates')
parser.add_argument('--ihsfactor', type=int, default=5, help='IHS folding factor')
parser.add_argument('--fastchisqinv', action='store_true', help='Use fast chi squared inversion methods')
parser.add_argument('--vectorMath', type=int, choices=[0, 1, 2], default=1, help='Use vector math 0 = none, 1 = SSE, 2 = AVX')
parser.add_argument('--markBadSFTs', action='store_true', help='Mark the non-Gaussian SFTs and remove them')
parser.add_argument('--assumeNScosi', type=float, help='Assume a particular cos(iota) for coherent SFT summing')
parser.add_argument('--assumeNSpsi', type=float, help='Assume a particular psi for coherent SFT summing')
parser.add_argument('--cosiSignCoherent', type=int, choices=[-1, 0, 1], default=0, help='For coherent analysis assume [-1,1] values (0), [0,1] values (1), or [-1,0] values (-1) for cosi (Note: unused when assumeNScosi is specified)')
parser.add_argument('--scoX1', action='store_true', help='Perform the restricted, Sco X-1 search')
parser.add_argument('--ulonly', action='store_true', help='Only do upper limits')
parser.add_argument('--skygridsetupProgram', required=True, type=str, help='Path to and program name of the skygridsetup program')
args = parser.parse_args()

os.mkdir('{}_{}'.format(args.dir, args.IFO.replace(',','')))
os.mkdir('{}_{}/in'.format(args.dir, args.IFO.replace(',','')))
os.mkdir('{}_{}/err'.format(args.dir, args.IFO.replace(',','')))
os.mkdir('{}_{}/output'.format(args.dir, args.IFO.replace(',','')))

condorfile = open('{}_{}/condor'.format(args.dir, args.IFO.replace(',','')),'w')
condorfile.write("""\
universe=vanilla
executable={}
arguments=@{}_{}/in/$(PID)
input=/dev/null
output=/dev/null
error={}_{}/err/err.$(PID)
log={}
request_memory={}
notification=Never
accounting_group={}
requirements= simd_max=="avx2"
queue
""".format(args.program, args.dir, args.IFO.replace(',',''), args.dir, args.IFO.replace(',',''), args.logfile, args.memRequest, args.accountingGroup))
condorfile.close()

directorynumber = 0

dagfile = open('{}_{}/dag'.format(args.dir, args.IFO.replace(',','')),'w')
for ii in range(0, args.numberBands):
    fmin = '{:4.3f}'.format(args.fstart+ii*args.fspan)

    if os.path.exists('{}_{}/skygrid-{}Hz-{}Hz'.format(args.dir, args.IFO.replace(',',''), fmin, args.fspan)): exit()

    if not args.scoX1:
        IFOList = args.IFO.split(',')
        subprocess.check_call(['{}'.format(args.skygridsetupProgram), '--fmin={}'.format(fmin), '--fspan={}'.format(args.fspan), '--Tsft={}'.format(args.Tsft), '--IFO={}'.format(IFOList[0]), '--SFToverlap={}'.format(args.SFToverlap), '--t0={}'.format(args.t0), '--Tobs={}'.format(args.Tobs), '--outfilename={}_{}/skygrid-{}Hz-{}Hz.dat'.format(args.dir, args.IFO.replace(',',''), fmin, args.fspan)])
    else:
         skyRegionFileName = '{}_{}/skygrid-{}Hz-{}Hz.dat'.format(args.dir, args.IFO.replace(',',''), fmin, args.fspan)
         skyRegionFile = open(skyRegionFileName, 'w')
         skyRegionFile.write('4.275699238500 -0.272973858335')
         skyRegionFile.close()

         scoX1P = 68023.70
         scoX1asini = 1.44
         scoX1Perr = 0.0432
         scoX1asinierr = 0.18
         args.Pmin = args.Tobs/(int(args.Tobs/scoX1P+0.5)+1.0)
         args.Pmax = args.Tobs/(int(args.Tobs/scoX1P+0.5)-1.0)
         args.dfmin = 2*math.pi*float(fmin)*(scoX1asini-3.0*scoX1asinierr)/(scoX1P+3.0*scoX1Perr)
         args.dfmax = 2*math.pi*(float(fmin)+args.fspan)*(scoX1asini+3.0*scoX1asinierr)/(scoX1P-3.0*scoX1Perr)

    with open('{}_{}/skygrid-{}Hz-{}Hz.dat'.format(args.dir, args.IFO.replace(',',''), fmin, args.fspan)) as skygridfile:
        numberSkyLocations = 0
        firstPointSet = 0
        startbin = int(float(fmin)*args.Tsft + 0.5)
        
        for line in skygridfile:
            if not firstPointSet and line!='':
                skyRegionFileName = '{}_{}/in/skygrid.{}'.format(args.dir, args.IFO.replace(',',''), directorynumber)
                skyRegionFile = open(skyRegionFileName, 'w')
                skyRegionFile.write(line)
                numberSkyLocations += 1
                firstPointSet = 1
            elif firstPointSet and numberSkyLocations<args.maxSkyLocPerJob and line!='':
                skyRegionFile.write(line)
                numberSkyLocations += 1

            if firstPointSet and numberSkyLocations==args.maxSkyLocPerJob:
                skyRegionFile.close()
                writeConfigFile(args, fmin, directorynumber)
                writeJobToDagFile(args, directorynumber)
                directorynumber += 1
                numberSkyLocations = 0
                firstPointSet = 0

        if not skyRegionFile.closed:
            skyRegionFile.close()
            writeConfigFile(args, fmin, directorynumber)
            writeJobToDagFile(args, directorynumber)
            directorynumber += 1
            numberSkyLocations = 0
            firstPointSet = 0
                
    
