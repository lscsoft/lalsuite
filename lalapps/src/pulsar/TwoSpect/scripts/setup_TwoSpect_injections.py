#!/usr/bin/python

import argparse, sys, os

parser = argparse.ArgumentParser(description='Script that sets up the condor/dag file for TwoSpect Monte Carlo injections')
parser.add_argument('--runScript', required=True, type=str, help='Path/file of run script')
parser.add_argument('--twospectExe', required=True, type=str, help='Path to TwoSpect executable')
parser.add_argument('--dir', required=True, type=str, help='Base directory where the job subdirectories are located')
parser.add_argument('--jobs', required=True, type=int, help='Number of jobs, also the base directory/subdirectory path')
parser.add_argument('--logfile', required=True, type=str, help='Path/filename of logfile')
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

os.mkdir(args.dir)
os.mkdir('{}/err'.format(args.dir))

dagfile = open('{}/dag'.format(args.dir),'w')
for ii in range(0, args.jobs):
    dagfile.writelines(['JOB A{} {}/condor\n'.format(ii, args.dir),'VARS A{} JOBNUM=\"{}\"\n'.format(ii, ii)])
    os.mkdir('{}/{}'.format(args.dir, ii))
dagfile.close()

condorfile = open('{}/condor'.format(args.dir),'w')
condorfile.write("""\
universe=vanilla
executable={}
input=/dev/null
output=/dev/null
error={}/err/err.$(JOBNUM)
log={}
request_memory=2500
notification=Never
""".format(args.runScript, args.dir, args.logfile))

condorfile.write('arguments=\"--dir={} --twospectExe={} --jobnum=$(JOBNUM) --seed={} --Tobs={} --Tsft={} --SFToverlap={} --fmin={} --fspan={} --injPol={} --t0={}'.format(args.dir, args.twospectExe, args.seed, args.Tobs, args.Tsft, args.SFToverlap, args.fmin, args.fspan, args.injPol, args.t0))

if args.h0val: condorfile.write(' --h0val={}'.format(args.h0val))
else: condorfile.write(' --h0min={} --h0max={} --h0dist={}'.format(args.h0min, args.h0max, args.h0dist))

if args.scox1: condorfile.write(' --scox1')
else:
    condorfile.write(' --injPmin={} --injPmax={} --periodDist={} --injDfmin={} --injDfmax={}'.format(args.injPmin, args.injPmax, args.periodDist, args.injDfmin, args.injDfmax))
    if args.injskyra and args.injskydec: condorfile.write(' --injskyra={} --injskydec={}'.format(args.injskyra, args.injskydec))
    if args.maxEcc>0.0: condorfile.write(' --minEcc={} --maxEcc={} --eccDist={}'.format(args.minEcc, args.maxEcc, args.eccDist))
    if args.minSpindown!=0.0 or args.maxSpindown!=0.0: condorfile.write(' --minSpindown={} --maxSpindown={} --spindownDist={}'.format(args.minSpindown, args.maxSpindown, args.spindownDist))

condorfile.write(' --IFO=')
for ii in range(0, numIFOs):
    condorfile.write(IFOList[ii])
    if ii<numIFOs-1: condorfile.write(',')

if args.timestampsFile:
    condorfile.write(' --timestampsFile=')
    for ii in range(0, numIFOs):
        condorfile.write(timestampFileList[ii])
        if ii<numIFOs-1: condorfile.write(',')
elif args.segmentFile:
    condorfile.write(' --segmentFile=')
    for ii in range(0, numIFOs):
        condorfile.write(segmentFileList[ii])
        if ii<numIFOs-1: condorfile.write(',')

if args.inputSFTs: condorfile.write(' --inputSFTs={}'.format(args.inputSFTs))
if args.gaussianNoiseWithSFTgaps: condorfile.write(' --gaussianNoiseWithSFTgaps')
if args.skylocations!=0: condorfile.write(' --skylocations={}'.format(args.skylocations))
if args.ihsfactor!=5: condorfile.write(' --ihsfactor={}'.format(args.ihsfactor))
if args.markBadSFTs: condorfile.write(' --markBadSFTs')
if args.ulonly: condorfile.write(' --ulonly')
if args.templateTest: condorfile.write(' --templateTest')
if args.ihsfar!=1.0: condorfile.write(' --ihsfar={}'.format(args.ihsfar))
if args.ihsfomfar!=1.0: condorfile.write(' --ihsfomfar={}'.format(args.ihsfomfar))
if args.tmplfar!=1.0: condorfile.write(' --tmplfar={}'.format(args.tmplfar))
if args.tmplLength!=500: condorfile.write(' --tmplLength={}'.format(args.tmplLength))
if args.injDfExpAllow!=1.0: condorfile.write(' --injDfExpAllow={}'.format(args.injDfExpAllow))
if args.keepOnlyTopNumIHS>-1: condorfile.write(' --keepOnlyTopNumIHS={}'.format(args.keepOnlyTopNumIHS))
if args.useCorrectNScosi: condorfile.write(' --useCorrectNScosi')
if args.useCorrectNSpsi: condorfile.write(' --useCorrectNSpsi')

condorfile.write('\"\n')
condorfile.write('queue')
condorfile.close()
