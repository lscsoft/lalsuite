#! /usr/bin/python
# inject.py , injects signals into a bunch of SFTs
import sys , commands , random , math
#from pylab import *
Tsft = 100000*10  
fmin = 104
Band = 5.0
h0 = 0.1
cosi = 0
psi = 0
phi0 = 0
F = 105
A = 0
D = 0
IFO = 'LHO'
t0 = 700010000
dur = 100000*10 
refTime = 700000000
Out= './SFTs' + str(round(random.uniform(0,9),2))
NSFT = './SFTs/*.sft'
Ephem = '/opt/lscsoft/lal/share/lal'
EphemYear = '00-04'
Sigma = 0.0
outfile = 'T' + str(round(random.uniform(0,9),2))
LOG = open('\archive\home\ppatel\CONDOR\resamp\runs\logfile','w')

LOG.write('\n' + commands.getoutput('pwd') + '\n')

def gen(first):
	if(first):
		S = '~/opt/lscsoft/lalapps/bin/lalapps_Makefakedata ' + ' --Tsft ' + str(Tsft) + ' --fmin ' + str(fmin)+ ' --Band ' + str(Band) + ' --h0 ' + str(h0) + ' --cosi ' + str(cosi) + ' --psi ' + str(psi) + ' --phi0 ' + str(phi0) + ' --Freq ' + str(F) + ' --Alpha ' +str(A) + ' --Delta ' + str(D) + ' --IFO ' + IFO + ' --startTime ' + str(t0) + ' --duration ' + str(dur) + ' --refTime ' + str(refTime) + ' --outSFTbname ' + Out + ' -E ' + Ephem + ' --noiseSigma ' + str(Sigma) + ' -t' + ' TS'
		return(S)
	else:
		S = '~/opt/lscsoft/lalapps/bin/lalapps_Makefakedata ' + ' --Tsft ' + str(Tsft) + ' --fmin ' + str(fmin)+ ' --Band ' + str(Band) + ' --h0 ' + str(h0) + ' --cosi ' + str(cosi) + ' --psi ' + str(psi) + ' --phi0 ' + str(phi0) + ' --Freq ' + str(F) + ' --Alpha ' +str(A) + ' --Delta ' + str(D) + ' --IFO ' + IFO + ' --startTime ' + str(t0) + ' --duration ' + str(dur) + ' --refTime ' + str(refTime) + ' --outSFTbname ' + Out + ' -E ' + Ephem + ' --noiseSFT ' + NSFT
		return(S)

def semi():
	S = '~/opt/lscsoft/lalapps/bin/lalapps_SemiAnalyticF ' + ' --h0 ' + str(h0) + ' --cosi ' + str(cosi) + ' --psi ' + str(psi) + ' --phi0 ' + str(phi0) + ' --Alpha ' +str(A) + ' --Delta ' + str(D) + ' --IFO ' + IFO + ' --duration ' + str(dur) + ' --ephemDir ' + Ephem + ' --ephemYear ' + EphemYear + ' -S ' + str(t0) + ' -N ' + str(1.0**0.5)
	return(S)

def ignore(IGNORED,STRNG):
	#DEBUG = open('debug','w')
	#DEBUG.write(IGNORED+'\n')
	TEMP = STRNG.split('\n')
	#DEBUG.write(STRNG + '\n')
	#print TEMP
	NEWSTR = ''
	for i in TEMP:
		if((i.split())[0] != IGNORED):
			#DEBUG.write((i.split())[0] + '\n')
			NEWSTR += (i + '\n')

	return(NEWSTR)

CHDIR = 'cd /usr1/ppatel/'
commands.getoutput(CHDIR)
MKDIR = 'mkdir -p temp'
commands.getoutput(MKDIR)
commands.getoutput('cd temp')
#MKDIR = 'mkdir -p runs'
#commands.getoutput(MKDIR)
#MKDIR = 'mkdir -p run_' + '${process}'
#print MKDIR
#commands.getoutput(MKDIR)
#CHDIR = 'cd resamp/runs/run_' +  '${process}'
MKDIR = 'mkdir -p ' + Out
commands.getoutput(MKDIR)
RMOLD = 'rm ' + Out + '/*'
commands.getoutput(RMOLD)

first = 1
Freqs = []
Fcalc1 = []
Fcalc2 = []
Diff1 = []
Diff2 = []
for i in range(1,31):
	#F = 7.22
	
	A = round(random.uniform(0,6.28),2)
	D = round(math.asin(2.0*random.uniform(0,1.0)-1.0),2)*0
	F = 106.0 #random.uniform(5.5,6.5) + 100
	#F = round(F,2)
	logline = 'Injection (' + str(i) + ') ' ' F = ' + str(F) + ' Alpha = ' + str(A) + ' Delta = ' + str(D) 
	print >> sys.stderr,('\n' + logline)
	S = gen(first)
	print >> sys.stderr,('\n' + S)
	K = commands.getoutput(S)
	print >> sys.stderr,('\n' + K)
	GenT = '~/CONDOR/resamp/proc ' + str(t0) + ' ' + str(t0+dur) + ' ' + 'H1 ' + str(5.0+100) + ' ' + str(7.0+100) + ' ./' + Out + '/*.sft ' + str(Tsft) + ' > ' + outfile
	print >> sys.stderr,('\n' + GenT)
	K = commands.getoutput(GenT)
	print >> sys.stderr,('\n' + K)
	Exec = '~/CONDOR/resamp/fstat ' + ' -a ' + str(A) + ' -d ' + str(D) + ' -f ' + str(F) + ' -i ' + outfile
	print >> sys.stderr,('\n' + Exec)
	dur = 86400*10-500
	print >> sys.stderr,('\n' + semi())
	FS = commands.getoutput(semi())
	FS = ignore('Condor:',FS)
	print >> sys.stderr,('\n' + FS)
	dur = 100000*10
	FScalc = commands.getoutput(Exec)
	print >> sys.stderr,('\n' + FScalc)
	Freqs.append(F)
	Split = FScalc.split()
	Ftemp = float(Split[3])/2.0
	Fcalc1.append(Ftemp)
	print A,D,F,Ftemp,FS

commands.getoutput('cd ..')
commands.getoutput('rm -rf temp')
	

