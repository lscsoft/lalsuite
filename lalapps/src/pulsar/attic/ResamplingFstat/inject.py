#! /usr/bin/python
# inject.py , injects signals into a bunch of SFTs
import sys , commands , random
Tsft = 10000.0
fmin = 100.0
Band = 12.0
h0 = 1
cosi = 0
psi = 0
phi0 = 0
F = 103.0
A = 0.5
D = 0.5
IFO = 'LHO'
t0 = 700000000
dur = 10000
refTime = 700000000
Out= './SFTs'
NSFT = './SFTs/*.sft'
Ephem = '${LAL_LOCATION}/share/lal'
Sigma = 0.0
outfile = 'T'
LOG = open('logfile','w')

def gen(first):
	if(first):
	#	S = 'lalapps_Makefakedata ' + ' --Tsft ' + str(Tsft) + ' --fmin ' + str(fmin)+ ' --Band ' + str(Band) + ' --h0 ' + str(h0) + ' --IFO ' + IFO + ' --startTime ' + str(t0) + ' --duration ' + str(dur) + ' --refTime ' + str(refTime) + ' --outSFTbname ' + Out + ' -E ' + Ephem + ' -t' + ' TS' + ' --lineFeature' + ' --Freq ' + str(F)
	        S = '~/opt/lscsoft/lalapps/bin/lalapps_Makefakedata ' + ' --Tsft ' + str(Tsft) + ' --fmin ' + str(fmin)+ ' --Band ' + str(Band) + ' --h0 ' + str(h0) + ' --cosi ' + str(cosi) + ' --psi ' + str(psi) + ' --phi0 ' + str(phi0) + ' --Freq ' + str(F) + ' --Alpha ' +str(A) + ' --Delta ' + str(D) + ' --IFO ' + IFO + ' --startTime ' + str(t0) + ' --duration ' + str(dur) + ' --refTime ' + str(refTime) + ' --outSFTbname ' + Out + ' -E ' + Ephem + ' --noiseSigma ' + str(Sigma) 
		return(S)
	else:
		 S = '~/opt/lscsoft/lalapps/bin/lalapps_Makefakedata ' + ' --Tsft ' + str(Tsft) + ' --fmin ' + str(fmin)+ ' --Band ' + str(Band) + ' --h0 ' + str(h0) + ' --cosi ' + str(cosi) + ' --psi ' + str(psi) + ' --phi0 ' + str(phi0) + ' --Freq ' + str(F) + ' --Alpha ' +str(A) + ' --Delta ' + str(D) + ' --IFO ' + IFO + ' --startTime ' + str(t0) + ' --duration ' + str(dur) + ' --refTime ' + str(refTime) + ' --outSFTbname ' + Out + ' -E ' + Ephem + ' --noiseSFT ' + NSFT
		 return(S)
	 	#S = 'lalapps_Makefakedata ' + ' --Tsft ' + str(Tsft) + ' --fmin ' + str(fmin)+ ' --Band ' + str(Band) + ' --h0 ' + str(h0) + ' --IFO ' + IFO + ' --startTime ' + str(t0) + ' --duration ' + str(dur) + ' --refTime ' + str(refTime) + ' --outSFTbname ' + Out + ' -E ' + Ephem + ' -t' + ' TS' + ' --lineFeature' + ' --Freq ' + str(F) + ' --noiseSFT ' + NSFT


RMOLD = 'rm ' + Out + '/*'
commands.getoutput(RMOLD)

first = 1
for i in range(1,2):
	F = 103.0 + 0.1*(i-1) #i * 0.1 + 4.0 + fmin
	A = 0
	D = 0 #random.uniform(0,0.3)
	#F = random.uniform(fmin+Band/4.0,fmin+Band*3.0/4.0)
	logline = 'Injection (' + str(i) + ') ' ' F = ' + str(F) + ' Alpha = ' + str(A) + ' Delta = ' + str(D) + '\n'
	LOG.write(logline)
	if(first):
		Band = 20.0
	
	S = gen(first)
	if(first):
		first = 0
		Band = 12.0
		
	K = commands.getoutput(S);
	print S,'\n',K,'\n',F,'\n'

#GenT = './proc ' + str(t0) + ' ' + str(t0+dur) + ' ' + 'H1 ' + str(1.0) + ' ' + str(11) + ' ./SFTs/*.sft > ' + outfile 
#print GenT 
#K =  commands.getoutput(GenT)
#commands.getoutput(' mv T ../resamp/')
#print K

	

