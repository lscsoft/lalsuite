#! /usr/bin/python
# This program runs makefakedata, analyzes it and spits out the required data.
# Author - Pinkesh Patel
import sys , commands , random, math, os


# Generates the command line arguments for Makefakedata, if it is not in your path, please replace the commands referring to it.
# The argument is if you are injecting more that one signal into the SFTs. Not really used here.
 
def gen(first):
	if(first):
	        S = 'lalapps_Makefakedata_v4 ' + ' --Tsft ' + str(Tsft) + ' --fmin ' + str(fmin)+ ' --Band ' + str(Band) + ' --h0 ' + str(h0) + ' --cosi ' + str(cosi) + ' --psi ' + str(psi) + ' --phi0 ' + str(phi0) + ' --Freq ' + str(F) + ' --Alpha ' +str(A) + ' --Delta ' + str(D) + ' --IFO ' + IFO + ' --startTime ' + str(t0) + ' --duration ' + str(dur) + ' --refTime ' + str(refTime) + ' --outSFTbname ' + Out + ' -E ' + Ephem + ' --f1dot ' + str(FDOT) + ' --noiseSqrtSh=' + str(Sh) + ' --ephemYear=' + str(EphemYear) #+ ' --lineFeature ' 
		return(S)
	else:
		 S = 'lalapps_Makefakedata_v4 ' + ' --Tsft ' + str(Tsft) + ' --fmin ' + str(fmin)+ ' --Band ' + str(Band) + ' --h0 ' + str(h0) + ' --cosi ' + str(cosi) + ' --psi ' + str(psi) + ' --phi0 ' + str(phi0) + ' --Freq ' + str(F) + ' --Alpha ' +str(A) + ' --Delta ' + str(D) + ' --IFO ' + IFO + ' --startTime ' + str(t0) + ' --duration ' + str(dur) + ' --refTime ' + str(refTime) + ' --outSFTbname ' + Out + ' -E ' + Ephem + ' --noiseSFT ' + NSFT + ' --f1dot ' + str(FDOT)
		 return(S)


# Generates commands for resamp 	 	
def resamp():
	S = 'time ./lalapps_ComputeFstatistic_resamp ' +  ' -b ' + str(Band) + ' -f ' + str(fmin) + ' --Alpha ' +str(A) + ' --Delta ' + str(D) + ' -E ' + Ephem + ' -t '+ str(Dterms) + ' -Z ' + str(F) + ' --f1dot ' + str(FDOT) + '  ' +  ' --outputFstat '+ str(OUTR) + ' --refTime ' + str(refTime) +  ' -S -F 0 ' + '  ' + ' --ephemYear=' + str(EphemYear) + ' -r ' + str(1.0/dur/10) + ' -D \"' + Out + '/*\" > plot1'
	return(S)

# Generates commands for ComputeFstat_v2	
def v2():
	S = 'lalapps_ComputeFstatistic_v2 ' +  ' -b ' + str(Band) + ' -f ' + str(fmin) + ' --Alpha ' +str(A) + ' --Delta ' + str(D) + ' -E ' + Ephem + ' -t '+ str(Dterms) + ' -D \"' + Out + '/*\" '  + ' --outputFstat ' + OUTV + ' -F 0 ' + ' -r ' + str(1.0/dur/10) + ' --f1dot ' + str(FDOT) + ' -S ' + ' --refTime ' + str(refTime) + ' --ephemYear=' + str(EphemYear) 
	return(S)

# Generates commands for SemiAnalyticF
def semi():
	S = 'lalapps_SemiAnalyticF ' + ' --h0 ' + str(h0) + ' --cosi ' + str(cosi) + ' --psi ' + str(psi) + ' --phi0 ' + str(phi0) + ' --Alpha ' +str(A) + ' --Delta ' + str(D) + ' --IFO ' + IFO + ' --duration ' + str(dur) + ' --ephemDir ' + Ephem + ' --ephemYear ' + EphemYear + ' -S ' + str(t0) + ' -N ' + str(Sh**0.5)
	return(S)

# mkdir subroutine
def mymkdir(newdir):
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        if tail:
            os.mkdir(newdir)

# SFT time baseline
Tsft = 1800
# Strength of signal
h0 = 0.1
# cosine of iota
cosi = 1
psi = 0.5
phi0 = 0
Dterms = 16
IFO = 'H2'
# Start Time
t0 = 820000000
# Duration of Analysis
dur = 1800*5
# Reference Time in SSB
refTime = 820000000
# Output Directory
Out= './SFTs'
# Noise Output, only needed if we are injecting multiple signals, NOT USED
NSFT = './SFTs/*.sft'
# Ephemeris Directory (Change as needed)
Ephem = '/Users/ppatel/home/opt/lscsoft/lal/share/lal'
# Ephemeris Year
EphemYear = '05-09'
# Noise Sigma
Sigma = 1*0
# Noise Sh
Sh = 1*0

# Make the OUT directory
mymkdir(Out)

# If running multiple time, delete all the old SFTs
RMOLD = 'rm ' + Out + '/*'
commands.getoutput(RMOLD)

# first is true, again not needed here
first = 1

# Loop , change the parameters you need
for j in range(1,2):
	# Spacing between chunks
	Spacing = 1307

	# Alpha
	A = 1.4596743833
	# Delta
	D =  3.141592/2

	# If you want to set it to random sky location 
	#A = round(random.uniform(0,6.28),2)
	#D = round(math.asin(2.0*random.uniform(0,1.0)-1.0),2)

	# F is the frequency of Signal
	# options to set random (edit as needed)
	#F = round(random.uniform(0,1.0),2)*4.0

	# Set FDOT.
	FDOT = -1e-8*0
	
	# Set Band for analysis
	Band = 1.0

	# Set F
	F =  59.54 + Band/2.0
	
	# Band is first used to create the SFTs, here the Band is large
	Band = Band*2.0
	
	# Minimum Frequency in the SFTs
	fmin = 59.0 

	# Generate Command for making SFTs
	Make = gen(first)

	# Print it if you want to see
	#print Make
	
	# Execute command to make SFTs
	for g in range(1,4):
		Temp = commands.getoutput(Make)
		t0 = t0 + dur + Spacing
		Make = gen(first)

	# Band for analysis.
	Band = Band/2.0

	# Minimum Frequency of analysis.
	fmin = 59.54

	# Output files for R-resampling and V-v2()
	OUTR = 'outputR'
	OUTV = 'outputV'

	# resamp() generates code for Resamp
	#print resamp()
	R = commands.getoutput(resamp())
	#print R

	# v2() generates code for v2()
	#print v2()
	V2 = commands.getoutput(v2())
	#print V2

	#print semi()
	Semi = commands.getoutput(semi())
	#print Semi

	#Plot the results
	if(1):
		P1 = commands.getoutput('./makeplot')
		P2 = commands.getoutput('./makeplotzoom')
		P3 = commands.getoutput('./makeplotzoom+')

	# Analyze the data, if statement only so that I can turn on or off
	# Essentially takes the twoF and Frequency values and spits them out.
	if(0):
		File1 = open('./outputV' , 'r')
		Contents = File1.read()
		Contents = Contents.split('\n')
		for i in Contents:
			X = i.split()
			if(len(X) > 0 and X[0] == 'Freq'):
				J = X[2]
				F1 = float(J[0:len(J)-1])
		
		       	if(len(X) > 0 and X[0] == 'twoF'):	
				J = X[2]
				D1 = float(J[0:len(J)-1])
				
		File2 = open('./outputR' , 'r')
	       	Contents = File2.read()
	       	Contents = Contents.split('\n')
	       	for i in Contents:
	       		X = i.split()
	       		if(len(X) > 0 and X[0] == 'Freq'):
	       			J = X[2]
	       			F2 = float(J[0:len(J)-1])
		
			if(len(X) > 0 and X[0] == 'twoF'):	
				J = X[2]
				D2 = float(J[0:len(J)-1])
		
		# Run semianalyticF()
		Semi = 2.0*float(commands.getoutput(semi()))
		# Print FDOT,Frequency,twoFvalue etc.
		print FDOT,F1,D1,F2,D2,FDOT*dur**2,(F2-F1)*dur,100.0*(D1/D2-1.0)
		# Print to stderr, in case I am looping and want to see progress
		print >> sys.stderr,FDOT,F1,D1,F2,D2,FDOT*dur**2,(F2-F1)*dur,100.0*(D1/D2-1.0)
		
	


