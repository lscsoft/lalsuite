#!/usr/bin/python
# Resampling Test Code. This code is much more readable and reusable than its predecessor. Author - Pinkesh Patel
import sys,random,commands,os,string,math
from numpy.fft import fft
from numpy.fft import ifft

def main():
    
    # Import the Configuration File into config
    if(len(sys.argv) < 2):
        print "Insufficient command line arguments "
        print " Usage is ",sys.argv[0]," <ConfigFile> {Job Number} {Cluster Number} "
        print " Exiting ........."
        sys.exit(1)
    
    try:
        config = __import__(sys.argv[1]);
    except:
        print "Cannot import Config file '",sys.argv[1],"' exiting ...."
        sys.exit(1)

    UniqueID = ''

    # Compute Unique ID Number
    if(len(sys.argv) >= 3):
        UniqueID += '_' + sys.argv[2]
        if(len(sys.argv) == 4):
               UniqueID += '_'+sys.argv[3]
   
    # Variables Structure/Dictionary
    Vars = {}

    # SFT Output Directory
    Vars['Out'] = '__Temp_SFT'+ UniqueID

    # SFT time baseline
    try:
        Vars['TSFT'] = config.TSFT
    except:
        print "TSFT cannot be read"
        sys.exit(1)

    # Strength of signal
    try:
        Vars['h0'] = config.h0
    except:
        print "h0 cannot be read"
        sys.exit(1)

    # Cosine of iota
    try:
        Vars['cosi'] = config.cosi
    except:
        try:
            Vars['cosi'] = random.uniform(config.cosi_min,config.cosi_max)
        except:
            print "Cannot read in cosi variable"
            sys.exit(1)
    if(Vars['cosi'] < -1):
        print "cosi out of bounds !!!"
        sys.exit(1)


    # Initial Phase
    try:
        Vars['phi0'] = config.phi0
    except:
        try:
            Vars['phi0'] = random.uniform(config.phi0_min,config.phi0_max)
        except:
            print "Cannot read in phi0 variable"
            sys.exit(0)
    if(Vars['phi0'] < 0):
        print "phi0 < 0 !!!"
        sys.exit(0)

    # Polarization Angle
    try:
        Vars['psi'] = config.psi
    except:
        try:
            Vars['psi'] = random.uniform(config.psi_min,config.psi_max)
        except:
            print "Cannot read in psi variable"
            sys.exit(0)
    if(Vars['psi'] < 0):
        print "psi < 0 !!!"
        sys.exit(0)

    # Number of Dirichlet terms used.
    try:
        Vars['Dterms'] = int(config.Dterms)
    except:
        print "Dterms cannot be read"
        sys.exit(0)

    # Interferometer
    try:
        Vars['IFO'] = config.IFO
    except:
        print "IFO cannot be read"
        sys.exit(0)

    # Start Time
    try:
        Vars['t0'] = float(config.t0)
    except:
        print "t0 cannot be read"
        sys.exit(0)

    # Reference Time in SSB
    try:
        Vars['refTime'] = float(config.refTime)
    except:
        print "refTime cannot be read"
        sys.exit(0)

    # Ephemeris Directory
    try:
        Vars['Ephem'] = config.Ephem
    except:
        print "Ephem cannot be read"
        sys.exit(0)

    # Ephemeris Year
    try:
        Vars['EphemYear'] = config.EphemYear
    except:
        print "EphemYear cannot be read"
        sys.exit(0)

    # Noise Sh
    try:
        Vars['Sh'] = float(config.Sh)
    except:
        print "Sh cannot be read"
        sys.exit(0)
 
    # Duration of Analysis
    try:
        Vars['TSpan'] = int(config.TSpan)
    except:
        print "TSpan cannot be read"
        sys.exit(0)

    # Number of SFTs to add
    try:
        Vars['NumSFTs'] = int(config.NumSFTs)
    except:
        print "NumSFTs cannot be read"
        sys.exit(0)

    # Number of Gaps to add
    try:
        Vars['NumGaps'] = int(config.NumGaps)
    except:
        print "NumGaps cannot be read"
        sys.exit(0)

    # Alpha (Right Ascension)
    try:
        Vars['Alpha'] = config.Alpha
    except:
        try:
            Vars['Alpha'] = random.uniform(config.Alpha_min,config.Alpha_max)
        except:
            print "Cannot read in Alpha variable"
            sys.exit(0)
    if(Vars['Alpha'] < 0):
        print "Alpha < 0 !!!"
        sys.exit(0)

    # Delta (Declination)
    try:
        Vars['Delta'] = config.Delta
    except:
        try:
            Vars['Delta'] = random.uniform(config.Delta_min,config.Delta_max)
        except:
            print "Cannot read in Delta variable"
            sys.exit(0)
    if(Vars['Delta'] < 0):
        print "Delta < 0 !!!"
        sys.exit(0)

    # Minimum Frequency
    try:
        Vars['Fmin'] = config.Fmin
    except:
        print "Fmin cannot be read"
        sys.exit(0)

    # Band of Analysis
    try:
        Vars['Band'] = config.Band
    except:
        print "Band cannot be read"
        sys.exit(0)

    # Injection Frequency
    try: 
        Vars['Finj'] = config.Finj
    except:
        try:
            Vars['Finj'] = random.uniform(config.Finj_min,config.Finj_max)
        except:
            print "Cannot read in Finj variable"
            sys.exit(0)
    if(Vars['Finj'] < 0):
        print "Finj < 0 !!!"
        sys.exit(0)

    # Spindown/ FDOT
    try: 
        Vars['FDot'] = config.FDot
    except:
        try:
            Vars['FDot'] = random.uniform(config.FDot_min,config.FDot_max)
        except:
            print "Cannot read in FDot variable"
            sys.exit(0)
    if(Vars['FDot'] < 0):
        print "FDot < 0 !!!"
        sys.exit(0)

    # Resolution
    try:
        Vars['Res'] = config.Res
        if(Vars['Res'] > 1.0/Vars['TSpan']):
            print "Resolution too low, set to 1/T"
            Vars['Res'] = 1.0/Vars['TSpan']
        if(Vars['Res'] < 0):
            print "Resolution < 0"
            sys.exit(0)
    except:
        Vars['Res'] = 1.0/Vars['TSpan']

    # Debug Check
    try:
        Vars['debug'] = config.debug
    except:
        Vars['debug'] = 1.0

    # F Threshold
    try:
        Vars['FThres'] = config.FThres
    except:
        Vars['FThres'] = 0

    # Optional OutputTimeSeries
    try:
        Vars['TimeSeriesOut'] = config.TimeSeriesOut
        TimeSeriesOut_Is_Set = True
    except:
        TimeSeriesOut_Is_Set = False

    # Optional Use Your Own TimeStampsFile
    try:
        Vars['TimeStampsFile'] = config.TimeStampsFile
        TimeStampsFile_Is_Set = True
    except:
        Vars['TimeStampsFile'] = 'TimeStampsFile' + UniqueID
        TimeStampsFile_Is_Set = False

    # Optional Output Variable from Resamp
    try:
        Vars['ResampOutput'] = config.ResampOutput
    except:
        Vars['ResampOutput'] = "MyTS"

    # Print out all the variables (if Debug is on)
    if(Vars['debug']):
        print "---------- Configuration Variables --------- "
        PrintValues(Vars)
        print "-------------------------------------------- \n \n"

    # Create the Time Stamp File if none was specified
    if(not(TimeStampsFile_Is_Set)):
        CreateTimeStampFile(Vars)

    # If running multiple time, delete all the old SFTs
    if(os.path.isdir(Vars['Out'])):
        RMOLD = 'rm -rf ' + Vars['Out'] + '/*'
        print RMOLD
        (status,output) = commands.getstatusoutput(RMOLD)
        if(status):
            print "Failed to delete old SFTs from ",Vars['Out']," Exiting ...."
            sys.exit(1)
    else:
        try:
            os.mkdir(Vars['Out'])
        except:
            print " Failed to create ",Vars['Out']," Exiting ....."
            sys.exit(1)

    # Generate Fake data string
    # If TimeSeries Output is set, then add it in 
    if(TimeSeriesOut_Is_Set):
        optionalstring = ' --TDDfile ' + Vars['TimeSeriesOut'] + ' '
        FakeDataString = GenFakeDataString(1,Vars,optionalstring)
    else:
        FakeDataString = GenFakeDataString(1,Vars)
    
    if(Vars['debug']):
        print "----------- Makefakedata String ------------"
        print FakeDataString
        print "--------------------------------------------\n\n"

    # Generate the data
    (status,output) = commands.getstatusoutput(FakeDataString)
    if(status):
        print "Tried to generate SFTs, Command - ",FakeDataString," failed"
        print "Output was ",output
        sys.exit(1)
    
    # Run Resamp
    OutputFile = "OutputR"
    startstring = "./lalapps_ComputeFStatistic_resamp --outputFstat " + OutputFile + " -F " + str(Vars['FThres']) + " "
    endstring = " -S " + " -t " + str(Vars['Dterms']) + " --outputTimeSeries " + str(Vars['ResampOutput'])
    RDataString = GenDataString(startstring,endstring,Vars)
    if(Vars['debug']):
        print "-------- Resamp String -----------"
        print RDataString
        print "----------------------------------\n\n"
        print "--------- Running Resamp ---------"
    
    (status,output) = commands.getstatusoutput(RDataString)
    if(status):
        print "Tried to run Resamp, Command - ",RDataString," failed"
        print "Output was ",output
        sys.exit(1)
    
    print output
        
    if(Vars['debug']):
        print "---------- Resamp Done -----------\n\n"

    if(Vars['debug']):
        print "---------- Reading in Resamp Time Series ------- \n\n"

    try:
        File = open(Vars['ResampOutput'],'r')
    except:
        print "Could not open ",Vars['ResampOutput'],"\n\n"
        exit(1)

    TimeR = []
    RealR = []
    ImagR = []
    Het_Freq = 0
    LinebyLine = string.split(File.read(),'\n')

    for line in LinebyLine:
        SplitLine = line.split()
        
        if(len(SplitLine) > 0):
            if(SplitLine[0] == '$'):
                Het_Freq = float(SplitLine[1])
        
        if(len(SplitLine) == 5):
            TimeR.append(float(SplitLine[0]))
            RealR.append(float(SplitLine[1]))
            ImagR.append(float(SplitLine[2]))
    
    
    if(Vars['debug']):
        print "---------- Done reading Resamp Time Series ------- \n\n"
        print "---------- Heterodyne Freq is ",Het_Freq,"-------------\n\n"
    
    if(Vars['debug']):
        print "---------- Reading in MFD Time Series ------- \n\n"

    try:
        File = open(Vars['TimeSeriesOut']+'.00','r')
    except:
        print "Could not open ",Vars['TimeSeriesOut']+ '.00',"\n\n"
        exit(1)

    Time = []
    HetData  = []
    Data = []
    
    LinebyLine = string.split(File.read(),'\n')

    for line in LinebyLine:
        SplitLine = line.split()
        if(len(SplitLine) == 2):
            Time.append(float(SplitLine[0]))
            Data.append(float(SplitLine[1]))
    
    
    if(Vars['debug']):
        print "---------- Done reading MFD Time Series ------- \n\n"
        print "---------- Heterodyning now ----------\n\n"

    StartTime = Time[0]
    for i in range(len(Data)):
        Time[i] = Time[i] - StartTime
        SinPhi = math.sin(2.0*math.pi*Het_Freq*Time[i])
        CosPhi = math.cos(2.0*math.pi*Het_Freq*Time[i])
        Real = Data[i]*CosPhi
        Imag = Data[i]*(SinPhi)*complex(0,1)
        HetData.append(Real + Imag)

    if(Vars['debug']):
        print "---------- Filtering ------- \n\n"
    
    fftdata = fft(HetData)
    fftdata[len(fftdata)/2:len(fftdata)] = 0
    FiltData = ifft(fftdata)

    if(Vars['debug']):
        print "---------- Writing Out the Time Series ------- \n\n"

    OutFile = open("FiltTS",'w-')
    for i in range(len(FiltData)):
        OutFile.write(str(Time[i]) + ' ' + str(FiltData[i].real) + ' ' + str(FiltData[i].imag)+ '\n')

    return(0)

def PrintValues(Dict):
    for key in Dict.keys():
        print key," = ",Dict[key]

def CreateTimeStampFile(Vars):
    try:
        File = open('./'+Vars['TimeStampsFile'],'w')
    except:
        print "Tried to open timestampsFile, failed"
        sys.exit(0)
    
    if(Vars['debug']):
        print "----------- Starting Random Time Stamp Creation -------------"

    t0 = Vars['t0']
    NumSFTs = Vars['NumSFTs']
    NumGaps = Vars['NumGaps']
    TSpan = Vars['TSpan']
    TSFT = Vars['TSFT']
    GapTime = TSpan - TSFT*NumSFTs
    Gaps = []
    Chunks = []

    GapTimeleft = GapTime
    SFTsleft = NumSFTs
    for i in range(NumGaps-1):
        Gap = round(random.uniform(0,GapTimeleft),0)
        GapTimeleft = GapTimeleft - Gap
        Gaps.append(Gap)

    Gaps.append(GapTimeleft) 

    for i in range(NumGaps):
        Chunk = int(round(random.uniform(0,SFTsleft),0))
        SFTsleft = SFTsleft - Chunk
        Chunks.append(Chunk)

    Chunks.append(SFTsleft)

    TimeNow = t0
    for i in range(NumGaps):
        if(Vars['debug']):
            print "Chunk number ",i+1," has ",Chunks[i]," SFTs "
            print Gaps[i]," is the gap after this chunk "

        for j in range(Chunks[i]):
            File.write(str(int(TimeNow + j*TSFT)))
            File.write(" 0 \n")
        TimeNow = TimeNow + Gaps[i] + TSFT*Chunks[i]

    if(Vars['debug']):
        print "Chunk number ",NumGaps+1," has ",Chunks[NumGaps]," SFTs "
        print "--------------------------------------------\n \n "

    for j in range(Chunks[NumGaps]):
        File.write(str(int(TimeNow + j*TSFT)))
        File.write(" 0 \n")
            

def GenFakeDataString(addtonoise,Vars,optionalstring = ' '):
    #CreationBand = Vars['Band']*2
    #CreationFmin = Vars['Fmin']-Vars['Band']/2
    CreationBand = Vars['Fmin'] + Vars['Band']*2
    CreationFmin = 0
    S = 'lalapps_Makefakedata_v4 ' + ' --Tsft ' + str(Vars['TSFT']) + ' --fmin ' + str(CreationFmin) + ' --h0 ' + str(Vars['h0']) + ' --Band ' + str(CreationBand) + ' --cosi ' + str(Vars['cosi']) + ' --psi ' + str(Vars['psi']) + ' --phi0 ' + str(Vars['phi0']) + ' --Freq ' + str(Vars['Finj']) + ' --Alpha ' + str(Vars['Alpha']) + ' --Delta ' + str(Vars['Delta']) + ' --IFO ' + str(Vars['IFO']) + ' --refTime ' + str(Vars['refTime']) + ' --outSFTbname ' + str(Vars['Out']) + ' --ephemDir ' + str(Vars['Ephem']) + ' --ephemYear ' + str(Vars['EphemYear']) + ' --f1dot ' + str(Vars['FDot']) + ' --noiseSqrtSh ' + str(Vars['Sh']**0.5) + ' --timestampsFile ' + Vars['TimeStampsFile'] + ' --generationMode 0 ' + optionalstring
    return(S)

                         
def GenDataString(beginstring,endstring,Vars):
    S = beginstring + ' --Freq ' + str(Vars['Fmin']) + ' --FreqBand ' + str(Vars['Band']) + ' --Alpha ' + str(Vars['Alpha']) + ' --Delta ' + str(Vars['Delta']) + ' --IFO ' + str(Vars['IFO']) + ' --refTime ' + str(Vars['refTime']) + ' --ephemDir ' + str(Vars['Ephem']) + ' --ephemYear ' + str(Vars['EphemYear']) + ' --f1dot ' + str(Vars['FDot']) + ' --dFreq ' + str(Vars['Res']) + ' --DataFiles \"' + str(Vars['Out']) + '/*" ' + endstring
    return(S)
 
if __name__ == '__main__':
    main()


















