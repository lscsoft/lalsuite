#!/usr/bin/python
# Resampling Test Code. This code is much more readable and reusable than its predecessor. Author - Pinkesh Patel
import sys,random,commands,math

def main():
    
    # Import the Configuration File into config
    config = __import__(sys.argv[1]);
   
    # Variables Structure/Dictionary
    Vars = {}

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
    if(math.fabs(Vars['cosi']) > 1):
        print "abs(cosi) > 1 !!!"
        sys.exit(1)


    # Initial Phase
    try:
        Vars['phi0'] = config.phi0
    except:
        try:
            Vars['phi0'] = random.uniform(config.phi0_min,config.phi0_max)
        except:
            print "Cannot read in phi0 variable"
            sys.exit(1)

    # Polarization Angle
    try:
        Vars['psi'] = config.psi
    except:
        try:
            Vars['psi'] = random.uniform(config.psi_min,config.psi_max)
        except:
            print "Cannot read in psi variable"
            sys.exit(1)

    # Number of Dirichlet terms used.
    try:
        Vars['Dterms'] = float(config.Dterms)
    except:
        print "Dterms cannot be read"
        sys.exit(1)

    # Interferometer
    try:
        IFOs = config.IFOs
        IFOs = IFOs.split(' ')
        NumofIFOs = len(IFOs)
        Vars['NumofIFOs'] = NumofIFOs
        Vars['IFO'] = []
        for i in range(NumofIFOs):
            Vars['IFO'].append(IFOs[i])
    except:
        print "IFO cannot be read"
        sys.exit(1)

    # Start Time
    try:
        t0 = config.t0
        t0 = t0.split(' ')
        if(len(t0) != NumofIFOs):
            print "Number of starttimes != Number of IFOs"
            sys.exit(1)
        
        Vars['t0'] = []
        for i in range(NumofIFOs):
            Vars['t0'].append(float(t0[i]))
    except:
        print "t0 cannot be read"
        sys.exit(1)

    # Reference Time in SSB
    try:
        Vars['refTime'] = float(config.refTime)
    except:
        print "refTime cannot be read"
        sys.exit(1)

    # Output Directory
    try:
        Vars['Out'] = config.Out
    except:
        print "Out cannot be read"
        sys.exit(1)

    # Ephemeris Directory
    try:
        Vars['Ephem'] = config.Ephem
    except:
        print "Ephem cannot be read"
        sys.exit(1)

    # Ephemeris Year
    try:
        Vars['EphemYear'] = config.EphemYear
    except:
        print "EphemYear cannot be read"
        sys.exit(1)

    # Noise Sh
    try:
        Vars['Sh'] = float(config.Sh)
    except:
        print "Sh cannot be read"
        sys.exit(1)
 
    # Duration of Analysis
    try:
        TSpan = config.TSpan
        TSpan = TSpan.split(' ')
        if(len(TSpan) != NumofIFOs):
            print "Number of TSpans != Number of IFOs"
            sys.exit(1)
        
        Vars['TSpan'] = []
        for i in range(NumofIFOs):
            Vars['TSpan'].append(float(TSpan[i]))
    except:
        print "TSpan cannot be read"
        sys.exit(1)

    # Number of SFTs to add
    try:
        NumSFTs = config.NumSFTs
        NumSFTs = NumSFTs.split(' ')
        if(len(NumSFTs) != NumofIFOs):
            print "Number of starttimes != Number of IFOs"
            sys.exit(1)
        
        Vars['NumSFTs'] = []
        for i in range(NumofIFOs):
            Vars['NumSFTs'].append(int(NumSFTs[i]))
    except:
        print "NumSFTs cannot be read"
        sys.exit(1)

    # Number of Gaps to add
    try:
        NumGaps = config.NumGaps
        NumGaps = NumGaps.split(' ')
        if(len(NumGaps) != NumofIFOs):
            print "Number of starttimes != Number of IFOs"
            sys.exit(1)
        
        Vars['NumGaps'] = []
        for i in range(NumofIFOs):
            Vars['NumGaps'].append(int(NumGaps[i]))
    except:
        print "NumGaps cannot be read"
        sys.exit(1)

    # Alpha (Right Ascension)
    try:
        Vars['Alpha'] = config.Alpha
    except:
        try:
            Vars['Alpha'] = random.uniform(config.Alpha_min,config.Alpha_max)
        except:
            print "Cannot read in Alpha variable"
            sys.exit(1)
    if(Vars['Alpha'] < 0 or Vars['Alpha'] > 2.0*math.pi):
        print "Alpha out of bounds !!!"
        sys.exit(1)

    # Delta (Declination)
    try:
        Vars['Delta'] = config.Delta
    except:
        try:
            Vars['Delta'] = random.uniform(config.Delta_min,config.Delta_max)
        except:
            print "Cannot read in Delta variable"
            sys.exit(1)
    if(math.fabs(Vars['Delta']) > math.pi/2.0):
        print "abs(Delta) > pi/2 !!!"
        sys.exit(1)

    # Minimum Frequency
    try:
        Vars['Fmin'] = config.Fmin
    except:
        print "Fmin cannot be read"
        sys.exit(1)

    # Band of Analysis
    try:
        Vars['Band'] = config.Band
    except:
        print "Band cannot be read"
        sys.exit(1)

    # Injection Frequency
    try: 
        Vars['Finj'] = config.Finj
    except:
        try:
            Vars['Finj'] = random.uniform(config.Finj_min,config.Finj_max)
        except:
            print "Cannot read in Finj variable"
            sys.exit(1)
    if(Vars['Finj'] < 0):
        print "Finj < 0 !!!"
        sys.exit(1)

    # Spindown/ FDOT
    try: 
        Vars['FDot'] = config.FDot
    except:
        try:
            Vars['FDot'] = random.uniform(config.FDot_min,config.FDot_max)
        except:
            print "Cannot read in FDot variable"
            sys.exit(1)

    # FDOTBand
    try: 
        Vars['FDotBand'] = config.FDotBand
    except:
        Vars['FdotBand'] = 0
        #print "Cannot read in FDot variable"
        #sys.exit(1)

    # dFDot
    try: 
        Vars['dFDot'] = config.dFDot
    except:
        Vars['dFDot'] = 10

    # Resolution
    try:
        Vars['Res'] = config.Res
        if(Vars['Res'] > 1.0/Vars['TSpan'][0]):
            print "Resolution too low, set to 1/T"
            Vars['Res'] = 1.0/Vars['TSpan'][0]
        if(Vars['Res'] < 0):
            print "Resolution < 0"
            sys.exit(1)
    except:
        Vars['Res'] = 1.0/Vars['TSpan'][0]
        current = 0
        for i in range(NumofIFOs):
            if(Vars['TSpan'][i] > Vars['TSpan'][current]):
                current = i
                Vars['Res'] = 1.0/Vars['TSpan'][i]
            

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

    # Create the Time Stamp File
    CreateTimeStampFile(Vars)

    # If running multiple time, delete all the old SFTs
    RMOLD = 'rm ' + Vars['Out'] + '/*'
    commands.getoutput(RMOLD)

    # Generate Fake data string
    for ifo in range(NumofIFOs):
        FakeDataString = GenFakeDataString(1,Vars,ifo)
        if(Vars['debug']):
            print "----------- Makefakedata String ------------"
            print FakeDataString
            print "--------------------------------------------\n\n"

        # Generate the data
        try:
            G = commands.getoutput(FakeDataString)
        except:
            print "Tried to generate SFTs, failed"
            sys.exit(1)

    # Run v2()
    OutputFile = "OutputV"
    startstring = "lalapps_ComputeFStatistic_v2 --outputFstat " + OutputFile + " -F " + str(Vars['FThres']) + "  "
    endstring = " "
    V2DataString = GenDataString(startstring,endstring,Vars)
    if(Vars['debug']):
        print "----------- V2 String ------------"
        print V2DataString
        print "----------------------------------\n\n"
        print "----------- Running V2 -----------"
    
    try:
        G = commands.getoutput(V2DataString)
    except:
        print "V2 failed"
        sys.exit(0)
    
    if(Vars['debug']):
        print "---------- V2 Done ---------------\n\n"

    # Run Resamp
    OutputFile = "OutputR"
    startstring = "./lalapps_ComputeFStatistic_resamp --outputFstat " + OutputFile + " -F " + str(Vars['FThres']) + "  "
    endstring = "  --outputTimeSeries "  + str(Vars['ResampOutput']) + " > plot1 "
    RDataString = GenDataString(startstring,endstring,Vars)
    if(Vars['debug']):
        print "-------- Resamp String -----------"
        print RDataString
        print "----------------------------------\n\n"
        print "--------- Running Resamp ---------"
    
    try:
        G = commands.getoutput(RDataString)
    except:
        print "Resamp failed"
        sys.exit(0)
    
    if(Vars['debug']):
        print "---------- Resamp Done -----------\n\n"
        print "---------- Output was ------------\n\n"
        print G,"\n"
    
    return(0)

def PrintValues(Dict):
    for key in Dict.keys():
        print key," = ",Dict[key]

def CreateTimeStampFile(Vars):
    for ifo in range(Vars['NumofIFOs']):
        ifotimestampfile = "./timestampsFile" + str(Vars['IFO'][ifo])
        try:
            File = open(ifotimestampfile,'w')
        except:
            print "Tried to open timestampsFile, failed"
            sys.exit(0)
            
        if(Vars['debug']):
            print "----------- Starting Random Time Stamp Creation -------------"
                
        t0 = Vars['t0'][ifo]
        NumSFTs = int(Vars['NumSFTs'][ifo])
        NumGaps = int(Vars['NumGaps'][ifo])
        TSpan = Vars['TSpan'][ifo]
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
            

def GenFakeDataString(addtonoise,Vars,ifo):
    if(Vars['Band'] > 1e-2):
        CreationBand = Vars['Band']*4
        CreationFmin = Vars['Fmin']-Vars['Band']*2
    else:
        CreationBand = 1
        CreationFmin = Vars['Fmin'] - 0.5

    S = 'lalapps_Makefakedata_v4 ' + ' --Tsft ' + str(Vars['TSFT']) + ' --fmin ' + str(CreationFmin) + ' --h0 ' + str(Vars['h0']) + ' --Band ' + str(CreationBand) + ' --cosi ' + str(Vars['cosi']) + ' --psi ' + str(Vars['psi']) + ' --phi0 ' + str(Vars['phi0']) + ' --Freq ' + str(Vars['Finj']) + ' --Alpha ' + str(Vars['Alpha']) + ' --Delta ' + str(Vars['Delta']) + ' --IFO ' + str(Vars['IFO'][ifo]) + ' --refTime ' + str(Vars['refTime']) + ' --outSFTbname ' + str(Vars['Out']) + ' --ephemDir ' + str(Vars['Ephem']) + ' --ephemYear ' + str(Vars['EphemYear']) + ' --f1dot ' + str(Vars['FDot']) + ' --noiseSqrtSh ' + str(Vars['Sh']**0.5) + ' --timestampsFile timestampsFile' + str(Vars['IFO'][ifo])
    return(S)
                         
def GenDataString(beginstring,endstring,Vars):
    S = beginstring + ' --Freq ' + str(Vars['Fmin']) + ' --FreqBand ' + str(Vars['Band']) + ' --Alpha ' + str(Vars['Alpha']) + ' --Delta ' + str(Vars['Delta'])  + ' --refTime ' + str(Vars['refTime']) + ' --ephemDir ' + str(Vars['Ephem']) + ' --ephemYear ' + str(Vars['EphemYear']) +  ' --dFreq ' + str(Vars['Res']) + ' --DataFiles \"' + str(Vars['Out']) + '/*" ' + ' --f1dotBand ' + str(Vars['FDotBand']) + ' --df1dot ' + str(Vars['dFDot']) +  endstring
    return(S)

#' --f1dot ' + str(Vars['FDot']) + 
if __name__ == '__main__':
    main()


















