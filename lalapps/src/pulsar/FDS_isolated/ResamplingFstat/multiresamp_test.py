#!/usr/bin/python
# Resampling Test Code. This code is much more readable and reusable than its predecessor. Author - Pinkesh Patel
import sys,random,commands,math,shutil,os,re,string

def main():
    
    # Import the Configuration File into config
    # Import the Configuration File into config
    if(len(sys.argv) < 2):
        print "Insufficient command line arguments "
        print " Usage is ",sys.argv[0]," <ConfigFile> {Job Number} "
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
        refString = ' --refTime ' + str(Vars['refTime'])
    except:
        refString = ' '

    # Output Directory
    try:
        Vars['Out'] = config.Out + UniqueID
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
        Vars['FDotBand'] = 0
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

     # Optional Use Your Own TimeStampsFile
    try:
        Vars['TimeStampsFile'] = config.TimeStampsFile
        TimeStampsFile_Is_Set = True
    except:
        Vars['TimeStampsFile'] = 'TimeStampsFile' + UniqueID
        TimeStampsFile_Is_Set = False

    # Print out all the variables (if Debug is on)
    if(Vars['debug']):
        print "---------- Configuration Variables --------- "
        PrintValues(Vars)
        print "-------------------------------------------- \n \n"

    # Create the Time Stamp File
    if(not(TimeStampsFile_Is_Set)):
        CreateTimeStampFile(Vars,UniqueID)

    # If running multiple time, delete all the old SFTs
    if(os.path.exists("./"+Vars['Out'])):
        try:
            shutil.rmtree(Vars['Out'])
        except:
            print "Could not delete old directory\n"
            sys.exit(1)

    try:
        os.mkdir(Vars['Out'])
    except:
        print "Something went wrong creating new directory\n"
        print "Permission trouble maybe\n"
        sys.exit(1)

    # Generate Fake data string
    for ifo in range(NumofIFOs):
        FakeDataString = GenFakeDataString(1,Vars,ifo,UniqueID,refString)
        if(Vars['debug']):
            print "----------- Makefakedata String ------------"
            print FakeDataString
            print "--------------------------------------------\n\n"

        # Generate the data
        try:
            G = commands.getoutput(FakeDataString)
            os.remove("timestampsFile"+Vars['IFO'][ifo]+UniqueID)
        except:
            print "Tried to generate SFTs, failed"
            sys.exit(1)

    # Run v2()
    OutputFileVLoudest = "OutputVLoudest" + UniqueID
    OutputFileVFstat = "OutputVFstat" + UniqueID
    startstring = "lalapps_ComputeFStatistic_v2 --outputFstat " + OutputFileVFstat + " --outputLoudest " + OutputFileVLoudest +  " -F " + str(Vars['FThres']) + "  "
    endstring = " " + refString
    V2DataString = GenDataString(startstring,endstring,Vars)
    if(Vars['debug']):
        print "----------- V2 String ------------"
        print V2DataString
        print "----------------------------------\n\n"
        print "----------- Running V2 -----------"
    
    
    (status,Voutput) = commands.getstatusoutput(V2DataString)
    if(status):
        print "V2 failed, Output was \n\n"
        print Voutput
        sys.exit(1)
    
    if(Vars['debug']):
        print "---------- V2 Done ---------------\n\n"

    # Run Resamp
    ResampLocation = "/home/ppatel/lalsuite/lalapps/src/pulsar/FDS_isolated/ResamplingFstat/"
    OutputFileRLoudest = "OutputRLoudest" + UniqueID
    OutputFileRFstat = "OutputRFstat" + UniqueID
    startstring = ResampLocation + "lalapps_ComputeFStatistic_resamp --outputFstat " + OutputFileRFstat + " --outputLoudest " + OutputFileRLoudest+ " -F " + str(Vars['FThres']) + "  "
    endstring = " " + refString
    RDataString = GenDataString(startstring,endstring,Vars)
    if(Vars['debug']):
        print "-------- Resamp String -----------"
        print RDataString
        print "----------------------------------\n\n"
        print "--------- Running Resamp ---------"
    
    (status,Routput) = commands.getstatusoutput(RDataString)
    if(status):
        print "Resamp failed, Output was\n\n"
        print Routput
        sys.exit(0)
    
    if(Vars['debug']):
        print "---------- Resamp Done -----------\n\n"

    if(Vars['debug']):
        print "---------- Deleting SFT Folder ---------\n\n"
    try:
        shutil.rmtree(Vars['Out'])
    except:
        print " Could not delete SFT folder \n\n"
        sys.exit(1)
        
    FreqOutput = AnalyzeFreq(OutputFileRFstat,OutputFileVFstat,Vars['Finj'],Vars['Res'])
    LoudestOutput = AnalyzeLoudest(OutputFileRLoudest,OutputFileVLoudest)
    os.remove(OutputFileRLoudest)
    os.remove(OutputFileVLoudest)
    #os.remove(OutputFileRFstat)
    #os.remove(OutputFileVFstat)
    print FreqOutput,LoudestOutput
    return(0)

def AnalyzeFreq(Filename1,Filename2,Freq,dF):
    File1 = open(Filename1,'r')
    File2 = open(Filename2,'r')
    File1lines = File1.readlines()
    File2lines = File2.readlines()
    Freq1 = 0
    dF1 = 0
    Freq2 = 0
    dF2 = 0
    
    exp = re.compile(r'^\d')

    for line in File1lines:
        if(exp.search(line)):
            linesplit = string.split(line)
            if(Freq1 and not(dF1)):
                dF1 = abs(float(linesplit[0])) - Freq1
            Freq1 = float(linesplit[0])
            if(abs(Freq1-Freq) < dF1/2.0):
                twoF1 = float(linesplit[6])
                Freq1store = Freq1
                storeline = line
    
    for line in File2lines:
        if(exp.search(line)):
            linesplit = string.split(line)
            if(Freq2 and not(dF2)):
                dF2 = abs(float(linesplit[0])) - Freq2
            Freq2 = float(linesplit[0])
            if(abs(Freq2-Freq) < dF2/2.0):
                twoF2 = float(linesplit[6])
                Freq2store = Freq2
                

    #print Freq,Freq1store,Freq2store,dF,dF1,dF2,twoF1,twoF2
    return(str(Freq) + " " + str(Freq1store) + " " + str(Freq2store) + " " + str(dF) + " " + str(dF1) + " " + str(dF2) + " " + str(twoF1) + " " + str(twoF2) + " ")

def AnalyzeLoudest(Filename1,Filename2):
    File1 = open(Filename1,'r')
    File2 = open(Filename2,'r')
    File1lines = File1.readlines()
    File2lines = File2.readlines()
    
    expression = re.compile(r'(\D*)(\d*.\d*|\d*)')

    for line in File1lines:
        if(expression.search(line)):
            linesplit = expression.search(line).groups()
            if(re.match('twoF',linesplit[0])):
                twoF1 = float(linesplit[1])
            
            if(re.match('Freq',linesplit[0])):
                Freq1 = float(linesplit[1])
    
    for line in File2lines:
        if(expression.search(line)):
            linesplit = expression.search(line).groups()
            if(re.match('twoF',linesplit[0])):
                twoF2 = float(linesplit[1])
            
            if(re.match('Freq',linesplit[0])):
                Freq2 = float(linesplit[1])    

    #print Freq1,twoF1,Freq2,twoF2
    return(str(Freq1) + " " + str(twoF1) + " " + str(Freq2) + " " + str(twoF2) + " ")
            
def PrintValues(Dict):
    for key in Dict.keys():
        print key," = ",Dict[key]

def CreateTimeStampFile(Vars,UniqueID):
    for ifo in range(Vars['NumofIFOs']):
        ifotimestampfile = "./timestampsFile" + str(Vars['IFO'][ifo])+UniqueID
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
            

def GenFakeDataString(addtonoise,Vars,ifo,UniqueID,endstring):
    if(Vars['Band'] > 1e-2):
        CreationBand = Vars['Band']*4
        CreationFmin = Vars['Fmin']-Vars['Band']*2
    else:
        CreationBand = 1
        CreationFmin = Vars['Fmin'] - 0.5

    CreationFmin -= 1
    CreationBand += 2
    S = 'lalapps_Makefakedata_v4 ' + ' --Tsft ' + str(Vars['TSFT']) + ' --fmin ' + str(CreationFmin) + ' --h0 ' + str(Vars['h0']) + ' --Band ' + str(CreationBand) + ' --cosi ' + str(Vars['cosi']) + ' --psi ' + str(Vars['psi']) + ' --phi0 ' + str(Vars['phi0']) + ' --Freq ' + str(Vars['Finj']) + ' --Alpha ' + str(Vars['Alpha']) + ' --Delta ' + str(Vars['Delta']) + ' --IFO ' + str(Vars['IFO'][ifo]) + ' --outSFTbname ' + str(Vars['Out']) + ' --ephemDir ' + str(Vars['Ephem']) + ' --ephemYear ' + str(Vars['EphemYear']) + ' --f1dot ' + str(Vars['FDot']) + ' --noiseSqrtSh ' + str(Vars['Sh']**0.5) + ' --timestampsFile timestampsFile' + str(Vars['IFO'][ifo]) + UniqueID + endstring
    return(S)
                         
def GenDataString(beginstring,endstring,Vars):
    S = beginstring + ' --Freq ' + str(Vars['Fmin']) + ' --FreqBand ' + str(Vars['Band']) + ' --Alpha ' + str(Vars['Alpha']) + ' --Delta ' + str(Vars['Delta'])  + ' --ephemDir ' + str(Vars['Ephem']) + ' --ephemYear ' + str(Vars['EphemYear']) +  ' --dFreq ' + str(Vars['Res']) + ' --DataFiles \"' + str(Vars['Out']) + '/*" ' + ' --f1dot ' + str(Vars['FDot']) +  endstring
    return(S)

#' --f1dot ' + str(Vars['FDot']) + 
if __name__ == '__main__':
    main()


















