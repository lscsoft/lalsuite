#!/usr/bin/python
# Resampling Test Code. This code is much more readable and reusable than its predecessor. Author - Pinkesh Patel
import sys,random,commands,math,shutil,os,re,string

def main():
    
    # Import the Configuration File into config
    if(len(sys.argv) < 4):
        print "Insufficient command line arguments "
        print " Usage is ",sys.argv[0]," <ConfigFile> Search_Name OutputFilename"
        print " Exiting ........."
        sys.exit(1)
    
    try:
        config = __import__(sys.argv[1]);
    except:
        print "Cannot import Config file '",sys.argv[1],"' exiting ...."
        sys.exit(1)

    UniqueID = ''

    # Compute Unique ID Number
    UniqueID +=  sys.argv[2]
    OutFile = sys.argv[3]
   
    try:
        OutFileObj = open(OutFile,'w')
    except:
        print "Cannot create file ",OutFile
        sys.exit(1)

    # Variables Structure/Dictionary
    Vars = {}

    try:
        Vars['SubFile'] = config.SubFile
    except:
        print "Submit file cannot be read"
        sys.exit(1)

    try:
        Vars['Alpha'] = config.Alpha
        if(Vars['Alpha'] > 2*math.pi or Vars['Alpha'] < 0):
            print "Alpha out of bounds\n"
            sys.exit(1)
    except:
        print "Cannot read Alpha\n"
        sys.exit(1)

    try:
        Vars['Delta'] = config.Delta
        if(Vars['Delta'] > math.pi/2.0 or Vars['Delta'] < -math.pi/2.0):
            print "Delta out of bounds\n"
            sys.exit(1)
    except:
        print "Cannot read Delta\n"
        sys.exit(1)

    try:
        Vars['FMin'] = config.FMin
    except:
        print "Cannot read FMin\n"
        sys.exit(1)

    try:
        Vars['FMax'] = config.FMax
    except:
        print "Cannot read FMax\n"
        sys.exit(1)

    if(Vars['FMin'] >= Vars['FMax']):
        print " Error: FMin >= FMax\n"
        sys.exit(1)

    try:
        Vars['GridType'] = config.GridType
    except:
        print "Cannot read GridType\n"
        sys.exit(1)
    
    try:
        Vars['MetricType'] = config.GridType
    except:
        print "Cannot read MetricType\n"
        sys.exit(1)
     
    try:
        Vars['Tau'] = config.Tau*3600.0*24*365
    except:
        print "Cannot read Spindown age \n"
        sys.exit(1)

    try:
        Vars['nMin'] = config.nMin
        Vars['nMax'] = config.nMax
    except:
        print "Cannot read Braking indices\n"
        sys.exit(1)

    try:
        Vars['Jobs'] = config.Jobs
    except:
        print "Cannot read Jobs\n"
        sys.exit(1)

    try:
        Vars['DataLocation'] = config.DataLocation
    except:
        print "No Data Location Found\n"
        sys.exit(1)

    # Number of Dirichlet terms used.
    try:
        Vars['Dterms'] = float(config.Dterms)
    except:
        print "Dterms cannot be read"
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

    # Reference Time in SSB
    try:
        Vars['refTime'] = float(config.refTime)
    except:
        print "refTime cannot be read"
        sys.exit(1)
    
    Step = (Vars['FMax'] - Vars['FMin'])/Vars['Jobs']

    # Generate command for each processor
    for i in range(Vars['Jobs']):
        # Calculate the Bounding Frequencies for each processor
        Vars['FMinNow'] = i*Step+Vars['FMin']
        Vars['FMaxNow'] = Vars['FMinNow'] + Step
        Vars['Band'] = Step

        # Calculate the maximum and minimum of FDot
        Vars['FDotMin'] = -Vars['FMaxNow']/Vars['Tau']
        Vars['FDotMax'] = -Vars['FMinNow']/6.0/Vars['Tau']
        Vars['FDotBand'] = Vars['FDotMax'] - Vars['FDotMin']

        # Calculate the maximum and minimum of F2Dot
        Vars['F2DotMin'] = 2.0 * min(Vars['FDotMin'],Vars['FDotMax'])**2.0/Vars['FMaxNow']
        Vars['F2DotMax'] = 7.0 * max(Vars['FDotMin'],Vars['FDotMax'])**2.0/Vars['FMinNow']
        Vars['F2DotBand'] = Vars['F2DotMax'] - Vars['F2DotMin']

        OutputLoudest = UniqueID + '_Loudest_' + str(i)
        OutputFStat = UniqueID + '_FStat_' + str(i)

        beginstring = ''
        endstring = ' --outputLoudest ' + OutputLoudest + ' -S --outputFstat ' + OutputFStat
        String = r'argList="' + GenDataString(beginstring,endstring,Vars) + r'"'
        OutFileObj.write("JOB " + UniqueID + str(i) + " " + Vars['SubFile'] + '\n')
        OutFileObj.write("VARS " + UniqueID + str(i) + " JobID=\"" + str(i) + "\" " + String + '\n')

    return(0)

def GenDataString(beginstring,endstring,Vars):
    S = beginstring + ' --Freq ' + str(Vars['FMinNow']) + ' --FreqBand ' + str(Vars['Band']) + ' --Alpha ' + str(Vars['Alpha']) + ' --Delta ' + str(Vars['Delta'])  + ' --refTime ' + str(Vars['refTime']) + ' --ephemDir ' + str(Vars['Ephem']) + ' --ephemYear ' + str(Vars['EphemYear']) + r' --DataFiles "' + str(Vars['DataLocation']) + r'/*" ' + ' --f1dot ' + str(Vars['FDotMin']) + ' --f2dot ' + str(Vars['F2DotMin']) + ' --f1dotBand ' + str(Vars['FDotBand'])  + ' --f2dotBand '+ str(Vars['F2DotBand']) +  endstring
    return(S)

if __name__ == '__main__':
    main()
