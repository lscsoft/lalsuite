#!/usr/bin/env python2.1
"""
Monte Carlo driver of the FDS pipeline.
"""

# Import standard modules to the python path
import string, os,sys, math, random, shutil, time, os.path
# Python path which runs the Xavier's scripts
pythonexe='/usr/bin/python2.1'

# This and Xavier's scripts are not compatible with an older version python.... 
if sys.version[:3] < '2.1':
    print 'python version >= 2.1 is required to run the script.'
    print 'You are using:'
    print sys.version
    sys.exit(1)

print os.uname()[1]
print time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime())


##---------------------------------------------------------------
# Class and function definition 
##---------------------------------------------------------------

# make a python configuration file from a parameters dictionary.
def makeinifile(filename='inifile.ini',paramsdictionary='paramsdict'):
    fp=open(filename,'w')
    for __thiskey in paramsdictionary.keys():
        fp.write('['+__thiskey+']\n')
        __dict_in_this_key=paramsdictionary[__thiskey]
        for __itemkey in __dict_in_this_key.keys():
            fp.write(__itemkey)
            fp.write('=')
            __item_in_itemkey=__dict_in_this_key[__itemkey]
            fp.write(str(__item_in_itemkey)+'\n')
        fp.write('\n')
    fp.close()

def myrm(directory):
    if directory.isspace() == 1:
        print 'Warning: wrong directory name', directory
        sys.exit(1)
    for c in string.whitespace:        
        if directory.find(c) != -1:
            print 'Warning: wrong directory name', directory
            sys.exit(1)
    if os.path.isdir(directory) and len(directory) >= 1:
        directory=os.path.normpath(directory)
        __com='rm -rf '+directory+'/*'
        print 'running:', __com
        os.system(__com)
    sys.stdout.flush()
    sys.stderr.flush()

##---------------------------------------------------------------
# General parameters, directories
##---------------------------------------------------------------
# mc_id is the seed for the random numbers
mc_id=sys.argv[1]
print mc_id

homedir='/home/yousuke/mchar_pipeline/'  #Where the result files should be copied to.
parentdir='/scratch/tmp/yousuke'         #Central working dir on each node.
codesdir=parentdir+'/codes'              #Where all the codes and timestamps shuld be.
ephemerisdir=parentdir+'/ephemeris'      
startingdir=parentdir+'/AllSky.'+mc_id   #Working dir. 

myrm(startingdir)
try: os.mkdir(startingdir)
except OSError, err:
    import errno
    print "Warning:", err
os.chdir(startingdir)


sys.stdout.flush()
sys.stderr.flush()

##---------------------------------------------------------------
# Copy necessary files to the working directory
##---------------------------------------------------------------
ifos=['LLO','LHO']

if len(ifos)>2:
    print "The scripts can handle only two IFOs at once."
    sys.exit(1)

# for a backward-compatibility
ifodict={'GEO':0,
         'LLO':1,
         'LHO':2}


necessities=['lalapps_ComputeFStatistic',
             'lalapps_SemiAnalyticF',
             'makefakedata_v2',
             'lalapps_FstatShapeTestLAL',
             'lalapps_polka',
             'lalapps_makeInvetofile',
             'TSLHO',
             'TSLLO',
             'allsky_pulsar_pipe_coincidence.py',
             'allsky_pulsar_pipe_search.py',
             'allsky_pulsar_pipe_upper_limit.py']
ephemerisfiles=['earth00-04.dat',
                'sun00-04.dat',
                'earth03.dat',
                'sun03.dat']

for code in necessities:
    shutil.copy(codesdir+'/'+code,startingdir)
for file in ephemerisfiles:
    shutil.copy(ephemerisdir+'/'+file,startingdir)


# Name of pipeline parameters file
paramsfile='allsky_pulsar_pipe_monte.ini'



##---------------------------------------------------------------
# loop
##---------------------------------------------------------------
mcparamsfile='mcparams.'+mc_id
resultsfile='results.'+mc_id
fres=open(resultsfile,'a')
Nmonte=20
random.seed(int(mc_id))
monte_counter=0
while monte_counter < Nmonte:
    print monte_counter
##---------------------------------------------------------------
# Make configuration file
##---------------------------------------------------------------

## randomly choose parameters
    fsignal    = random.uniform(  160.0,       460.0 )
    fdotsignal = random.uniform( -1.0e-12,     1.0e-12 )
##    fdotsignal = 0.0
    asignal    = random.uniform(  0.0,         2.0*math.pi )
    dsignal    = math.asin( 2.0*random.random() - 1.0 )
    psisignal  = random.uniform( -math.pi/4.0, math.pi/4.0 )
    phisignal  = random.uniform(  0.0,         2.0*math.pi )
    cosiota    = random.uniform( -1.0,         1.0 )
# We use lalapps_Makefakedata in this script....
    h0         = random.uniform(  4.0e-23,      4.0e-22 )

    aPlussignal = h0*(1.0+cosiota**2)/2.0
    aCrossignal = h0*cosiota
# Multiply 1e19 because we use makefakedata_v2 in the Xavier's scripts...
    mch0=2*h0*1e19
    dh0=0.4*h0*1e19





# for mc makefakedata injection
    mfd_fband=8.0
    mfd_fstart=math.floor(fsignal-mfd_fband/2.0)

## for search scripts
    fband=0.1
    fstart=fsignal-fband/2.0
    skyband=0.28
    randalpha=random.uniform(-0.01,0.01)
    randdelta=random.uniform(-0.01,0.01)
    if math.cos(dsignal)/skyband < 1.0/2.0/math.pi:
        aband=2*math.pi
        dband=skyband
    else:
        aband=skyband/math.cos(dsignal)
        dband=skyband
    astart=asignal-aband/2.0 + randalpha
    dstart=dsignal-dband/2.0 + randdelta

        
    phisignals=[phisignal]
    fsignals=[fsignal]
        
    ic=0
    sftdir=[]
    ifoname=[]
    datadir=[]
    nifo=[]
    sifo=[]
    timestamps=[]
    gps_start=[]
    nanogps_start=[]
    for ifo in ifos:
        datadir.append(parentdir+'/sfts/S2-'+ifo+'/')
        sftdir.append(parentdir+'/sfts/tmp_'+ifo+mc_id+'/')
        myrm(sftdir[ic])
        try: os.mkdir(sftdir[ic])
        except OSError, err:
            import errno
            print "Warning:", err
        timestamps.append('TS'+ifo)
        nifo.append(ifodict[ifo])
        sifo.append(ifo)
        fp=open(timestamps[ic],'r')
        gpspair=string.split( fp.readline() )
        gps_start.append( gpspair[0] )
        nanogps_start.append( gpspair[1] )
        fp.close()
        ic=ic+1
        

    timedifference=float(gps_start[1])-float(gps_start[0])
    nanotimedifference=float(nanogps_start[1])-float(nanogps_start[0])
    phidelay=math.fmod( fsignal*timedifference + fdotsignal/2.0*timedifference**2 , 1)
    phinanodelay=math.fmod( fsignal*nanotimedifference + fdotsignal/2.0*nanotimedifference**2 , 10.0**9)
   
    phisignals.append( phisignal+2*math.pi*(phidelay+phinanodelay) )
                
    fsdelay=math.fmod( fdotsignal*timedifference , 1)
    fsnanodelay=math.fmod( fdotsignal*nanotimedifference , 10.0**9)
    fsignals.append( fsignal + fsdelay+ fsnanodelay)
    
    allskysearchini={'fstat-params':{'start_freq': fstart,
                                     'freq_band' : fband,
                                     'df' :  0.00000347222222222222,
                                     'ifo1' : nifo[0],
                                     'ifo2' : nifo[1],
                                     'a_search' : string.join(['-a', str(astart),'-z', str(aband), '-l 0.02 --gridType 1']),
                                     'd_search' : string.join(['-d',str(dstart),'-c',str(dband), '-g 0.02']),
                                     'data1' : sftdir[0],
                                     'data2' : sftdir[1],
                                     'Fth' : 20.0
                                     },
                     'polka-params':{'freq_window' : 0.001,
                                     'alpha_window' : 0.02,
                                     'delta_window' : 0.02,
                                     'coincidence_band' : 1.0
                                     },
                     'chisq-params':{'Fth_chisq' : 50.0,
                                     'chisq_points' : 50,
                                     'ts1': timestamps[0],
                                     'ts2': timestamps[1],
                                     'sifo1' : sifo[0],
                                     'sifo2' : sifo[1]
                                     },
                     'mc-params':{'Ninj':1000,
                                  'h0':mch0,
                                  'dh0':dh0,
                                  'c0':0.95,
                                  'tol':0.03,
                                  'gps_start1': gps_start[0],
                                  'gps_start2': gps_start[1]
                                  }
                     }
    
    makeinifile(filename=paramsfile,paramsdictionary=allskysearchini)
                
##---------------------------------------------------------------
## Generation a signal
##---------------------------------------------------------------

                
    ic=0
    for ifo in ifos:
        sftdirname=sftdir[ic]
        sftbase=sftdirname+'/SFT_TMP'+ifo
        makefakedataargslist=map(str,[codesdir+'/lalapps_Makefakedata',
                                      '--outSFTbname',sftbase,
                                      '--detector',ifo,
                                      '--ephemDir',ephemerisdir,
                                      '--ephemYear','03',
                                      '--timestampsFile',timestamps[ic],
                                      '--Tsft',1800,
                                      '--nTsft',20,
                                      '--fmin',mfd_fstart,
                                      '--Band',mfd_fband,
                                      '--longitude',asignal,
                                      '--latitude',dsignal,
                                      '--aPlus',aPlussignal,
                                      '--aCross',aCrossignal,
                                      '--psi',psisignal,
                                      '--phi0',phisignals[ic],
                                      '--f0',fsignals[ic],
                                      '--f1dot',fdotsignal,
                                      '--noiseSFTs','\"'+datadir[ic]+"CAL_SFT.*\""
                                      ])
        makefakedataargs=string.join(makefakedataargslist)
        os.system(makefakedataargs)
        ic=ic+1




##---------------------------------------------------------------
## Save the signal parameters and perfect matched F.
##---------------------------------------------------------------

## Perfect matched F stat
    perfectmatchF='perfectmatch.'+mc_id
    fp=open(perfectmatchF,'a')
    for ic in range(len(ifos)):
        fstatoutfile='fstatout'+ifos[ic]
        cfsargs = string.join([startingdir+'/lalapps_ComputeFStatistic',
                               '-f', str(fsignals[ic]),
                               '-a', str(asignal),
                               '-d', str(dsignal),
                               '-I', ifos[ic],
                               '-E', ephemerisdir,
                               '-D', sftdir[ic],
                               '--outputFstat', fstatoutfile
                               ])
        os.system(cfsargs)        
        fin=open(fstatoutfile,'r')
        fp.write(fin.read())
        fin.close()
    fp.close()
                               
    params=(int(mc_id),monte_counter,fsignal, fdotsignal, asignal, dsignal, h0, cosiota, psisignal, phisignal)
    fdparams=open(mcparamsfile,'a')
    format='%5d %5d %15.9f %10.6g %10.6f %10.6f %10.6g %10.6f %10.6f %10.6f'
    print >> fdparams, format %  params
    fdparams.close()

        

##jobid=0 because we determine the starting frequency directly in the .ini file. 
    jobid=0
##---------------------------------------------------------------
## Run the search script
##---------------------------------------------------------------
    print 'running: search script'
    localworkdir=startingdir+'/search'
    myrm(localworkdir)
    allskysearch = string.join([pythonexe,
                                './allsky_pulsar_pipe_search.py',
                                '--job-id',str(jobid),
                                '--starting-dir',startingdir,
                                '--local-work-dir',localworkdir,
                                '--params-file',paramsfile])
    
    os.system(allskysearch)
    sys.stdout.flush()
    sys.stderr.flush()
    

##---------------------------------------------------------------
## cat and rename the resulting file 
##---------------------------------------------------------------
    
#    print "running: rename the resulting files"
    filename=startingdir+'/Fstats-1-'+str(fstart)+'.gz'
    try: os.rename(filename, startingdir+'/FstatsJOINED-1-'+str(fstart)+'.gz')
    except OSError, err:
        import errno
        print "Warning:", err, filename
        sys.exit(1)
    filename=startingdir+'/Fstats-2-'+str(fstart)+'.gz'
    try: os.rename(filename, startingdir+'/FstatsJOINED-2-'+str(fstart)+'.gz')
    except OSError, err:
        import errno
        print "Warning:", err, filename
        sys.exit(1)

    sys.stdout.flush()
    sys.stderr.flush()

##---------------------------------------------------------------
## Run the coincidence script 
##---------------------------------------------------------------

    print 'running: coincidence script'
    
    localworkdir=startingdir+'/coincidence'
    myrm(localworkdir)
    allskycoincidence = string.join([pythonexe,
                                     './allsky_pulsar_pipe_coincidence.py',
                                     '--job-id',str(jobid),
                                     '--starting-dir',startingdir,
                                     '--local-work-dir',localworkdir,
                                     '--params-file',paramsfile])
    
    os.system(allskycoincidence)

    filename=startingdir+'/results_out-'+str(fstart)
    if os.path.isfile(filename):
        fin=open(filename,'r')
        for line in fin.xreadlines():
            if line[0]=='%': continue
            print >> fres, '%5d %5d \t' % (int(mc_id), monte_counter), 
            fres.write(line)
        fin.close()

    os.system('rm -f polka_out-*')
    os.system('rm -f results_out-*')
    os.system('rm -f FstatsJOINED-*')
    
    sys.stdout.flush()
    sys.stderr.flush()
    monte_counter=monte_counter+1

fres.close()
shutil.copy(resultsfile, homedir)
shutil.copy(mcparamsfile, homedir)
shutil.copy(perfectmatchF, homedir)

for ic in range(len(ifos)):
    myrm(sftdir[ic])
os.chdir(homedir)    
myrm(startingdir)

print time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime())
