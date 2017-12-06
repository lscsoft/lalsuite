close all
clear all

%%%%%%%%%%%%
%Basic Setup
%%%%%%%%%%%%
BASENAME = 'crab_H1H2L1_5'; %Name to be appended to all the output files
STARTINGDIR = '/archive/home/josephb/searches2/crab_H1H2L1_5/'; %Directory on the cluster where the code will be excuted from and results saved
CODEDIR = '/archive/home/josephb/searches/';  %Directory containing the ComputeFStatistic_v2 code
CODEFILENAME = 'ComputeFStatistic_v2'; %Name of the ComputeFStatistic_v2 executable
LOGDIR = '/usr1/josephb/'; %Location to save the condor log files

%%%%%%%%%%%%%%%%%%%%%%%
%Known Parameters of source
%%%%%%%%%%%%%%%%%%%%%%%
F0 = 29.7778867428; %Frequency for the pulsar - from Radio Ephemeris
F1 = -372940.45e-15; %1st Spindown (dot F) for the pulsar - from Radio Ephemeris 
F2 = 1.2426e-20; %2nd Spindown (dot dot F) for the pulsar - from Radio Ephemeris

PEPOCH = 813369613; %Period epoch - time from which the Radio Ephemeris comes in GPS time

POSEPOCH = 497318410;  %Positional Epoch - time from which the Crab position comes
Current_EPOCH =  815201292; %Time at which the search begins (start time of first sft)
pmra = -1.9972e-15; %Change in Right Ascension over time - rad/sec
pmdec = 1.0754e-15; %Change in declincation over time - rad/sec

Current_Freq = F0 + F1 * (Current_EPOCH - PEPOCH) + 0.5 * F2 * (Current_EPOCH - PEPOCH)^2;  %Calculated frequency at the time of the start of the search
Current_Spindown = F1 + F2 * (Current_EPOCH - PEPOCH); %Calculated 1st Spindown (dot F) at the time of the start of the search 

dRA = pmra*(Current_EPOCH-POSEPOCH);  %Calculated difference in Right Ascension from POSEPOCH to start of search
dDEC = pmdec*(Current_EPOCH-POSEPOCH); %Calculated difference in declination from POSEPOCH to start of search

RA = 1.45967501813819 + dRA; %Final Right Ascension in radians
DEC = 0.38422487307437 + dDEC; %Final Declination in radians

%%%%%%%%%%%%%%%%%%%%
%Details of data
%%%%%%%%%%%%%%%%%%%%
NUMBEROFSFT = 10005  %Total number of SFTs for H1 = 8967, Total number of SFTs for L1 = 6828, Total number of SFTs for H2 = 10005 - use the largest number in the search

RUNTIME = 840336581 - 815201292; %Seconds between start time of first SFT and end time of last SFT

SFTLENGTH = 1800;

COHTIME = (SFTLENGTH*NUMBEROFSFT);

SFTDIRECTORY = '/archive/home/josephb/SFTs/crabSFTsMerged_C03_Ht_removed'; %Directory containing SFT files or links to the SFT files
SFTNAME = '/*'; %Pattern for the SFT files
SFTPATTERN = strcat([SFTDIRECTORY,SFTNAME]); %Combined Directory and Pattern indicating the SFTs
EPHEMERISDIR = '/archive/home/josephb/share/lal'; %Directory containing the emphemeris data
EPHEMERISYEAR = '03-06'; %Year of the emphemeris data
DAGFILENAME = strcat(BASENAME,'.dag'); %Name to be used for the dag file
SUBFILENAME = strcat(BASENAME,'.sub'); %Name to be used for the sub file
FSTATFILENAME = strcat(BASENAME,'_result'); %Name for the file containing all the Fstat values above the TWOFCUTOFF
LOUDESTFILENAME = strcat(BASENAME,'_loudest'); %File to put the single largest template in.
NODEFILENAME = strcat(BASENAME,'_nodeA'); %Name for the condor .out and .err files
INIFILENAME = strcat(BASENAME,'_MC.ini'); %Name for MC injection initialization file

%%%%%%%%%%%%%%%%%%%%
%Details of Search
%%%%%%%%%%%%%%%%%%%%
TOTALJOBS = 100; %Must be greater than 1

TWOFCUTOFF = '20'; %Smallest 2F value to be kept

BANDSIZE = 1.2e-2; %Overall width in Frequency space to be searched
STARTFREQ = 2*Current_Freq - BANDSIZE/2; %Frequency at which the search will start and move upwards (I.e. from STARTFREQ to STARTFREQ + BANDSIZE)

MISMATCH = 0.15


%The following are derived from "Parameter space metric for combined diurnal and orbital motion" by Ian Jones, Ben Owen, and David Whitbeck
%which can be found at: http://www.lsc-group.phys.uwm.edu/cgi-bin/enote.pl?nb=puls5directed&action=view&page=5
%and also from "Monte-Carlo tests of ComputeFStatistic for 4D parameter-spaces" by Reinhard Prix
%which can be found at: http://www.lsc-group.phys.uwm.edu/cgi-bin/enote.pl?nb=puls3knownpulsardemod&action=view&page=2

dFreq = 2*sqrt(3*MISMATCH/(pi^2*RUNTIME^2)); %Expected Spacing between templates in Frequency space
dfdot = 2*sqrt(4*5*9*MISMATCH/(pi^2*RUNTIME^4)); %Expected Spacing between templates in 1st Spindown space
dfdotdot = 2*2*sqrt(MISMATCH/((1/28)*pi^2*RUNTIME^6)); %Expected Spacing between templates in 2nd Spindown space

dFreq_step = BANDSIZE/TOTALJOBS; %Size of steps in frequency space for each job
dFreq_step_band = dFreq_step + dFreq; %Size of band in frequency space for each job (step size plus 1 template spacing ensures between 1 and 2 frequency templates worth of overlap)

FDOTBAND = 3.0e-13; %Overall width in 1st Spindown space to be searched
fdotStart = 2*Current_Spindown-FDOTBAND/2; %1st Spindown starting point
fdotEnd = fdotStart + FDOTBAND;

FDOTDOTBAND = 1e-23; %Overall width in 2nd Spindown space to be searched
fdotdotStart = 2*F2-FDOTDOTBAND/2; %2nd Spindown starting point
fdotdotEnd = fdotdotStart + FDOTDOTBAND;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define line format for .dag file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

firstline = strcat(['# DAG for base_name = ',BASENAME,' \n']); %First line

secondline = strcat(['JOB A%d ',SUBFILENAME,'\n']); %1st Job line

outputformat = strcat(['VARS A%d JobID="%d" argList=" -a %15.10f -d %15.10f -f %15.10e -b %15.10e -s %15.10e -m %15.10e --f2dot %15.10e --f2dotBand %15.10e --dFreq ',num2str(dFreq), ' --df1dot ',num2str(dfdot), ' --df2dot ',num2str(dfdotdot) ,' --DataFiles=',SFTPATTERN, ' -E ',EPHEMERISDIR,' -y ',EPHEMERISYEAR,' -F ', TWOFCUTOFF,' -o ',NODEFILENAME,'%d --outputLoudest ',LOUDESTFILENAME,'%d --outputFstat ',FSTATFILENAME,'%d"']); %2nd job line - defines most of the search details





%%%%%%%%%%%%%%%%%%%%%%%%%
%Create Condor.dag file
%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(DAGFILENAME) ~= 0 %Delete any existing file with same name 
    delete(DAGFILENAME)
end

fid = fopen(DAGFILENAME,'w'); %Open the file

%%%%%%%%%%%%%%%%%
%Write .dag file
%%%%%%%%%%%%%%%%%

fprintf(fid,firstline)

job = 0;
for y = 0:(TOTALJOBS - 1)

    startfreq_job = STARTFREQ + dFreq_step*y; %Calculate starting frequency for this job
    
    fprintf(fid,secondline,job);
        
    fprintf(fid,outputformat,[job,job,RA,DEC,startfreq_job,dFreq_step_band,fdotStart,FDOTBAND,fdotdotStart,FDOTDOTBAND,job,job,job]);
    fprintf(fid,'\n\n');
    job = job + 1;
end

fclose(fid) %close the file

if exist(SUBFILENAME) ~= 0
    delete(SUBFILENAME)
end

%%%%%%%%%%%%%%%%%%%%%%%%
%Create Condor .sub file
%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(SUBFILENAME,'w');
line1 = 'universe   = standard\n';
line2 = strcat(['executable = ',CODEDIR,CODEFILENAME,'\n']);
line3 = 'output     = nodeA.out.$(JobID)\n';
line4 = 'error      = nodeA.err.$(JobID)\n';
line5 = strcat(['log        = ',LOGDIR,BASENAME,'.log\n']); 
line6 = '\n';
line7 = 'arguments  = $(arglist)\n';
line8 = 'queue\n';

fprintf(fid,line1);
fprintf(fid,line2);
fprintf(fid,line3);
fprintf(fid,line4);
fprintf(fid,line5);
fprintf(fid,line6);
fprintf(fid,line7);
fprintf(fid,line8);

fclose(fid)

%%%%%%%%%%%%%%%%%%%%%%%%
%Details of Monte Carlo Search
%%%%%%%%%%%%%%%%%%%%%%%%
NIFO = 3;
IFO1 = 'H1';
IFO2 = 'L1';
IFO3 = 'H2';
DATA1 = strcat([SFTDIRECTORY,'/H*H1*']);
DATA2 = strcat([SFTDIRECTORY,'/L*L1*']);
DATA3 = strcat([SFTDIRECTORY,'/H*H2*']);

ADD_BAND = 0.3;
JOB_GROUP_STEP = 0;
JOB_BAND = BANDSIZE;
LOUDEST_DIR = STARTINGDIR;
NFILES = TOTALJOBS;
INJBAND = 0.3;

freq_window = 5*dFreq;
fdot_window = 5*dfdot;
fdotdot_window = FDOTDOTBAND;

NINJ = 50;
NINJMAX = 500;
H0 = 1.5e-24;
DH0 = 5e-25;
C0 = 0.95;


%%%%%%%%%%%%%%%%%%%%%%%%%
%Create Monte Carlo ini file
%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen(INIFILENAME,'w');
iline(1) = {'[fstat-params]'};
iline(2) = {strcat(['alpha = ',num2str(RA,15)])};
iline(3) = {strcat(['delta = ',num2str(DEC,15)])};
iline(4) = {strcat(['freq = ',num2str(STARTFREQ,15)])};
iline(5) = {'[data-params]'};
iline(6) = {strcat(['add_band = ', num2str(ADD_BAND,15)])};
iline(7) = {strcat(['job_group_step = ',num2str(JOB_GROUP_STEP,15)])};
iline(8) = {strcat(['job_band = ',num2str(JOB_BAND,15)])};
iline(9) = {'[files-names]'};
iline(10) = {strcat(['loudest_dir = ',STARTINGDIR])};
iline(11) = {strcat(['loudest_files = ',LOUDESTFILENAME])};
iline(12) = {strcat(['Nfiles = ',num2str(NFILES,15)])};
iline(13) = {'[ifo-params]'};
iline(14) = {strcat(['Nifo = ',num2str(NIFO,15)])};
iline(15) = {strcat(['ifo1 = ',IFO1])};
iline(16) = {strcat(['ifo2 = ',IFO2])};
iline(17) = {strcat(['ifo3 = ',IFO3])};
iline(18) = {strcat(['data1 = ',DATA1])};
iline(19) = {strcat(['data2 = ',DATA2])};
iline(20) = {strcat(['data3 = ',DATA3])};
iline(21) = {'[inj-params]'}
iline(22) = {strcat(['ref_time = ',num2str(Current_EPOCH,15)])};
iline(23) = {strcat(['fdot_start = ',num2str(fdotStart,15)])};
iline(24) = {strcat(['fdot_end = ',num2str(fdotEnd,15)])};
iline(25) = {strcat(['fdotdot_start = ',num2str(fdotdotStart,15)])};
iline(26) = {strcat(['fdotdot_end = ',num2str(fdotdotEnd,15)])};
iline(27) = {strcat(['inj_band = ', num2str(INJBAND,15)])};
iline(28) = {'[search-params]'};
iline(29) = {strcat(['dFreq = ',num2str(dFreq,15)])};
iline(30) = {strcat(['dfdot = ',num2str(dfdot,15)])};
iline(31) = {strcat(['dfdotdot = ',num2str(dfdotdot,15)])};
iline(32) = {strcat(['freq_window = ',num2str(freq_window,15)])};
iline(33) = {strcat(['fdot_window = ',num2str(fdot_window,15)])};
iline(34) = {strcat(['fdotdot_window = ',num2str(fdotdot_window,15)])};
iline(35) = {'[mc-params]'}
iline(36) = {strcat(['Ninj = ',num2str(NINJ,15)])};
iline(37) = {strcat(['Ninjmax = ',num2str(NINJMAX,15)])};
iline(38) = {strcat(['h0 = ',num2str(H0,15)])};
iline(39) = {strcat(['dh0 = ',num2str(DH0,15)])};
iline(40) = {strcat(['c0 = ',num2str(C0,15)])};

for a = 1:40
    line = strcat([iline{a},'\n']);
    fprintf(fid,line);
end

fclose(fid)
%%%%%%%%%%%%%%%%%%%%%%%%%
%Create matlab Log file
%%%%%%%%%%%%%%%%%%%%%%%%%
LOGFILENAME = strcat(['save ',BASENAME,'_log.mat']);
eval(LOGFILENAME)

%%%%%%%%%%%%%%%%%%%%%%%%%
%
