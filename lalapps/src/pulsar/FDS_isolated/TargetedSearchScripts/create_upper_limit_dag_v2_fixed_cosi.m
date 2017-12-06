close all
clear all

%%%%%%%%%%%%
%Basic Setup
%%%%%%%%%%%%

BASENAME = 'crab_H1H2L1_5_ul_fixed_cosi_psi'; %Name to be appended to all the output files
STARTINGDIR = '/archive/home/josephb/searches2/crab_H1H2L1_5_ul_fixed_cosi_psi/';  %Directory where the starting files and executables are
WORKINGDIR = '/usr1/josephb/upper_limit'; %Directory on the individual node where temporary files will be kept
CODEFILENAME = strcat([STARTINGDIR,'/','Upper_Limit_V2_cosi_psi.py']);  %Name of the python script to run
INIFILENAME = strcat([STARTINGDIR,'/','crab_H1H2L1_5_MC.ini']); %Name of the initialization file which sets the parameters the python script uses
LOGDIR = '/usr1/josephb/'; %Directory where the condor log files are placed

DAGFILENAME = strcat(BASENAME,'.dag'); %Name to be used for the dag file
SUBFILENAME = strcat(BASENAME,'.sub'); %Name to be used for the sub file

TOTALJOBS = 100; %Number of jobs to split into

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define line format for .dag file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

firstline = strcat(['# DAG for base_name = ',BASENAME,' \n']); 

secondline = strcat(['JOB A%d ',SUBFILENAME,'\n']);

outputformat = strcat(['VARS A%d JobID="%d" argList=" -j %d -S ',STARTINGDIR,' -W ',WORKINGDIR,' -p ',INIFILENAME,'"']);

%%%%%%%%%%%%%%%%%%%%%%%%%
%Create Condor .dag file
%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(DAGFILENAME) ~= 0 %Delete any existing file with the same name
    delete(DAGFILENAME)
end

fid = fopen(DAGFILENAME,'w'); %Open the file
    
%%%%%%%%%%%%%%%%
%Write .dag file
%%%%%%%%%%%%%%%%

fprintf(fid,firstline)
job = 0;
for x = 1:TOTALJOBS
    
    fprintf(fid,secondline,job);
    
    fprintf(fid,outputformat,[job,job,job]);
    fprintf(fid,'\n\n');
    job = job + 1;
end

fclose(fid) %Close the file

if exist(SUBFILENAME) ~= 0
    delete(SUBFILENAME);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%Create Condor .sub file
%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(SUBFILENAME,'w');
line1 = 'universe   = vanilla\n';i
line2 = strcat(['executable = ',CODEFILENAME,'\n']);
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