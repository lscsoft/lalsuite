% 
% $Id$
%
%  Matlab script to count the number of templates used in a 1 Hz band
%  for the full sky search
%  This will change for the different frequency bands.
 
% Remember to add the path where this file is located: 
% addpath /scratch/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type Count1HzTemplates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Directory = '/scratch/sintes/S2/Driver_allsky/';
Detector = 'L1';
%Detector = 'H1';
DirName = '/skypatch_';
prefix ='/freq_';
sufix = '_histo';

NskyPatches = 23;

bandnumber =0;  % the band of interest

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutDir = strcat(Directory,Detector);
Ntemplates =0;
  
  for skypatch = 1:NskyPatches;
     ThisDir = strcat(DirName, int2str( skypatch ));
     ThisDir = strcat(OutDir,ThisDir);
     file = strcat(prefix, int2str( bandnumber ));
     file = strcat(file, sufix);
     file = strcat(ThisDir,file);

     output=load(file);
     partialSum = sum( output(:,2));
     Ntemplates =Ntemplates + partialSum;
  end
  
Ntemplates
