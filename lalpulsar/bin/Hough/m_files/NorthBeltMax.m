%
%
%  Matlab script to get the maximum numbercount values of the Hough map results
%  for the different frequency bins.
%  It will look into all the different patches around the north belt

% Remember to add the path where this file is located:
% addpath /scratch/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type NorthBeltMax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Directory = '/scratch/sintes/S2/Driver_allsky/';
%Detector = 'L1';
Detector = 'H1';
DirName = '/skypatch_';
prefix ='/freq_';
sufix = '_stats';

NskyPatches = 23;
NspinValues = 11;

% %%%%%%%%%%%%%%%For the north belt
Skypatch_ini=2;
Skypatch_end=8;

Nbands = 199;
%Nbands = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutDir = strcat(Directory,Detector);

FreqValues =[];
NoMaxValues=[];

for filenumber = 0:Nbands; %the  current frequency band

  filenumber

%%%% The 1st sky patch. Initializing
  %%%skypatch = 1;
  skypatch = Skypatch_ini;

     ThisDir = strcat(DirName, int2str( skypatch ));
     ThisDir = strcat(OutDir,ThisDir);
     file = strcat(prefix, int2str( filenumber ));
     file = strcat(file, sufix);
     file = strcat(ThisDir,file);

     output=load(file);
     xFreq = output(:,8);
     xMax = output(:,4);
     clear output;

     Length = length(xMax);
     freq = xFreq(1:NspinValues:Length)';
     XX =  reshape(xMax, NspinValues, Length/NspinValues );
     TheMax= max(XX);

%%%% The remaining sky patches


  for skypatch = (Skypatch_ini +1):Skypatch_end;
     ThisDir = strcat(DirName, int2str( skypatch ));
     ThisDir = strcat(OutDir,ThisDir);
     file = strcat(prefix, int2str( filenumber ));
     file = strcat(file, sufix);
     file = strcat(ThisDir,file);

     output=load(file);
     xMax = output(:,4);
     Length = length(xMax);
     XX =  reshape(xMax, NspinValues, Length/NspinValues );
     maxi= max(XX);
     TheMax= max(TheMax, maxi);
  end

  FreqValues = [FreqValues,freq];
  NoMaxValues  = [NoMaxValues,TheMax];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  plot(FreqValues, NoMaxValues, 'g.');
  xlabel ('Frequency [Hz]');
  ylabel ( 'Maximum number count');
  grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save NorthMaxResults FreqValues NoMaxValues
