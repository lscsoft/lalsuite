%
% $Id$
%
%  Matlab script to get C(h_0)
%
% Remember to add the path where this file is located:
% addpath /scratch/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type  H1getCofH_coarse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Detector = 'H1';

fileoutput = strcat(Detector,'kkoutput');
fid = fopen(fileoutput, 'w');

Directory1 = '/scratch/sintes/S2_clean/Driver_allsky/';
Directory2 = '/scratch/sintes/S2_clean/MC_allsky/';
subdir = 'H1_coarse/';

sufix = 'bandlist';
	%sufix2='.m.dat';

file = '/AllSkyMaxResults.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the  FreqValues and MaxValues

OutDir1 = strcat(Directory1,Detector);
file = strcat(OutDir1,file);
load(file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getting the different bands

file2 = strcat(Detector,sufix);
file2 = strcat(Directory2,file2);
BandList = load(file2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtaining the loudest event per band

Nbands = length(BandList);

for bandnumber = 1:Nbands; %the  current frequency band
   fmin = BandList(bandnumber, 1);
   fmax = BandList(bandnumber, 2);

   indices = find (FreqValues >= fmin & FreqValues <=fmax);
   Nmax(bandnumber) = max( MaxValues(indices));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear FreqValues  MaxValues fmin fmax indices bandnumber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = strcat(Directory2,subdir);
filepre = strcat(filepath,'/MC_');
filepre = strcat(filepre,Detector);
filepre = strcat(filepre,'_');

  nh0=10;
  
for bandnumber=0:Nbands-1
  bn=bandnumber+1;
  mystring = strcat(filepre, int2str( bandnumber ) );
  mystring = strcat(mystring, '_nc');
  Ncount0 = load(mystring);
  nMonteCarlos=length(Ncount0);

   for h0num=1:nh0
     x=Ncount(:, h0num+1);
     kkcount = find(x>Nmax(bn));
     CH(h0num) = length(kkcount)/nMonteCarlos;
   end
   CH
   
   fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n',bn-1, ...
           BandList(bn,1),BandList(bn, 2), Nmax(bn), CH);
end

 fclose(fid);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
