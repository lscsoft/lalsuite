%
% $Id$
%
%  Matlab script to get C(h_0)
%  
% Remember to add the path where this file is located:
% addpath /scratch/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type  getCofH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Detector = 'L1';

fileoutput = strcat(Detector,'kkoutput');
fid = fopen(fileoutput, 'w');

Directory1 = '/scratch/sintes/S2/Driver_allsky/';
Directory2 = '/scratch/sintes/S2/MC_allsky/';
subdir = 'L1_test/';

sufix = 'bandlist';
sufix2='.m.dat';

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
  addpath(filepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getting info from the first one

%MC_L1_0_000
%   nMonteCarlos=length(Ncount);
%   x=Ncount(:,2);
%   kkcount = find(x>Nmax(bandnumber));
%   CH = length(kkcount)/nMonteCarlos
%   h0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepre = strcat(filepath,'/dat/MC_');
filepre = strcat(filepre,Detector);
filepre = strcat(filepre,'_');

  nh0=10; 
  
for bandnumber=0:2
  filepre2 = strcat(filepre, int2str( bandnumber ) );

   for h0num=1:nh0
     string1 = sprintf('_00%d.m.dat',h0num-1);
     string2 = strcat(filepre2,string1);
     Ncount= load(string2);
     nMonteCarlos=length(Ncount);
     x=Ncount(:,2);
     kkcount = find(x>Nmax(bandnumber+1));
     CH(h0num) = length(kkcount)/nMonteCarlos;
   end
   CH
   bn=bandnumber+1;
   fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n',bn-1, ...
           BandList(bn,1),BandList(bn, 2), Nmax(bn), CH);
end 
 
 fclose(fid);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  
%%%%%%%%%%%%%%%%%%%%%%%
%load('./L1_test/MC_L1_0_009.m', '-mat')
%load('./L1_test/MC_L1_0_009.m', '-ascii')
%currentfile = './L1_test/MC_L1_0_003.m'; 
%[kk] = textread(currentfile,'%s')
% kk2= kk(6:length(kk)-1);
