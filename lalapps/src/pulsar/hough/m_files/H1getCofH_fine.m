% $Id$
%
%  Matlab script to get C(h_0)
%
% Remember to add the path where this file is located:
% addpath /scratch/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Detector = 'H1';

fileoutput = strcat(Detector,'outCH2');
fid = fopen(fileoutput, 'w');

fileinput = strcat(Detector,'h0newband');

Directory1 = '/scratch/sintes/S2/Driver_allsky/';
Directory2 = '/scratch/sintes/S2/MC_allsky/';

subdir = 'H1_fine/';
sufix2='.m.dat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the different bands, loudest event per band
% and new h0 interval to refine

nh0=10;
h0inFile = load(fileinput);

fmin=h0inFile(:,2);
fmax=h0inFile(:,3);
Nmax=h0inFile(:,4);
h0min=h0inFile(:,5);
h0max=h0inFile(:,6);

Nbands =length(fmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = strcat(Directory2,subdir);
  addpath(filepath);

filepre = strcat(filepath,'/dat/MC_');
filepre = strcat(filepre,Detector);
filepre = strcat(filepre,'_');

for bandnumber=0:Nbands-1

  bn=bandnumber+1;  
  steph0= (h0max(bn)-h0min(bn))/(nh0-1);
  h0vect = h0min(bn):steph0:h0max(bn);
  
  filepre2 = strcat(filepre, int2str( bandnumber ) );

  fprintf(fid,'%d %d %d %d ', bn-1, fmin(bn), fmax(bn), Nmax(bn) );
  
   for h0num=1:nh0
     string1 = sprintf('_00%d.m.dat',h0num-1);
     string2 = strcat(filepre2,string1);
     Ncount= load(string2);
     nMonteCarlos=length(Ncount);
     x=Ncount(:,2);
     kkcount = find(x>Nmax(bn));
     CH(h0num) = length(kkcount)/nMonteCarlos;
     fprintf(fid,' %d %d ', h0vect(h0num), CH(h0num) );
   end
     fprintf(fid,' \n');
   CH
   
  %%%  fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n',bn-1, ...
  %           BandList(bn,1),BandList(bn, 2), Nmax(bn), CH);
end

 fclose(fid);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
