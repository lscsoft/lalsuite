% $Id$
%
%  Matlab script to get C(h_0)
%
% Remember to add the path where this file is located:
% addpath /local_data/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Detector = 'H1';

fileoutput = strcat(Detector,'outCH2');
fid = fopen(fileoutput, 'w');

fileinput = strcat(Detector,'h0newband');

DirectoryMC = '/local_data/badkri/S2-clean/MC_allsky/';

subdir = strcat(Detector, '_fine');
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

filepath = strcat(DirectoryMC,subdir);
  addpath(filepath);

filepre = strcat(filepath,'/MC_');
filepre = strcat(filepre,Detector);
filepre = strcat(filepre,'_');



for bandnumber=0:Nbands-1

  bn=bandnumber+1;  
  steph0= (h0max(bn)-h0min(bn))/(nh0-1);
  h0vect = h0min(bn):steph0:h0max(bn);
 
  basestring = strcat(filepre, int2str( bandnumber ) );

  mystring = strcat(basestring, '_nc');
  parstring = strcat(basestring, '_par');

  Ncount0 = load(mystring);
  parvals = load(parstring);  

  MC_FreqVals = parvals(1:length(Ncount0(:,1)),2);
  clear parvals

  % use doppler wing to veto frequencies 
  dopp = fmax(bn)*0.0001;
  vetoindices = find((mod(MC_FreqVals,1) > dopp & mod(MC_FreqVals,1) < 0.25-dopp) | (mod(MC_FreqVals,1) > 0.25+dopp & mod(MC_FreqVals,1) < 0.50-dopp) | (mod(MC_FreqVals,1) > 0.50+dopp &mod(MC_FreqVals,1) < 0.75-dopp) | (mod(MC_FreqVals,1) > 0.75+dopp & mod(MC_FreqVals,1) < 1.00-dopp)); 
  clear MC_FreqVals

  nMonteCarlos=length(vetoindices);

 
  fprintf(fid,'%d %d %d %d ', bn-1, fmin(bn), fmax(bn), Nmax(bn) );
  
   for h0num=1:nh0
     x=Ncount0(vetoindices, h0num+1);   
     kkcount = find(x>Nmax(bn));
     CH(h0num) = length(kkcount)/nMonteCarlos;
     fprintf(fid,' %d %d ', h0vect(h0num), CH(h0num) );
   end
     fprintf(fid,' \n');
   CH
   
end

 fclose(fid);
