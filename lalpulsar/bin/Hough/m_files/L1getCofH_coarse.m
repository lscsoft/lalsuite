%
%
%  Matlab script to get C(h_0)
%
% Remember to add the path where this file is located:
% addpath /local_data/sintes/CVSDIR/waves/people/sintes/PULSAR/CODES/m_files/
%   To run simply type  L1getCofH_coarse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Detector = 'H1';

fileoutput = strcat(Detector,'kkoutput');
fid = fopen(fileoutput, 'w');

% directory with results of driver
Directory1 = '/local_data/badkri/S2-clean/Driver_allsky/';
% file with driver output
file = '/AllSkyMaxResults.mat';


% directory with results of coarse MC run
Directory2 = '/local_data/badkri/S2-clean/MC_allsky/';
subdir = strcat(Detector, '_coarse');
sufix = 'bandlist';
%sufix2='.m.dat';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the  FreqValues and MaxValues
% look into the output of the driver
% this loads {Directory1}/{detector}/{file}

OutDir1 = strcat(Directory1,Detector);
file = strcat(OutDir1,file);
load(file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the different bands
% this loads {Directory2}/{detector}bandlist

file2 = strcat(Detector,sufix);
file2 = strcat(Directory2,file2);
BandList = load(file2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtaining the loudest event per band and away from 0.25Hz lines

Nbands = length(BandList);

for bandnumber = 1:Nbands; %the  current frequency band
   fmin = BandList(bandnumber, 1);
   fmax = BandList(bandnumber, 2);

   %maximum doppler shift for this band
   dopp_wing(bandnumber)  = fmax*0.0001;
   dopp = dopp_wing(bandnumber);
% this is not very elegant but should work for now
   indices = find (FreqValues >= fmin & FreqValues <=fmax & ((mod(FreqValues,1) > dopp & mod(FreqValues,1) < 0.25-dopp) | (mod(FreqValues,1) > 0.25+dopp & mod(FreqValues,1) < 0.50-dopp) | (mod(FreqValues,1) > 0.50+dopp &mod(FreqValues,1) < 0.75-dopp) | (mod(FreqValues,1) > 0.75+dopp & mod(FreqValues,1) < 1.0-dopp) ));
   Nmax(bandnumber) = max( MaxValues(indices));
   bandnumber
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear FreqValues  MaxValues fmin fmax indices bandnumber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct filename "{Directory2}/{Detector}_coarse/MC_{Detector}_"
filepath = strcat(Directory2,subdir);
filepre = strcat(filepath,'/MC_');
filepre = strcat(filepre,Detector);
filepre = strcat(filepre,'_');

  nh0=10;

for bandnumber=0:Nbands-1
  bn=bandnumber+1;
  basestring = strcat(filepre, int2str( bandnumber ) );
  mystring = strcat(basestring, '_nc');
  Ncount0 = load(mystring);
  parstring = strcat(basestring, '_par');
  parvals = load(parstring);


  MC_FreqVals = parvals( 1:length(Ncount0(:,1)) , 2);


  clear parvals

  %use doppler wing calculated earlier--note that here the band index
  %starts from 0 and not from 1 as before
  dopp = dopp_wing(bandnumber + 1);
  vetoindices = find((mod(MC_FreqVals,1) > dopp & mod(MC_FreqVals,1) < 0.25-dopp) | (mod(MC_FreqVals,1) > 0.25+dopp & mod(MC_FreqVals,1) < 0.50+dopp) | (mod(MC_FreqVals,1) > 0.50+dopp &mod(MC_FreqVals,1) < 0.75-dopp) | (mod(MC_FreqVals,1) > 0.75+dopp & mod(MC_FreqVals,1) < 1.00-dopp));

  clear MC_FreqVals


  nMonteCarlos=length(vetoindices);



   for h0num=1:nh0
     x=Ncount0(vetoindices, h0num+1);
     kkcount = find(x>Nmax(bn));
     CH(h0num) = length(kkcount)/nMonteCarlos;
   end
   bandnumber
   CH

   fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n',bn-1, ...
           BandList(bn,1),BandList(bn, 2), Nmax(bn), CH);
end

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
